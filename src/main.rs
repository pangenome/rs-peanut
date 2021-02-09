use std::fs::File;
use std::io::{self, BufRead, BufReader};

use bstr::ByteSlice;

use gfa::{
    cigar::CIGAR,
    gafpaf::parse_gaf,
    optfields::{OptFieldVal, OptFields, OptionalFields},
};

use clap::{App, Arg};

type GAF = gfa::gafpaf::GAF<OptionalFields>;

fn main() -> io::Result<()> {
    let arguments = App::new("peanut")
        .version("0.1.0")
        .author("Simon Heumos <simon.heumos@qbic.uni-tuebingen.de")
        .about("Evaluate GAF alignment quality")
        .arg(
            Arg::with_name("GAF")
                .short("g")
                .long("gaf")
                .required(true)
                .takes_value(true)
                .help("Input GAF file of which to evaluate the alignment quality."),
        )
        .get_matches();

    let gaf_filename = arguments.value_of("GAF").unwrap();
    let gaf_file_exists = std::path::Path::new(gaf_filename).exists();
    if !gaf_file_exists {
        eprintln!(
            "[peanut::main::error]: GAF file {} does not exist!",
            gaf_filename
        );
        std::process::exit(1);
    }

    let file = File::open(gaf_filename).unwrap();
    let mut lines = BufReader::new(file);
    let mut line: Vec<u8> = Vec::new();

    let mut cur_seq_len: usize = 0;
    let mut cur_seq_name: Vec<u8> = vec![0; 0];
    let mut first_line_seen: bool = false;
    let mut seq_name: Vec<u8>;
    let mut seq_len: usize = 0;
    let mut _map_start: usize = 0;

    let mut total_seq_len_qsc: usize = 0;
    let mut total_map_len_qsc: usize = 0;
    let mut nuc_bv_qsc: Vec<bool> = vec![false; 0];
    let mut nuc_bv_multi_qsc: Vec<bool> = vec![false; 0];

    loop {
        line.clear();
        let bytes_read = lines.read_until(0xA, &mut line);
        if bytes_read.is_err() || bytes_read.unwrap() == 0 {
            break;
        }
        let fields: bstr::Split = line[0..line.len()].split_str(b"\t");
        if let Some::<GAF>(gaf) = parse_gaf(fields) {
            let opt_fields = gaf.optional;
            let cigar = get_cigar(&opt_fields).unwrap();

            seq_name = gaf.seq_name;
            seq_len = gaf.seq_len;
            _map_start = gaf.seq_range.0;

            if !first_line_seen {
                first_line_seen = true;
                cur_seq_name = seq_name;
                cur_seq_len = seq_len;
                nuc_bv_qsc = vec![false; cur_seq_len];
                nuc_bv_multi_qsc = vec![false; cur_seq_len];
                eval_cigar(&cigar, &_map_start, &mut nuc_bv_qsc, &mut nuc_bv_multi_qsc);
            } else {
                if seq_name != cur_seq_name {
                    // finish the current one
                    total_seq_len_qsc += cur_seq_len;
                    total_map_len_qsc += nuc_bv_qsc.iter().filter(|&b| *b == true).count().clone();

                    nuc_bv_qsc = vec![false; seq_len];
                    nuc_bv_multi_qsc = vec![false; seq_len];

                    cur_seq_len = seq_len;
                    cur_seq_name = seq_name;
                    eval_cigar(&cigar, &_map_start, &mut nuc_bv_qsc, &mut nuc_bv_multi_qsc);
                } else {
                    eval_cigar(&cigar, &_map_start, &mut nuc_bv_qsc, &mut nuc_bv_multi_qsc);
                }
            }
        } else {
            eprintln!("Error parsing GAF line {}", line.as_bstr());
        }
    }

    // we have to add the last step
    total_seq_len_qsc += cur_seq_len;
    total_map_len_qsc += nuc_bv_qsc.iter().filter(|&b| *b == true).count().clone();

    let final_ratio_qsc: f64 = total_map_len_qsc as f64 / total_seq_len_qsc as f64;
    print!("{}\n", final_ratio_qsc);

    Ok(())
}

fn get_cigar<T: OptFields>(opts: &T) -> Option<CIGAR> {
    let cg = opts.get_field(b"cg")?;
    if let OptFieldVal::Z(cg) = &cg.value {
        CIGAR::from_bytestring(&cg)
    } else {
        None
    }
}

fn eval_cigar(
    cigar: &gfa::cigar::CIGAR,
    map_start: &usize,
    nuc_bv_qsc: &mut Vec<bool>,
    nuc_bv_multi_qsc: &mut Vec<bool>,
) {
    let cigar_iter = cigar.iter();
    let mut idx: usize = 0;
    for (len, op) in cigar_iter {
        // we have seen a match!
        use gfa::cigar::CIGAROp as Op;
        if matches!(op, Op::E) {
            // did we already mark this position?
            for offset in 0..len {
                nuc_bv_qsc[(idx + map_start + offset as usize)] = true;
            }
        }
        if matches!(op, Op::E | Op::M | Op::X | Op::I) {
            idx += len as usize;
        }
    }
}
