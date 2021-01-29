use std::fs::File;
use std::io::{self, BufReader};

use bstr::{io::*, ByteSlice};

use gfa::{
    cigar::{CIGAR},
    gafpaf::{parse_gaf},
    optfields::{OptFieldVal, OptFields, OptionalFields},
};

use clap::{App, Arg};

type GAF = gfa::gafpaf::GAF<OptionalFields>;

fn get_cigar<T: OptFields>(opts: &T) -> Option<CIGAR> {
    let cg = opts.get_field(b"cg")?;
    if let OptFieldVal::Z(cg) = &cg.value {
        CIGAR::from_bytestring(&cg)
    } else {
        None
    }
}

fn eval_cigar(cigar: & gfa::cigar::CIGAR, nuc_v: &mut Vec<bool>, map_start: & usize, nuc_o: &mut usize) {
    let cigar_iter = cigar.iter();
    let mut idx: usize = 0;
    for val in cigar_iter {
        // we have seen a match!
        use gfa::cigar::CIGAROp as Op;
        if matches!(val, Op::E) {
            // did we already mark this position?
            let nuc_b = nuc_v[(idx + map_start)];
            if nuc_b {
                *nuc_o += 1;
            } else {
                nuc_v[(idx + map_start)] = true;
            }
        }
        if matches!(val, Op::E | Op::M | Op::X | Op::I) {
            idx += 1;
        }
    }
}

fn main() -> io::Result<()> {
    let arguments = App::new("peanut")
        .version("0.1.0")
        .author("Simon Heumos <simon.heumos@qbic.uni-tuebingen.de")
        .about("Evaluate GAF alignment quality")
        .arg(Arg::with_name("GAF")
            .short("g")
            .long("gaf")
            .required(true)
            .takes_value(true)
            .help("Input GAF file of which to evaluate the alignment quality."))
        .get_matches();

    let gaf_filename = arguments.value_of("GAF").unwrap();
    let gaf_file_exists = std::path::Path::new(gaf_filename).exists();
    if !gaf_file_exists {
        eprintln!("[peanut::main::error]: GAF file {} does not exist!", gaf_filename);
        std::process::exit(1);
    }

    let file = File::open(gaf_filename).unwrap();
    let lines = BufReader::new(file).byte_lines().map(|l| l.unwrap());

    let mut cur_seq_len: usize = 0;
    let mut cur_seq_name: String = String::from("");
    let mut first_line_seen: bool = false;

    let mut total_seq_len: usize = 0;
    let mut total_map_len: usize = 0;

    let mut seq_name: String;
    let mut seq_len: usize = 0;
    let mut _map_start: usize = 0;

    let mut nuc_v: Vec<bool> = vec![false; 0];
    let mut nuc_o: usize = 0;

    for (i, line) in lines.enumerate() {
        let fields: bstr::Split = line.split_str(b"\t");
        if let Some::<GAF>(gaf) = parse_gaf(fields) {
            let opt_fields = gaf.optional;
            let cigar = get_cigar(&opt_fields).unwrap();
            
            // TODO Is there a faster way to do this? Seems super ugly.
            let chars: Vec<char> = gaf.seq_name.as_bytes().chars().collect();
            seq_name = chars.into_iter().collect();
            seq_len = gaf.seq_len;
            _map_start = gaf.seq_range.0;

            if !first_line_seen {
                first_line_seen = true;
                cur_seq_name = seq_name;
                cur_seq_len = seq_len;
                nuc_v = vec![false; cur_seq_len];
                eval_cigar(& cigar, &mut nuc_v, & _map_start, &mut nuc_o);
            } else {
                if seq_name != cur_seq_name {
                    // finish the current one
                    total_seq_len += cur_seq_len;
                    total_seq_len += nuc_o;
                    total_map_len += nuc_o;

                    total_map_len += nuc_v.iter().filter(|&b| *b == true).count().clone();

                    nuc_o = 0;
                    nuc_v = vec![false; seq_len];
                    cur_seq_len = seq_len;
                    cur_seq_name = seq_name;
                    eval_cigar(& cigar, &mut nuc_v, & _map_start, &mut nuc_o);
                } else {
                    eval_cigar(& cigar, &mut nuc_v, & _map_start, &mut nuc_o);
                }
            }
        } else {
            eprintln!("Error parsing GAF line {}", i);
        }
    }

    // we have to add the last step
    total_seq_len += seq_len;
    total_seq_len += nuc_o;
    total_map_len += nuc_o;
    total_map_len += nuc_v.iter().filter(|&b| *b == true).count();

    let final_ratio: f64 = total_map_len as f64 / total_seq_len as f64;
    println!("{}", final_ratio);

    Ok(())
}