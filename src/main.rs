use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

extern crate clap;
use clap::{App, Arg};

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

    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines(gaf_filename) {
        let mut cur_seq_len: u64 = 0;
        let mut cur_seq_name: String = String::from("");
        let mut first_line_seen: bool = false;

        let mut total_seq_len: u64 = 0;
        let mut total_map_len: f64 = 0.0;

        let mut seq_name: &str;
        let mut seq_len: u64;
        let mut map_start: f64 = 0.0;
        let mut map_end: f64 = 0.0;
        let mut map_id: f64 = 0.0;

        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(ip) = line {
                let gaf: String = String::from(ip);
                let v: Vec<&str> = gaf.split("\t").collect();
                seq_name = v[0];
                seq_len = v[1].parse().unwrap();
                map_start = v[2].parse().unwrap();
                map_end = v[3].parse().unwrap();
                let map_id_v: Vec<&str> = v[15].split(":").collect();
                map_id = map_id_v[2].parse().unwrap();
                // println!("{}\t{}\t{}\t{}\t{}", seq_name, seq_len, map_start, map_end, map_id);
                if !first_line_seen {
                    first_line_seen = true;
                    cur_seq_name = seq_name.to_owned();
                    cur_seq_len = seq_len;
                } else {
                    // do we have a change of chromosome?
                    if seq_name != cur_seq_name {
                        total_seq_len += cur_seq_len;
                        // println!("cur_seq_len: {}", cur_seq_len);
                    } 
                    let map_len: f64 = (map_end - map_start) * map_id;
                    total_map_len += map_len;
                    // println!("map_len: {}", map_len);
                    cur_seq_name = seq_name.to_owned();
                    cur_seq_len = seq_len;
                }
            }
        }

        // we have to add the last step
        let map_len: f64 = (map_end - map_start) * map_id;
        total_map_len += map_len;
        total_seq_len += cur_seq_len;

        // println!("total_map_len: {}", total_map_len);
        // println!("total_seq_len: {}", total_seq_len);

        let final_ratio: f64 = total_map_len / total_seq_len as f64;
        println!("{}", final_ratio);
    }

    Ok(())
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}