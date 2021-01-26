use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;


fn main() {
    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines("/home/heumos/Downloads/pgge/aln.gaf") {
        // Consumes the iterator, returns an (Optional) String
        let mut i: u64 = 0;
        let mut cur_seq_len: u64 = 0;
        let mut cur_seq_name: &str = "";
        let mut first_line_seen: bool = false;
        for line in lines {
            if let Ok(ip) = line {
                // println!("{}", ip);
                let gaf: String = String::from(ip);
                let v: Vec<&str> = gaf.split("\t").collect();
                //println!("{}", v[0..1]);
                // println!("{:?}", &v[..]);
                i = i + 1;
                //println!("{:?}", i);
                let seq_name: &str = v[0];
                let seq_len: u64 = v[1].parse().unwrap();
                let map_start: u64 = v[2].parse().unwrap();
                let map_end: u64 = v[3].parse().unwrap();
                // println!("{}", v[15]);
                let map_id_v: Vec<&str> = v[15].split(":").collect();
                let map_id: f64 = map_id_v[2].parse().unwrap();
                println!("{}\t{}\t{}\t{}\t{}", seq_name, seq_len, map_start, map_end, map_id);
                if !first_line_seen {
                    first_line_seen = true;
                    // FIXME
                    // cur_seq_name = seq_name;
                    cur_seq_len = seq_len;
                } else {
                    // do we have a change of chromosome?
                    // FIXME the following creates an out of scope error
                    if seq_name != cur_seq_name {
                        println!("OH YEAH");
                    } 
                    // FIXME
                    // cur_seq_name = seq_name;
                    // cur_seq_len = seq_len;
                }
            }
        }
    }
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}