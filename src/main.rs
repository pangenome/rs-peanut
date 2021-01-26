use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;


fn main() {
    // File hosts must exist in current path before this produces output
    if let Ok(lines) = read_lines("/home/heumos/Downloads/pgge/aln.gaf") {
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
                        println!("cur_seq_len: {}", cur_seq_len);
                    } 
                    let map_len: f64 = (map_end - map_start) * map_id;
                    total_map_len += map_len;
                    println!("map_len: {}", map_len);
                    cur_seq_name = seq_name.to_owned();
                    cur_seq_len = seq_len;
                }
            }
        }

        // we have to add the last step
        let map_len: f64 = (map_end - map_start) * map_id;
        total_map_len += map_len;
        total_seq_len += cur_seq_len;

        println!("total_map_len: {}", total_map_len);
        println!("total_seq_len: {}", total_seq_len);

        let final_ratio: f64 = total_map_len / total_seq_len as f64;
        println!("{}", final_ratio);
    }
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}