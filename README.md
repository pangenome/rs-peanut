# peanut

GAF alignment evaluation tool.

## metric
### query sequence match (qsm)
<!--- 
https://jsfiddle.net/8ndx694g/
--->
<!---
\frac{(query_1\_matches + query_1\_multi\_matches) + (query_2\_matches + query_2\_multi\_matches) + \dots + (query_n\_matches + query_n\_multi\_matches)}{(query_1\_len + query_1\_multi\_matches) + (query_2\_len + query_2\_multi\_matches) + \dots + (query_n\_len + query_n\_multi\_matches)}
--->
<img src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B(query_1%5C_matches%20%2B%20query_1%5C_multi%5C_matches)%20%2B%20(query_2%5C_matches%20%2B%20query_2%5C_multi%5C_matches)%20%2B%20%5Cdots%20%2B%20(query_n%5C_matches%20%2B%20query_n%5C_multi%5C_matches)%7D%7B(query_1%5C_len%20%2B%20query_1%5C_multi%5C_matches)%20%2B%20(query_2%5C_len%20%2B%20query_2%5C_multi%5C_matches)%20%2B%20%5Cdots%20%2B%20(query_n%5C_len%20%2B%20query_n%5C_multi%5C_matches)%7D">

- `query_matches` are the number of matches (`=` symbol) in the CIGAR of a GAF line. Nucleotide positions with matches in multiple alignments are only counted once.
- `query_multi_matches` are the number of matches of multiple alignments of nucleotide positions.
- `query_len` is the length of the query in nucleotides.

## usage
### building

`git clone https://github.com/subwaystation/rs-peanut.git`

`cd rs-peanut`

`cargo build --release`

### example

`./target/release/peanut --gaf aln.gaf`

## TODOs
- [ ] query sequence alignment match mismatch (qsamm)