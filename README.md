# peanut

## GAF alignment evaluation tool.
_`peanut`_ calculates alignment metrics of a given [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) file from _[GraphAligner](https://github.com/maickrau/GraphAligner)_ evaluating the [CIGAR](https://metacpan.org/pod/Bio::Cigar#CIGAR-operations) string.
It outputs two [metrics](#metrics): 

1. [qsm](#query-sequence-match-(qsm)) 
2. [qsamm](#query-sequence-alignment-match-mismatch-(qsamm))

## metrics
### query sequence match (qsm)
<!--- 
https://jsfiddle.net/8ndx694g/
--->
<!---
\large \frac{(query_1\_\!\!#\!\!E + query_1\_multi\_\!\!#\!\!E) + \dots + (query_n\_\!\!#\!\!E + query_n\_multi\_\!\!#\!\!E)}{(query_1\_len + query_1\_multi\_\!\!#\!\!E) + \dots + (query_n\_len + query_n\_multi\_\!\!#\!\!E)}
--->
<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20%5Cfrac%7B(query_1%5C_%5C!%5C!%23%5C!%5C!E%20%2B%20query_1%5C_multi%5C_%5C!%5C!%23%5C!%5C!E)%20%2B%20%5Cdots%20%2B%20(query_n%5C_%5C!%5C!%23%5C!%5C!E%20%2B%20query_n%5C_multi%5C_%5C!%5C!%23%5C!%5C!E)%7D%7B(query_1%5C_len%20%2B%20query_1%5C_multi%5C_%5C!%5C!%23%5C!%5C!E)%20%2B%20%5Cdots%20%2B%20(query_n%5C_len%20%2B%20query_n%5C_multi%5C_%5C!%5C!%23%5C!%5C!E)%7D">

- `query_#E` are the number of sequence matches (`=` or `E` symbol) in the CIGAR of all GAF lines of a query. Nucleotide positions with sequence matches in multiple alignments are only counted once.
- `query_multi_#E` are the number of sequence matches of overlapping, multiple alignments of nucleotide positions.
- `query_len` is the length of the query in nucleotides.

### query sequence alignment match mismatch (qsamm)
<!--
\large \frac{(query_1\_\!\!#\!\!E\!\!#\!\!M\!\!#\!\!X + query_1\_multi\_\!\!#\!\!E\!\!#\!\!M\!\!#\!\!X) + \dots + (query_n\_\!\!#\!\!E\!\!#\!\!M\!\!#\!\!X + query_n\_multi\_\!\!#\!\!E\!\!#\!\!M\!\!#\!\!X)}{(query_1\_len + query_1\_multi\_\!\!#\!\!E\!\!#\!\!M\!\!#\!\!X) + \dots + (query_n\_len + query_n\_multi\_\!\!#\!\!E\!\!#\!\!M\!\!#\!\!X)}
--->
<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20%5Cfrac%7B(query_1%5C_%5C!%5C!%23%5C!%5C!E%5C!%5C!%23%5C!%5C!M%5C!%5C!%23%5C!%5C!X%20%2B%20query_1%5C_multi%5C_%5C!%5C!%23%5C!%5C!E%5C!%5C!%23%5C!%5C!M%5C!%5C!%23%5C!%5C!X)%20%2B%20%5Cdots%20%2B%20(query_n%5C_%5C!%5C!%23%5C!%5C!E%5C!%5C!%23%5C!%5C!M%5C!%5C!%23%5C!%5C!X%20%2B%20query_n%5C_multi%5C_%5C!%5C!%23%5C!%5C!E%5C!%5C!%23%5C!%5C!M%5C!%5C!%23%5C!%5C!X)%7D%7B(query_1%5C_len%20%2B%20query_1%5C_multi%5C_%5C!%5C!%23%5C!%5C!E%5C!%5C!%23%5C!%5C!M%5C!%5C!%23%5C!%5C!X)%20%2B%20%5Cdots%20%2B%20(query_n%5C_len%20%2B%20query_n%5C_multi%5C_%5C!%5C!%23%5C!%5C!E%5C!%5C!%23%5C!%5C!M%5C!%5C!%23%5C!%5C!X)%7D">

- `query_#E_#M_#X` are the number of sequence matches (`=` or `E` symbol), the number of alignment matches (`M` symbol), and the number of sequence mismatches (`X` symbol) in the CIGAR of all GAF lines of a query. Nucleotide positions with sequence matches, alignment matches, or sequence mismatches in multiple alignments are only counted once.
- `query_multi_#E_#M_#X` are the number of sequence matches, the number of alignment matches, and the number of sequence mismatches of overlapping, multiple alignments of nucleotide positions.
- `query_len` is the length of the query in nucleotides.

## usage
### building

```
git clone https://github.com/subwaystation/rs-peanut.git
cd rs-peanut
cargo build --release
```

### example
_`peanut`_ requires as an input a GAF file `-g`.
```
./target/release/peanut -g aln.gaf
```
The output is written to stdout in a tab-delimited format.
```
0.9927656450276804	0.9968190451624734
```
The first number is the `qsm`, the second number is the `qsamm`.
## TODOs
- [x] query sequence alignment match mismatch (qsamm)
- [ ] describe `qsc`

## limits
So far, it has not been tested if _`peanut`_ also works with GAF files not originating from GraphAligner.
