# peanut

## GAF alignment evaluation tool.
_`peanut`_ calculates alignment metrics of a given [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) file from _[GraphAligner](https://github.com/maickrau/GraphAligner)_ evaluating the [CIGAR](https://metacpan.org/pod/Bio::Cigar#CIGAR-operations) string.
It outputs four [metrics](#metrics): 

1. [qsc](#query-sequence-containment-(qsc))
2. [uniq](#unique-query-sequence-matches-(uniq)) 
3. [multi](#multi-query-sequence-matches-(multi))
4. [nonaln](#non-query-sequence-matches-(nonaln))

## metrics

### query sequence containment (qsc)
<!--- 
https://jsfiddle.net/8ndx694g/
--->
<!---
\large qsc=\frac{#\!\!E}{query\_lens}
--->
<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20qsc%3D%5Cfrac%7B%23%5C!%5C!E%7D%7Bquery%5C_lens%7D">

- `#E` are the number of sequence matches (`=` or `E` symbol) in the GAF file. Nucleotide positions with sequence matches in multiple alignments are only counted once.
- `query_lens` is the length of all queries in the GAF in nucleotides.

### unique query sequence matches (uniq)
<!--- 
https://jsfiddle.net/8ndx694g/
--->
<!---
\large uniq=\frac{uniq\_\!\!#\!\!E}{query\_lens}
--->
<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20uniq%3D%5Cfrac%7Buniq%5C_%5C!%5C!%23%5C!%5C!E%7D%7Bquery%5C_lens%7D">

- `uniq_#E` are the number of unique sequence matches in the GAF file. 
- `query_lens` is the length of all queries in the GAF in nucleotides.

### multi query sequence matches (multi)
<!--
\large multi=\frac{multi\_\!\!#\!\!E}{query\_lens}
--->
<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20multi%3D%5Cfrac%7Bmulti%5C_%5C!%5C!%23%5C!%5C!E%7D%7Bquery%5C_lens%7D">

- `multi_#E` are the number of multiple sequence matches in the GAF file. Nucleotide positions with more than one multiple sequence matches are only counted once.
- `query_lens` is the length of all queries in the GAF in nucleotides.

### non query sequence matches (nonaln)
<!--
\large nonaln=\frac{nonaln\_\!\!#\!\!E}{query\_lens}
--->
<img src="https://render.githubusercontent.com/render/math?math=%5Clarge%20nonaln%3D%5Cfrac%7Bnonaln%5C_%5C!%5C!%23%5C!%5C!E%7D%7Bquery%5C_lens%7D">

- `non_#E` are the number of non-sequence matches in the GAF file.
- `query_lens` is the length of all queries in the GAF in nucleotides.

## usage
### building

```
git clone https://github.com/pangenome/rs-peanut.git
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
0.992910744238371	0.9926967987671109	0.00021394547126006352	0.007089255761628998
```
The first number is the `qsc`, the second number is the `uniq`, and the third number is the `multi`, and the fourth number is the `nonaln`.
## TODOs
- [x] Add query sequence alignment match mismatch (qsamm).
- [x] Describe `qsc`.
- [x] Remove non-helping metrics `qsamm` and `qsm`.
- [x] Add 3 new metrics: number of `uniq`ue query base alignments, number of `multi`ple query base alignments, and number of `nonaln` query bases.

## limits
So far, it has not been tested if _`peanut`_ also works with GAF files not originating from GraphAligner.
