# dase

conditional-specific Differential Alternative Splicing variants Estimetion methd without reference genome sequence

## Description

DASE is a method for conditional-specific alternative splicing estimation based on de novo transcriptome assembly.

## Features

DASE does not require reference genome. Instead, DASE requires the output sequence file of Trinity and the expression files of each contig under all different conditions calculated by some tool, Kallisto for example.

## Requirement 

Python (>= 2.7)
MAFFT
NumPy
Scikit-learn

## Usage

    python dase.py -ex ex1.dat,ex2.dat -mafft /path/to/mafft -seq sequences.fasta -ef expression_format -th threshold_of_log(expression) -o output_file.dat -nc min_contigs -gap gap_penalty

## Author

Kouki Yonezawa (k_yonezawa@nagahama-i-bio.ac.jp)

## Reference of DASE

Kouki Yonezawa, Shunichi Shigeno, Tsukasa Mori, Atsushi Ogura. DASE: conditional-specific Differential Alternative Splicing variants Estimation method without reference genome sequence, and its application to non-model organisms. Proc. The IEEE International Conference on Bioinformatics and Biomedicine (BIBM2016), Shenzhen, 2016.

