#!/bin/bash

blastn -query Table/02_16Smerge_clus_seq.fasta -subject Table/STD_prok.fasta -outfmt 6 -perc_identity 0.99 > Table/03_blastResult.txt

