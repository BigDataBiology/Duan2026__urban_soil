# ARGs table
Predicted ORF by fARGene, the gene class and sequence of the ORF; the contig and sample where they were found; and the gene name from ResFinder if there were any matches. 

## Columns:

- orf: header of the predicted ORF by fARGene in the format `sample@@@class@@@contig\_suffix. A suffix is added to avoid redundancy if more than one ARG from the same class is found in a single contig.
- centroid: the centroid to where the ORF belongs after clustering all ORFs at 95% identity at amino acid level.
- sample: the sample to where the ORF belongs to.
- contig: contig within the sample to where the ORF belongs to.
- class: *aac(2')*, *aac(3)*, *aac(6')*, *aph(2'')*, *aph(3')*, *aph(6)*, beta-lactamase A, beta-lactamase B1-B2, beta-lactamase B3, beta-lactamase C, beta-lactamase D, *erm*, *mph*, *qnr*, tet efflux (tetracycline efflux pump), tet enzyme (tetracycline inactivation enzyme), tet RPG (tetracycline Ribosomal Protection Gene).
- hclass: aminoglycoside, beta-lactamase, macrolide, quinolone, tetracycline.
- sequence: amino acid sequence of the ORF.
- resfinder: the gene name from ResFinder if the sequence in nucleotide format is detected by ResFinder.
- start: start of the alignment of the ORF in nucleotide format in its corresponding contig.
- end: end of the alignment of the ORF in nucleotide format in its corresponding contig.




