This script sorts through fastq files generated from a GTseq library run and grabs reads that match the forward sequence that's hardcoded in the script. It then counts allele1 and allele2 (also hard-coded) and looks at the read ratio between them to determine the genotypes at the locus. 
input files: fastq files, text file with name of sample in column 1 and phenotypic sex in column 2. The text file should be in same directory as script. The output file will be in same directory as script.
hard-coded values in code: sequence of allele 1 ("GGTTTTTTTTC"), sequence of allele 2 ("GGTTTTTTTTTC"), sequence of forward primer ("GCTGTCAGATAAATGTAGTGAAACAAA"), name of file with phenotypic calls, name of output file
output file: text file with 10 columns: sample name, allele 1 counts, allele 2 counts, forward primer counts, sum of allele1 and allele 2, allele1/allele2, Genotype (ZZ,ZW,WW,TooLow,WeirdRatio), Genotypic Sex (ZZ=male, else female), Phenotypic sex, Concordance between genotypic and phenotypic sex
Total read counts <10 = no call ("Too Low")
Allele1/Allele2 counts >= 10, called as WW genotype
Allele1/Allele2 counts <= 0.2, called as ZW genotype
0.70 <= Allele1/Allele2 counts <= 3, called as ZW genotype
Any other value = "Weird Ratio"
It compares phenotypic sex calls (an input text file) with the genotypic sex calls 
