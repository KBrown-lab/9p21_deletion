## Data description:
1. Fasta files: 
Sequences were extracted from the human genome reference assembly hg38 (https://github.com/broadinstitute/gatk/blob/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz)
    - Reference file (Reference_hg38.fasta): 524,288 bp (chr9:21756689-22280977)
    - 9p21 deletion file (9p21_deletion_masked.fasta): same sequence coordinates as the reference file, but masked with Ns for the deletion (chr9:22180210-22280977) using bedtools maskfasta function

2. Bam file for ClinSV analysis: available upon request at dbGaP (phs004649.v1)
