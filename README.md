# miRNA trimming and oligo-tailing analysis
Description of the methods and R scripts used to analyze miRNA isoforms trimming and tailing (as well as nucleotide composition of non-templated tails).

## **Analysis of miRNA sequencing data and tail composition:**

The small RNA sequencing data were analyzed using an in-house pipeline. Briefly, adaptors were removed, reads were mapped using Bowtie and visualized using IGV. More detailed study of the isomiR profile was done using [QuagmiR](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA#references). This software uses a unique algorithm to pull specific reads and aligns them against a consensus sequence in the middle of a miRNA, allowing mismatches on the ends to capture 3’ isomiRs. The reports included tabulated analysis of miRNA expression, length, number of nucleotides trimmed and tail composition at individual read level. Customized R scripts were used to calculate percentages of canonical miRNA (defined as the most abundant templated read) and 3’ isomiRs, a well as percentages of tailing and trimming. Long tail composition was calculated by counting the number of non-templated nucleotides present in the tail of each isomiR read. Reads with equal number of non-templated nucleotides in the tail were added together and cumulative distribution was calculated for all the oligo-tailed isomiRs going from ones with longer to shorter tails.

  * [Long tail composition](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA/tree/master/Long%20Tail%20Composition)
  * [Descriptive example of the analysis performed](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA#minimum-number-of-n-nucleotide-in-tail)

## **Analysis of isomiR profiles on AGO1 and AGO2 from TCGA:**

Tumoral samples from TCGA bearing genomic mutations in either AGO1 or AGO2 leading to missense and synonymous amino acid changes were identified from Genomic Data Commons Data Portal (https://portal.gdc.cancer.gov/, accessed during May 2019). GDC uses combined reports from several variant callers (mutect2, varscan, muse and somaticsniper).

Selected Case ID were: P295L TCGA-53-A4EZ, R315M TCGA-HU-A4G8 and E299K TCGA-Z6-A8JE (AGO2), F310L TCGA-94-7033 (AGO1). The analysis of selected patient samples was also performed using QuagmiR56, with a previous conversion of the bam files to fastq files by Picard Sam-to-Fastq, using Amazon cloud instances through the Seven Bridges Genomics implementation of the NCI Cancer Genomics Cloud. Mutations were plotted into the PDB structures of AGO1 and AGO2 using pymol.

  * AGO2
    * [GDC Portal AGO2 missense mutations](https://portal.gdc.cancer.gov/exploration?facetTab=mutations&filters=%7B%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22ENSG00000123908%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22ssms.consequence.transcript.consequence_type%22%2C%22value%22%3A%5B%22missense_variant%22%5D%7D%7D%5D%2C%22op%22%3A%22and%22%7D&searchTableTab=mutations)
    * [AGO2 scripts](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA/tree/master/AGO2%20mutants)
    
  * AGO1
    * [GDC Portal AGO1 missense mutations](https://portal.gdc.cancer.gov/exploration?facetTab=mutations&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22ENSG00000092847%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22ssms.consequence.transcript.consequence_type%22%2C%22value%22%3A%5B%22missense_variant%22%5D%7D%7D%5D%7D&searchTableTab=mutations)
    * [AGO1 scripts](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA/tree/master/AGO1%20mutants)



## **Minimum number of "N" nucleotide in tail**

**miRBase reference**
```
>hsa-miR-7-5p MIMAT0000252
UGGAAGACUAGUGAUUUUGUUGUU
```

**Example long tailed read**
```
<--templated-----------><--non-templated-->
UGGAAGACUAGUGAUUUUGUUGUUUUUUUUUAAUUUUGUCUUU
........................UUUUUUUAAUUUUGUCUUU

Number of U in tail: 15
Number of A in tail: 2
Number of G in tail: 1
Number of C in tail: 1
```

## **References:**
* [QuagmiR: A Cloud-based Application for IsomiR Big Data Analytics.](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty843/5123434)
Bofill-De Ros X, Chen K, Chen S, Tesic N, Randjelovic D, Skundric N, Nesic S, Varjacic V, Williams EH, Malhotra R, Jiang M, Gu S. Bioinformatics. 2018 Oct 8. doi: 10.1093/bioinformatics/bty843.([Pubmed link](https://www.ncbi.nlm.nih.gov/pubmed/30295744))
