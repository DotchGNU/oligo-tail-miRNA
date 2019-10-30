# miRNA trimming and oligo-tailing analysis
Description of the methods and R scripts used to analyze miRNA isoforms trimming and tailing (as well as nucleotide composition of non-templated tails).

## **Dataset availability**
[GEO accession GSE139567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139567)
* Currently only available under reviewer's token (10/29/2019)

## **System requirements**
This code was tested under:
* MacBook Pro (15-inch, 2016)
* Processor: 2.7 GHz Intel Core i7
* Memory: 16 GB 2133 MHz LPDDR3

### **R and RStudio**
* R version 3.5.1 (2018-07-02)
  * [Download and installation](https://www.r-project.org/)
* RStudio version 1.1.456 – © 2009-2018
  * [Download and installation](https://rstudio.com/)

### **Cloud computing tools**
All the cloud computing tools can be found in the [Cancer Genomics Cloud (CGC)](www.cancergenomicscloud.org).

The Cancer Genomics Cloud (CGC), powered by Seven Bridges, is one of three systems funded by the National Cancer Institute to explore the paradigm of colocalizing massive public datasets, like The Cancer Genomics Atlas (TCGA), alongside secure and scalable computational resources to analyze them.

* [QuagmiR](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA#references)
* [Picard Sam-to-Fastq](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.3.0/picard_sam_SamToFastq.php)

“The Seven Bridges Cancer Genomics Cloud has been funded in whole or in part with Federal funds from the National Cancer Institute, National Institutes of Health, Contract No. HHSN261201400008C and ID/IQ Agreement No. 17X146 under Contract No. HHSN261201500003I.”

## **Analysis of miRNA sequencing data and tail composition:**

The small RNA sequencing data were analyzed using an in-house pipeline. Briefly, adaptors were removed, reads were mapped using Bowtie and visualized using IGV. More detailed study of the isomiR profile was done using [QuagmiR](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA#references). This software uses a unique algorithm to pull specific reads and aligns them against a consensus sequence in the middle of a miRNA, allowing mismatches on the ends to capture 3’ isomiRs. The reports included tabulated analysis of miRNA expression, length, number of nucleotides trimmed and tail composition at individual read level. Customized R scripts were used to calculate percentages of canonical miRNA (defined as the most abundant templated read) and 3’ isomiRs, a well as percentages of tailing and trimming. Long tail composition was calculated by counting the number of non-templated nucleotides present in the tail of each isomiR read. Reads with equal number of non-templated nucleotides in the tail were added together and cumulative distribution was calculated for all the oligo-tailed isomiRs going from ones with longer to shorter tails.

  * [Long tail composition scripts](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA/tree/master/Long%20Tail%20Composition)
  * [Descriptive examples of the analysis performed](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA#descriptive-example-of-the-analysis-performed)

## **Analysis of isomiR profiles on AGO1 and AGO2 from TCGA:**

Tumoral samples from TCGA bearing genomic mutations in either AGO1 or AGO2 leading to missense and synonymous amino acid changes were identified from Genomic Data Commons Data Portal (https://portal.gdc.cancer.gov/, accessed during May 2019). GDC uses combined reports from several variant callers (mutect2, varscan, muse and somaticsniper).

Selected Case ID were: P295L TCGA-53-A4EZ, R315M TCGA-HU-A4G8 and E299K TCGA-Z6-A8JE (AGO2), F310L TCGA-94-7033 (AGO1). The analysis of selected patient samples was also performed using QuagmiR56, with a previous conversion of the bam files to fastq files by Picard Sam-to-Fastq, using Amazon cloud instances through the Seven Bridges Genomics implementation of the NCI Cancer Genomics Cloud. Mutations were plotted into the PDB structures of AGO1 and AGO2 using pymol.

  * AGO2
    * [GDC Portal AGO2 missense mutations](https://portal.gdc.cancer.gov/exploration?facetTab=mutations&filters=%7B%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22ENSG00000123908%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22ssms.consequence.transcript.consequence_type%22%2C%22value%22%3A%5B%22missense_variant%22%5D%7D%7D%5D%2C%22op%22%3A%22and%22%7D&searchTableTab=mutations)
    * [AGO2 scripts](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA/tree/master/AGO2%20mutants)
    
  * AGO1
    * [GDC Portal AGO1 missense mutations](https://portal.gdc.cancer.gov/exploration?facetTab=mutations&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22genes.gene_id%22%2C%22value%22%3A%5B%22ENSG00000092847%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22ssms.consequence.transcript.consequence_type%22%2C%22value%22%3A%5B%22missense_variant%22%5D%7D%7D%5D%7D&searchTableTab=mutations)
    * [AGO1 scripts](https://github.com/Gu-Lab-RBL-NCI/oligo-tail-miRNA/tree/master/AGO1%20mutants)



## **Descriptive example of the analysis performed**

*The examples shown here are just to illustrate the logic implemented in the analysis and calculations used in the R scripts.*

**miRBase reference**
```
>hsa-miR-7-5p MIMAT0000252 (mature miRNA)
UGGAAGACUAGUGAUUUUGUUGUU

>hsa-mir-7-1 MI0000263 (pri-miRNA paralog 1)
<--mature-miRNA--------><---------templated (genomic reference)---------------------->
UGGAAGACUAGUGAUUUUGUUGUUUUUAGAUAACUAAAUCGACAACAAAUCACAGUCUGCCAUAUGGCACAGGCCAUGCCUCUACAG

>hsa-mir-7-2 MI0000264 (pri-miRNA paralog 2)
<--mature-miRNA--------><---------templated (genomic reference)--------------->
UGGAAGACUAGUGAUUUUGUUGUUGUCUUACUGCGCUCAACAACAAAUCCCAGUCUACCUAAUGGUGCCAGCCAUCGCA

>hsa-mir-7-3 MI0000265 (pri-miRNA paralog 3)
<--mature-miRNA--------><---------templated (genomic reference)---------------->
UGGAAGACUAGUGAUUUUGUUGUUCUGAUGUACUACGACAACAAGUCACAGCCGGCCUCAUAGCGCAGACUCCCUUCGAC
```


**Minimum number of "N" nucleotide in tail**
```
Example long tailed read:
<--templated-----------><--non-templated-->
UGGAAGACUAGUGAUUUUGUUGUUUUUUUUUAAUUUUGUCUUU
........................UUUUUUUAAUUUUGUCUUU

Number of U in tail: 15
Number of A in tail: 2
Number of G in tail: 1
Number of C in tail: 1
```

**Weighted Average of the Minimum number of U in oligo-tail**
```
Example long tailed reads:
<--templated-----------><--non-templated-->  U_in_tail  Counts  Fraction  Weighted_U_in_tail
UGGAAGACUAGUGAUUUUGUUGUU                     
........................UUUUUUUAAUUUUGUCUUU  15         100     0.2       3
........................UUUAUUU              6          100     0.2       1.2
........................UUUUUUU              7          100     0.2       1.4
........................UUU                  3          100     0.2       0.6
........................UU                   2          100     0.2       0.4

Weighted Average of the Minimum number of U in oligo-tail: 6.6
```

## **References:**
* [QuagmiR: A Cloud-based Application for IsomiR Big Data Analytics.](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty843/5123434)Bofill-De Ros X, Chen K, Chen S, Tesic N, Randjelovic D, Skundric N, Nesic S, Varjacic V, Williams EH, Malhotra R, Jiang M, Gu S. Bioinformatics. 2018 Oct 8. doi: 10.1093/bioinformatics/bty843.([Pubmed link](https://www.ncbi.nlm.nih.gov/pubmed/30295744))
* [The Cancer Genomics Cloud: Collaborative, Reproducible, and Democratized—A New Paradigm in Large-Scale Computational Research.](https://cancerres.aacrjournals.org/content/77/21/e3.long) Lau JW, Lehnert E, Sethi A, Malhotra R, Kaushik G, Onder Z, Groves-Kirkby N, Mihajlovic A, DiGiovanna J, Srdic M, Bajcic D, Radenkovic J, Mladenovic V, Krstanovic D, Arsenijevic V, Klisic D, Mitrovic M, Bogicevic I, Kural D, Davis-Dusenbery B; Seven Bridges CGC Team. The Cancer Genomics Cloud:
Collaborative, Reproducible, and Democratized-A New Paradigm in Large-Scale
Computational Research. Cancer Res. 2017 Nov 1;77(21)([Pubmed link](https://www.ncbi.nlm.nih.gov/pubmed/29092927))
