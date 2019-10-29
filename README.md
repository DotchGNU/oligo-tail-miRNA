# miRNA trimming and oligo-tailing analysis
R scripts used to analyze miRNA isoforms trimming and tailing (as well as nucleotide composition of non-templated tails).

### Minimum number of "N" nucleotide in tail

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

### References:
* [QuagmiR: A Cloud-based Application for IsomiR Big Data Analytics.](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty843/5123434)
Bofill-De Ros X, Chen K, Chen S, Tesic N, Randjelovic D, Skundric N, Nesic S, Varjacic V, Williams EH, Malhotra R, Jiang M, Gu S. Bioinformatics. 2018 Oct 8. doi: 10.1093/bioinformatics/bty843.([Pubmed link](https://www.ncbi.nlm.nih.gov/pubmed/30295744))
