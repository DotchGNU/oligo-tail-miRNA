# oligo-tail-miRNA
R scripts used to analyze miRNA isoforms trimming and tailing (as well as nucleotide composition of non-templated tails).

### Minimum number of "X" nt in tail

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
```
