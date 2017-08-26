# Gene Set Enrichment Analysis GSEA.py

Python implementation of Gene Set Enrichment Analysis (GSEA) method as is
described in this [PNAS article](http://www.pnas.org/content/102/43/15545.abstract).

## Instalation

## Using as a standalone tool
To perform the GSEA with the standalone program `gsea.py` one needs to have two
files 
 - gene expression profiles
 - a list of gene sets

Example data can be found [here](https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999)) 
 
Example call of the program
```
python3 gsea.py leukemmia.txt pathways.txt
```
 
## Using as a Library

``` python
import gsea
result = gsea.gsea(profiles, gene_sets)
```
