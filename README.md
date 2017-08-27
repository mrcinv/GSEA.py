# Gene Set Enrichment Analysis GSEA.py [![Build Status](https://travis-ci.org/mrcinv/GSEA.py.svg?branch=master)](https://travis-ci.org/mrcinv/GSEA.py)

Python implementation of Gene Set Enrichment Analysis (GSEA) method as is
described in this [PNAS article](http://www.pnas.org/content/102/43/15545.abstract).

## Introduction to GSEA
The method is used for the interpretation of the genome wide mRNA expression
profiles for a collection of samples that can be classified into two separate classes.

The method assume, that genes are ordered in a list *L* according to
their differential expressions between the two classes of samples. To interpret
the list *L*, the method GSEA uses prior biological knowledge in the form of
different sets of genes *S* that are obtained from publicly published information
from previous experiments.

For a given set of genes *S*, the method try to determine whether the genes in
*S* are predominately at the top, the bottom or randomly distributed across the
list *L*.

### Input data

 - an *N x k* table *D* of expression data with *N* genes and *k* samples
 - ranking procedure to produce ranked gene list *L*
 - an exponent *p* to control the weight of the step
 - a gene set *S* (an independently obtained set of *N_H* genes)


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
