# Gene Set Enrichment Analysis GSEA.py [![Build Status](https://travis-ci.org/mrcinv/GSEA.py.svg?branch=master)](https://travis-ci.org/mrcinv/GSEA.py)

Python implementation of Gene Set Enrichment Analysis (GSEA) method as is
described in this [PNAS article](http://www.pnas.org/content/102/43/15545.abstract).


 
## Instalation

Instalation is eacy using *pip*

```
pip install https://github.com/mrcinv/GSEA.py/archive/master.zip
```

Or clone the repo and run `setup.py`

```
git clone https://github.com/mrcinv/GSEA.py.git
cd GSEA.py
python3 setup.py install
```

## Using as a standalone tool
To perform the GSEA with the standalone program `gsea.py` one needs to have two
files 
 - gene expression profiles
 - a list of gene sets

 
Example call of the program
```
python3 gsea.py leukemmia.txt pathways.txt
```

should produce the following results:

```
Gene set, normalized enrichment score (NES), p-value

```
The results can also be piped into a csv file and viewed in a spreadsheet program:

```
python3 gsea.py leukemmia.txt pathways.txt > es.csv
```

## Using as a Library

``` python
import gsea
result = gsea.gsea(profiles, gene_sets)
```

## Sample data

You can download the sample data that has been used in the original article 
[at the Broad Institute download page](http://software.broadinstitute.org/gsea/downloads.jsp).

Example expression profiles can also be found [here](https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999)) 

## Introduction to GSEA
The method is used for the interpretation of the genome wide mRNA expression
profiles for a collection of samples that can be classified into two separate classes.

The genes are ordered in a list *L* according to their correlation of expressions for two classes. 
To interpret the list *L*, the method GSEA uses prior biological knowledge in the form of
different sets of genes *S* that are obtained from publicly published information
from previous experiments.

For a given set of genes *S*, the method calculates *Enrichment Score (ES)*, which tells if the gene from the gene set *S* are at the top, bottom or randomly distributed in the ordered list *L*.

### Input data

 - an *N x k* table *D* of expression data with *N* genes and *k* samples
 - samples are grouped into two classes
 - an exponent *p* to control the weight of the step
 - a collection of gene sets *S* (an independently obtained sets of genes of different sizes)

### Output 

 - Normalized Gene Enrichment Scores for all gene sets
 - P-values for gene sets
