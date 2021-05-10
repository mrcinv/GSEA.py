# Gene Set Enrichment Analysis GSEA.py [![Build Status](https://travis-ci.org/mrcinv/GSEA.py.svg?branch=master)](https://travis-ci.org/mrcinv/GSEA.py)

Python implementation of Gene Set Enrichment Analysis (GSEA) method as
described by A. Subramanian et al. ([PNAS article](http://www.pnas.org/content/102/43/15545.abstract)).


## Installation

Installation is best done using *pip*

```
pip install https://github.com/mrcinv/GSEA.py/archive/master.zip
```

Or clone the repo and run `setup.py`

```
git clone https://github.com/mrcinv/GSEA.py.git
cd GSEA.py
python3 setup.py install
```

## Using GSEA.py as a standalone tool
To perform the GSEA with the standalone program `gsea.py` two
files have to be provided
 - gene expression profiles
 - a list of gene sets.

 
Example call of the program
```
gsea leukemmia.txt pathways.txt
```

should produce results similar to this:

```
gene set	NES	p-value
MAP00480_Glutathione_metabolism  -1.84   0.0040404
DNA_DAMAGE_SIGNALLING    1.83    0.0020080
MAP00193_ATP_synthesis   -1.82   0.0020040
MAP00195_Photosynthesis  -1.82   0.0020040
MAP03070_Type_III_secretion_system       -1.82   0.0020040
Cell_Cycle       1.74    0.0294695
ANTI_CD44_UP     -1.74   0.0038835
ST_Ga13_Pathway  -1.73   0.0100000
CR_CELL_CYCLE    1.72    0.0187891
...
```

The results can also be piped into a tsv file and viewed in a spreadsheet program:

```
gsea leukemmia.txt pathways.txt > es.tsv
```

## Using GSEA.py as a Library

``` python
import gsea
order, NES, p_value = gsea.gsea(D, C, S_sets)
```

`gsea.gsea` also supports parallelization with Python `multiprocessing`. To use multiple cores, pass the `n_jobs` keyword argument.
To use all available cores, pass `n_jobs = -1`.

See [example Jupyter notebook](Leukemia.ipynb) for an ilustration.

## Sample data

The sample data that has been used in the original article 
can be downloaded [from the Broad Institute download page](http://software.broadinstitute.org/gsea/downloads.jsp).

Example expression profiles can also be found [here](https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999)). 

## Introduction to GSEA
The method is used for the interpretation of the genome-wide mRNA expression
profiles for a collection of samples from two different classes (phenotypes).

The method GSEA uses prior biological knowledge in the form of
different sets of genes *S* that are obtained from publicly available information
from previous experiments.

For a given set of genes *S*, the method calculates *Enrichment Score ES(S)*,
which shows if the genes from the gene set *S* are at the top or bottom in the
ordered list of genes *L* or randomly distributed throughout *L*.
The genes in *L*
are sorted according to the correlation of their expression values with phenotypes.


### Input data

 - an *N x k* table *D* of expression data with *N* genes and *k* samples
 - samples are grouped into two classes
 - a collection of gene sets *S* (an independently obtained sets of genes of different sizes)

### Output 

 - Normalized Gene Enrichment Scores for all input gene sets
 - P-values for all input gene sets
