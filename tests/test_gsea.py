# Dummy test
import numpy as np
from gsea import *
from numpy.testing import assert_almost_equal

def test_rank_genes():
    D = np.array([[-1,1],[1,-1]])
    C = [0,1]
    L,r = rank_genes(D,C)
    assert_almost_equal(L, [0,1])
    assert_almost_equal(r, [1,-1])

def test_enrichment_score():
    L = [1,0]
    r = [-1,1]
    S = [0,1]
    ES = enrichment_score(L,r,S)
    assert_almost_equal(ES,1)

    L = [0,1,2]
    r = [-1,0,1]
    assert_almost_equal(enrichment_score(L,r,[0]),1)
    assert_almost_equal(enrichment_score(L,r,[1]),-1)
