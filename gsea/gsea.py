# MIT License

# Copyright (c) 2017 Martin Vuk

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import numpy as np

def enrichment_score(L, r, S, p_exp=1):
    """Calculates enrichment score (ES) for a given gene expression data.

    Arguments:
    ---------
    L: an ordered list of indexes for N genes

    r: a list of correlations of a gene expression data with phenotypes for N genes

    C: a list of classes(phenotypes) for k samples

    S: a list of gene indexes that belong to a gene set

    p_exp: exponent parameter to control the weight of the step (defaults to 1)

    Returns:
    -------

    ES: enrichment score for the given gene set S
    """

    N = len(L)
    S_mask = np.zeros(N)
    S_mask[S] = 1
    # reorder gene set mask
    S_mask = S_mask[L]
    N_R = sum(abs(r*S_mask)**p_exp)
    P_hit = np.cumsum(abs(r*S_mask)**p_exp)/N_R if N_R!=0 else np.zeros_like(S_mask)
    N_H = len(S)
    P_mis = np.cumsum((1-S_mask))/(N-N_H) if N!=N_H else np.zeros_like(S_mask)
    idx = np.argmax(abs(P_hit - P_mis))
    return P_hit[idx] - P_mis[idx]

def rank_genes(D,C):
    """Ranks genes in expression dataset according to the correlation with a
    phenotype.

    Arguments:
    ----------
    D: a 2D array of expression data for N genes and k samples

    C: a list of values 1, if a i-th sample is in the phenotype or 0 otherwise

    Returns:
    --------
    L: an ordered list of gene indexes

    r: a ordered list of correlation coefficients

    """
    N, k = D.shape
    rL = []
    for i in range(N):
        rL.append(np.corrcoef(D[i,:],C)[0,1])
    rL = sorted(enumerate(rL), key=lambda x: x[1])
    r = [x[1] for x in rL]
    L = [x[0] for x in rL]
    return L, r

# Multiple Hypothesis testing

def gsea(D, C, S_sets, p_exp=1, random_sets=1000):
    """Performs Multiple Hypotesis Testing.

    Arguments:
    ----------
    D: a 2D array of expression data for N genes and k sample

    C: a list of values 1, if a i-th sample is in the phenotype or 0 otherwise

    S_sets: a list of variable length of gene indexes that belong to a gene set

    p_exp: exponent parameter to control the weight of the step (defaults to 1)

    random_sets: number of randomly generated gene sets


    Returns:
    --------

    order: list of gene indexes, ordered by the greatest Normalized Enrichment Scores

    NES: a list of Normalized Enrichment Scores (ordered by absolute value)

    p_value: a list of p-values for the NES
    """
    N, k = D.shape
    l = len(S_sets)
    # generate random gene sets
    Pi_sets = []
    p_value = np.zeros(l)
    NES = np.zeros(l)
    ES = np.zeros(l)
    ES_pi = np.zeros((random_sets,l))
    L, r = rank_genes(D, C)
    # enrichment scores for S_i
    for i in range(l):
        ES[i] = enrichment_score(L,r,S_sets[i],p_exp)
    for i in range(random_sets):
        pi = np.array([np.random.randint(0,2) for i in range(k)])
        L, r = rank_genes(D,pi)
        ES_pi[i,:] = [enrichment_score(L,r,S_sets[j],p_exp) for j in range(l)]

    # calculate normalized enrichment scores and p-values
    for i in range(l):
        # normalize separately positive and negative values
        ES_plus = ES_pi[i,:][ES_pi[i,:]>0]
        ES_minus = ES_pi[i,:][ES_pi[i,:]<0]
        mean_plus = np.mean(ES_plus)
        mean_minus = np.mean(ES_minus)
        if ES[i]>0:
            NES[i] = ES[i]/mean_plus
            p_value[i] = sum(ES_plus>ES[i])/len(ES_plus)
        elif ES[i]<0:
            NES[i] = -ES[i]/mean_minus
            p_value[i] = sum(ES_minus<ES[i])/len(ES_minus)
    NES_sort = sorted(enumerate(NES),key=lambda x: -abs(x[1]))
    order = [x[0] for x in NES_sort]
    NES = [x[1] for x in NES_sort]
    return order,NES,p_value[order]

