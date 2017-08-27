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

def enrichment_score(D, C, S, p_exp=1):
    """Calculates enrichment score (ES) for a given gene expression data.

    Arguments:
    ---------
    D: a 2D array with expression data for N genes and k samples

    C: a list of classes(phenotypes) for k samples

    S: a list of gene indexes

    p_exp: exponent parameter to control the weight of the step (defaults to 1)

    Returns:
    -------

    ES: enrichment score
    """

    N, k = D.shape
    L, r = rank_genes(D,C)
    S_mask = np.zeros(N)
    S_mask[S] = 1
    N_R = sum(abs(r*S_mask)**p)
    P_hit = cumsum(abs(r*S_mask)**p)/N_R
    N_H = len(S)
    P_mis = cumsum((1-S-mask))/(N-N_H)
    ES = np.amax(P_hit - P_mis)
    return ES

def rank_genes(D,C):
    """Ranks genes in expression dataset according to the correlation with the
    profile of interest."""
    return None

def generate_ES_null(D, S, p_exp=1, samples=1000):
    """Generates distribution of enrichment scores (ES) for randomly assigned
    phenotipes.
    """
    N, k = D.shape
    ES_null = []
    for i in range(samples):
        C = [np.random.randint(0,2) for i in range(k)]
        ES_null.append(enrichment_score(D,C,S,p_exp))
    ES_null.sort()
    return ES_null

def estimate_P(ES, ES_null):
    """Estimates P-value for a given value of enrichment score (ES)
    by comparing it with a set of scores ES_null, which was calculated for
    randomly generated gene sets."""

    n = len(ES_null)
    location = np.searchsorted(ES_null, ES)
    P = location/N
    return min(P,1-P)

# Multiple Hypothesis testing

def multiple_hypotesis_testing(D, C, S_sets, p_exp=1, random_sets=1000):
    """Performs Multiple Hypotesis Testing."""
    N, k = D.shape
    # generate random gene sets
    Pi_sets = []
    ES = []
    for i in range(random_sets):
        Pi_sets.append(tuple(np.random.randint(0,2) for i in range(k)))
    # calculate enrichment scores
    for i in range(len(S_sets)):
        ES.append(enrichment_score(D,C,S, p_exp))
        for pi in Pi_sets:
            ES_pi[i].append(enrichment_score(D,pi,S,p_exp))
    # normalize ES

    # calculate FDRs

