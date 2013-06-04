# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

from cmath import *

_g = 7
_p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
     771.32342877765313, -176.61502916214059, 12.507343278686905,
     -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]


def lngamma(z):
    """
    Lanchos approximation of log((z-1)!)
    Reference: http://en.wikipedia.org/wiki/Lanczos_approximation
    """
    z = complex(z)
    z -= 1
    x = _p[0]
    for i in range(1, _g+2):
        x += _p[i]/(z+i)
    t = z + _g + 0.5
    return 0.9189385332046727 + (z + 0.5) * log(t) - t + log(x) # 0.9189385332046727 = log(sqrt(2*pi))

def lnfactorial(n):
    return 0 if n < 1 else lngamma(n + 1)


def lncombination(n, k):
    return lnfactorial(n) - lnfactorial(k) - lnfactorial(n - k)


def hypergeometric_probability(k, n, K, N):
    return exp(lncombination(K, k) + lncombination(N - K, n - k) - lncombination(N, n)).real


def hypergeometric_test(k, n, K, N):
    um = min(n, K)
    lm = max(0, n + K - N)
    eps = 1e-10

    if um == lm:
        return(1.0, 1.0, 1.0)

    prob = hypergeometric_probability(k, n, K, N)
    
    two_tail = 0
    for i in range(lm, um + 1):
        p = hypergeometric_probability(i, n, K, N)
        if p < prob + eps:
            two_tail += p

    two_tail = two_tail if two_tail < 1.0 else 1.0
    return two_tail

if __name__ == '__main__':
    pass