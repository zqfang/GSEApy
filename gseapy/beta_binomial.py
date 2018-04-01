# Implements a beta binomial distribution using the scipy discrete random
# variable API.
# Requires scipy 0.14.0 or newer; recommend scipy 0.15.0 or newer
# Distribution functions unapologetically adapted from Wikipedia
# http://en.wikipedia.org/wiki/Beta-binomial_distribution

from numpy import exp, floor, ceil, sqrt, r_, maximum, where, sum, zeros, linspace
from numpy.lib.function_base import vectorize
from scipy.stats._distn_infrastructure import rv_discrete, rv_frozen
from scipy.stats import beta, binom
from scipy.special import betaln, gammaln, bdtr, bdtrc, btdtri, bdtrik, btdtr
try:
    from scipy.special import entr # scipy[>=0.15.0]
except ImportError:
    from scipy.special import xlogy
    entr = lambda x: -xlogy(x, x)

class betabinom_gen(rv_discrete):
    def __init__(self, **kwargs):
        super(betabinom_gen, self).__init__(**kwargs)
        self._sf = vectorize(self._sf_single, otypes='d')
        self._sf.nin = self.numargs + 1

    def _rvs(self, n, a, b):
        p = self._random_state.beta(a, b, self._size)
        return self._random_state.binomial(n, p, self._size)

    def _argcheck(self, n, a, b):
        self.b = n
        return (n >= 0) & beta._argcheck(a, b)

    def _logpmf(self, x, n, a, b):
        k = floor(x)
        combiln = gammaln(n + 1) - (gammaln(k + 1) + gammaln(n - k + 1))
        return combiln + betaln(k + a, n - k + b) - betaln(a, b)

    def _pmf(self, x, n, a, b):
        return exp(self._logpmf(x, n, a, b))

    def _cdf_single(self, x, n, a, b):
        k = floor(x)
        p = linspace(0, 1, num = 10001)
        bta = btdtr(a, b, p)
        p_med = (p[:-1] + p[1:]) / 2
        bta_med = bta[1:] - bta[:-1]
        vals = (bdtr(k, n, p_med) * bta_med).sum(axis = -1)
        return vals

    def _sf_single(self, x, n, a, b):
        k = floor(x)
        p = linspace(0, 1, num = 10001)
        bta = btdtr(a, b, p)
        p_med = (p[:-1] + p[1:]) / 2
        bta_med = bta[1:] - bta[:-1]
        vals = (bdtrc(k, n, p_med) * bta_med).sum(axis = -1)
        return vals

    def _ppf(self, q, n, a, b):
        # This brute-force solution is far from optimal and will be slow for
        # large n.
        k = 0
        d = n
        while (d > 1).any():
            d = ceil(d / 2)
            temp = self._cdfvec(k, n, a, b)
            k += where(temp <= q, d, -d)
        temp = self._cdfvec(k, n, a, b)
        vals1 = maximum(k - 1, 0)
        return where(temp <= q, k, vals1)

    def _stats(self, n, a, b, moments = 'mv'):
        e_p = float(a) / (a + b)
        e_q = 1 - e_p
        mu = n * e_p
        var = n * (a + b + n) * e_p * e_q / (a + b + 1)
        g1, g2 = None, None
        if 's' in moments:
            g1 = 1.0 / sqrt(var)
            g1 *= (a + b + 2 * n) * (b ** 2 - a ** 2)
            g1 /= (a + b + 2)
        if 'k' in moments:
            g2 = float(a + b)
            g2 *= (a + b - 1 + 6 * n)
            g2 += 3 * a * b * (n - 2)
            g2 += 6 * n ** 2
            g2 -= float(3 * e_p * b * n * (6 - n))
            g2 -= float(18 * e_p * e_q * n ** 2)
            g2 *= (a + b) ** 2 * (1 + a + b)
            g2 /= (n * a * b * (a + b + 2) * (a + b + 3) * (a + b + n))
        return mu, var, g1, g2

    def _entropy(self, n, a, b):
        k = r_[0:n + 1]
        vals = self._pmf(k, n, a, b)
        return sum(entr(vals), axis = 0)

betabinom = betabinom_gen(name="betabinom")