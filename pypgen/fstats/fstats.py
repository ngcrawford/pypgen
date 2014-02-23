#!/usr/bin/env python
# encoding: utf-8


import math

## UTIlITY STATS:
def de_NaN_list(l):
    """Remove NaN values from list.

    Returns list without NaNs.

    Parameters:

        a : array_like

    Returns

        m : array_like"""

    return [i for i in l if math.isnan(i) != True]


def _mean_(l):
    """Compute the arithmetic mean along the specified axis.

    Returns the average of the array elements.

    Parameters

        a : array_like

    Returns

        mean : array_like"""

    l = de_NaN_list(l)
    mean = sum(l) * 1.0 / len(l)
    return mean


def _mean_variance_(l):
    """Compute the variance of the 1D list or array.

    Returns the variance of the array elements, a measure of the spread of a
    distribution.

    Parameters

        a : array_like

    Returns

        variance : array_like

    Note: Doc string derived from numpy.var()"""

    mean = _mean_(l)
    var = _mean_(map(lambda x: (x - mean) ** 2, l))
    return var


def _stdev_(l):
    """Compute the standard deviation along the specified axis.

    Returns the standard deviation, a measure of the spread of a distribution,
    of the array elements.

    Parameters

        l : array_like
            Calculate the standard deviation of these values.

    Returns

        standard_deviation : array_like

    Note: Doc string derived from numpy.std()"""

    var = _mean_variance_(l)
    stdev = math.sqrt(var)
    return stdev


def harmonic_mean(l):
    """Calculates harmonic mean from list of integers

    Parameters

        l : array_like
            Calculate the harmonic mean of these values.

    Returns

        harmonic mean : array_like"""

    fractional_counts = sum([1.0 / v for v in l])
    harmonic_mean = float(len(l)) / fractional_counts

    return harmonic_mean


def harmonic_mean_chao(l):
    """Calculates the harmonic mean following the method suggested by
    Anne Chao. The formula is: 1/[(1/A)+var(D)(1/A)**3]. Used for
    calculating multilocus Dest.

    Parameters

        l : array_like
            Calculate Chao's harmonic mean of these values.

    Returns

        Chao's harmonic mean : array_like"""

    count = float(len(l))
    if count == 0.0:
        return 0.0

    elif sum(l) / count == 0:
        return 0.0

    else:
        A = sum(l) / count
        varD = sum([(v - A) ** 2 for v in l]) / count
        harmonic_mean_chao = 1 / ((1 / A) + (varD) * pow((1 / A), 3))
        return harmonic_mean_chao


def Hs_prime_est(allele_freqs, n):
    """Calculate corrected Hs: the mean within-subpopulation
    heterozygosity (Nei and Chesser 1983).

    Parameters

        allele_freqs : array_like
            These values contain the allele freqeuncies

        n : int,float
            The number of populations

    Returns

        H's : array_like"""

    n = float(n)
    Hj = [1.0 - sum([freq ** 2 for freq in pop]) for pop in allele_freqs]
    Hs_prime_est = (1 / n) * sum(Hj)
    return Hs_prime_est


def Hs_est(Hs_prime_est, harm_mean):
    """Basic Equation: ((2*N_harmonic)/(2*N_harmonic-1))*Hs"""

    Hs = Hs_prime_est
    Hs_est = ((2.0 * harm_mean) / (2.0 * harm_mean - 1.0)) * Hs
    return Hs_est


def Ht_prime_est(allele_freqs, n):
    """Calculate corrected Ht: the heterozygosity of the pooled
    subpopulations (Nei and Chesser 1983)"""

    n = float(n)
    inner = [((1 / n) * sum(allele_list)) ** 2 for allele_list in zip(*allele_freqs)]
    Ht_prime_est = 1.0 - sum(inner)
    return Ht_prime_est


def Ht_est(Ht_p_est, Hs_est, harm_mean, n):
    """Basic Equation: Ht+Hs_est/(2*N_harmonic*n)"""

    n = float(n)
    Ht = Ht_p_est
    Ht_est = Ht + Hs_est / (2.0 * harm_mean * n)
    return Ht_est


def Gst_est(Ht_est, Hs_est):
    """Basic Equation: Gst = (Ht-Hs)/Ht"""

    if Ht_est == 0.0:
        return 0.0
    else:
        Gst_est = (Ht_est - Hs_est) / Ht_est
        return Gst_est


def G_prime_st_est(Ht_est, Hs_est, Gst_est, n):
    """Basic Equation: G'st = Gst/Gst max = (Gst*(k-1+Hs))/((k-1)*(1-Hs))
       from Hedrick 2005)"""

    n = float(n)
    if (n - 1.0) * (1.0 - Hs_est) == 0:
        return 0.0

    else:
        G_prime_st = (Gst_est * (n - 1.0 + Hs_est)) / ((n - 1.0) * (1.0 - Hs_est))
        return G_prime_st


def G_double_prime_st_est(Ht_est, Hs_est, n):
    """Basic Equation: G''st = k*(HT-HS)/((k*HT-HS)*(1-HS)
        from Meirmans and Hedrick (2011)"""

    n = float(n)
    if (n * Ht_est - Hs_est) * (1 - Hs_est) == 0.0:
        return 0.0
    else:
        G_double_prime_st_est = n * (Ht_est - Hs_est) / ((n * Ht_est - Hs_est) * (1 - Hs_est))
        return G_double_prime_st_est


def D_est(Ht_est, Hs_est, n):
    """Basic Equation: ((Ht-Hs)/(1.0-Hs))*(n/(n-1))"""

    n = float(n)
    if ((1.0 - Hs_est)) * (n / (n - 1)) == 0.0:
        return 0.0
    else:
        D_est = ((Ht_est - Hs_est) / (1.0 - Hs_est)) * (n / (n - 1))
        return D_est

###########################################
#
#     MULTILOCUS FUNCTIONS:
#


def multilocus_Gst_est(Ht_est, Hs_est):
    """Averages Ht_est and Hs_est across loci before calculating Gst."""

    stdev = _stdev_(map(Gst_est, Ht_est, Hs_est))

    Ht_est = sum(Ht_est) / float(len(Ht_est))
    Hs_est = sum(Hs_est) / float(len(Hs_est))

    if Ht_est == 0.0:
        return (0.0, 0.0)
    else:
        ml_Gst_est = (Ht_est - Hs_est) / Ht_est
        return (ml_Gst_est, stdev)


def multilocus_G_prime_st_est(Ht_est, Hs_est, n):
    """Averages across loci before calculating G'st."""

    Gst_est_list = map(Gst_est, Ht_est, Hs_est)
    stdev = _stdev_(map(G_prime_st_est, Ht_est, Hs_est, Gst_est_list, [float(n)] * len(Ht_est)))

    n = float(n)
    # Ht_est_ = sum(Ht_est) / float(len(Ht_est))
    Hs_est_ = sum(Hs_est) / float(len(Hs_est))
    G_est = multilocus_Gst_est(Ht_est, Hs_est)[0]

    if (n - 1.0) * (1.0 - Hs_est_) == 0:
        return (0.0, 0.0)
    else:
        G_prime_st = (G_est * (n - 1.0 + Hs_est_)) / ((n - 1.0) * (1.0 - Hs_est_))
        return (G_prime_st, stdev)


def multilocus_G_double_prime_st_est(Ht_est, Hs_est, n):
    """Averages across loci before calculating G''st."""

    stdev = _stdev_(map(G_double_prime_st_est, Ht_est, Hs_est, [float(n)] * len(Ht_est)))

    n = float(n)
    Ht_est = sum(Ht_est) / float(len(Ht_est))
    Hs_est = sum(Hs_est) / float(len(Hs_est))

    # G_est = multilocus_Gst_est(Ht_est,Hs_est)
    if (n * Ht_est - Hs_est) * (1.0 - Hs_est) == 0:
        return (0.0, 0.0)
    else:
        ml_G_double_prime_st_est = n * (Ht_est - Hs_est) / ((n * Ht_est - Hs_est) * (1.0 - Hs_est))

        return (ml_G_double_prime_st_est, stdev)


def multilocus_D_est(Ht_est, Hs_est, n):
    """Calculate Dest values using Anne Chao's harmonic mean."""

    pairs = zip(Ht_est, Hs_est)
    Dest_values = [D_est(pair[0], pair[1], 2) for pair in pairs]
    ml_D_est = harmonic_mean_chao(Dest_values)
    stdev = _stdev_(Dest_values)
    return (ml_D_est, stdev)
