#!/usr/bin/env python
# encoding: utf-8

import math


def de_NaN_list(l):
    return [i for i in l if math.isnan(i) != True]


def harmonic_mean(values):
    """calculates harmonic mean from list of integers"""

    fractional_counts = sum([1.0 / v for v in values])
    harmonic_mean = float(len(values)) / fractional_counts

    return harmonic_mean


def harmonic_mean_chao(values):
    """Calculates the harmonic mean following the method suggested by
    Anne Chao. The formula is: 1/[(1/A)+var(D)(1/A)**3]. Used for
    calculating multilocus Dest."""

    count = float(len(values))
    if count == 0.0:
        return 0.0
    elif sum(values) / count == 0:
        return 0.0
    else:
        A = sum(values) / count
        varD = sum([(v - A) ** 2 for v in values]) / count
        harmonic_mean_chao = 1 / ((1 / A) + (varD) * pow((1 / A), 3))
        return harmonic_mean_chao


def Hs_prime_est(allele_freqs, n):
    """Calculate corrected Hs: the mean within-subpopulation
    heterozygosity (Nei and Chesser 1983)."""

    n = float(n)
    Hj = [1.0 - sum([freq ** 2 for freq in pop[:4]]) for pop in allele_freqs]
    Hs_prime_est = (1 / n) * sum(Hj)
    return Hs_prime_est


def Hs_est(Hs_prime_est, harm_mean):
    """ Basic Equation: ((2*N_harmonic)/(2*N_harmonic-1))*Hs"""

    Hs = Hs_prime_est
    Hs_est = ((2.0 * harm_mean) / (2.0 * harm_mean - 1.0)) * Hs
    return Hs_est


def Ht_prime_est(allele_freqs, n):
    """Calculate corrected Ht: the heterozygosity of the pooled
    subpopulations (Nei and Chesser 1983)"""

    n = float(n)
    inner = [((1 / n) * sum(allele_list)) ** 2 for allele_list in zip(*allele_freqs)[:4]]
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

    n = float(n)
    if (n - 1.0) * (1.0 - Hs_est) == 0:
        return 0.0

    else:
        G_prime_st = (Gst_est * (n - 1.0 + Hs_est)) / ((n - 1.0) * (1.0 - Hs_est))
        return G_prime_st


def G_double_prime_st_est(Ht_est, Hs_est, n):
    """Basic Equation: G''st = k*(HT-HS)/((k*HT-HS)*(1-HS)"""

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

    Ht_est = sum(Ht_est) / float(len(Ht_est))
    Hs_est = sum(Hs_est) / float(len(Hs_est))

    if Ht_est == 0.0:
        return 0.0
    else:
        Gst_est = (Ht_est - Hs_est) / Ht_est
        return Gst_est


def multilocus_G_prime_st_est(Ht_est, Hs_est, n):
    """Averages across loci before calculating G'st."""

    n = float(n)
    # Ht_est_ = sum(Ht_est) / float(len(Ht_est))
    Hs_est_ = sum(Hs_est) / float(len(Hs_est))
    G_est = multilocus_Gst_est(Ht_est, Hs_est)

    if (n - 1.0) * (1.0 - Hs_est_) == 0:
        return 0.0
    else:
        G_prime_st = (G_est * (n - 1.0 + Hs_est_)) / ((n - 1.0) * (1.0 - Hs_est_))
        return G_prime_st


def multilocus_G_double_prime_st_est(Ht_est, Hs_est, n):
    """Averages across loci before calculating G''st."""

    n = float(n)
    Ht_est = sum(Ht_est) / float(len(Ht_est))
    Hs_est = sum(Hs_est) / float(len(Hs_est))

    # G_est = multilocus_Gst_est(Ht_est,Hs_est)
    if (n * Ht_est - Hs_est) * (1.0 - Hs_est) == 0:
        return 0.0
    else:
        G_double_prime_st_est = n * (Ht_est - Hs_est) / ((n * Ht_est - Hs_est) * (1.0 - Hs_est))
        return G_double_prime_st_est


def multilocus_D_est(Ht_est, Hs_est, n):
    """Calculate Dest values using Anne Chao's harmonic mean."""

    pairs = zip(Ht_est, Hs_est)
    Dest_values = [D_est(pair[0], pair[1], 2) for pair in pairs]
    D_est_ = harmonic_mean_chao(Dest_values)
    return D_est_
