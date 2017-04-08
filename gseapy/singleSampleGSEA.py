'''
This file contains script for emulating the ssGSEA algorithm from Barbie et al. 2009
'''

import numpy as np
import random
from scipy.stats import norm

def calculate_enrichment_score(gene_set, expressions, omega):
    """
    Given a gene set, a map of gene names to expression levels, and a weight omega, returns the ssGSEA
    enrichment score for the gene set as described by *D. Barbie et al 2009*
    :requires: every member of gene_set is a key in expressions
    :param gene_set: a set of gene_names in the set
    :type gene_set: set
    :param expressions: a dictionary mapping gene names to their expression values
    :type expressions: dict
    :param omega: the weighted exponent on the :math:`P^W_G` term.
    :type omega: float
    :returns: an array representing the intermediate Enrichment Scores for each step along the sorted gene list.
              To find the total enrichment score, take the sum of all values in the array.
    """

    #first sort by absolute expression value, starting with the highest expressed genes first
    keys_sorted = sorted(expressions, key=expressions.get, reverse=True)

    #values representing the ECDF of genes in the geneset
    P_GW_numerator = 0
    P_GW_denominator = 0

    #determining denominator value
    i = 1 #current rank stepping through listing of sorted genes
    for gene in keys_sorted:
        if gene in gene_set:
            P_GW_denominator += i ** omega
        i += 1

    P_GW = lambda : P_GW_numerator / P_GW_denominator

    #values representing the ECDF of genes not in the geneset
    P_NG_numerator = 0
    P_NG_denominator = len(expressions) - len(gene_set)
    P_NG = lambda : P_NG_numerator / P_NG_denominator

    #integrate different in P_GW and P_NG
    i = 1 #current index in the traversal of sorted genes
    scores = []
    for gene in keys_sorted:
        if gene in gene_set:
            P_GW_numerator += i ** omega
        else:
            P_NG_numerator += 1

        scores.append(P_GW() - P_NG())
        i += 1

    return scores


'''
This module contains methods for simulating data given a hybrid data model
'''

def simulate_data(models, master_gene_name, sample_profiles, n):
    """
    Given a list of gaussian mixture models for genes, a master gene, and a list of sample profiles,
    returns n simulated phenotypes for each sample profile. This is done by selection of sample profiles
    without replacement n times. Note this and the actual simulation is a stochastic process which looks
    at the combined likehood for each class according to the mixture model.
    :requires: n is positive and less than or equal to len(sample_profiles)
    :param models: a mapping of gene names to the gaussian mixture model (GMM) for each gene.
    :type models: dict
    :param master_gene_name: the gene to be used as the master gene for phenotype simulation
    :type master_gene_name: str
    :param sample_profiles: a list of samples (data_models.sample)
    :type sample_profiles: list
    :param n: the number of phenotypes to simulate
    :type n: int
    :return: a dict mapping sample id's to phenotype classes. This phenotype does not necessarily correspond
             to any real phenotype but rather the components of the mixture model of the master gene. class 0 represents
             the component in the model with a smaller expression mean while class 1 represents the one with the larger.
    """

    master_gene_model = models[master_gene_name]

    #for selecting random front n values
    random.shuffle(sample_profiles)

    classifications = {}
    for i in range(0, n):
        cur_profile = sample_profiles[i]
        id = cur_profile.id
        cur_intensity = cur_profile.profiles[master_gene_name].intensity

        #classification of phenotype for this profile
        proportions = master_gene_model.weights_
        mus = master_gene_model.means_
        sigmas = [x[0] ** 0.5 for x in master_gene_model.covars_]

        #setting 0/1 to the proper low/high peaks. 0 = smaller mu
        if(mus[0] < mus[1]):
            mu0 = sum(mus[0])
            sigma0 = sigmas[0]
            proportions0 = proportions[0]
            mu1 = sum(mus[1])
            sigma1 = sigmas[1]
            proportions1 = proportions[1]
        else:
            mu0 = sum(mus[1])
            sigma0 = sigmas[1]
            proportions0 = proportions[1]
            mu1 = sum(mus[0])
            sigma1 = sigmas[0]
            proportions1 = proportions[0]

        #figuring out probability for each sample being in each peak
        density0 = norm.pdf(cur_intensity, loc=mu0, scale=sigma0)
        density1 = norm.pdf(cur_intensity, loc=mu1, scale=sigma1)

        probability1 = proportions1 * density1 / (proportions1 * density1 + proportions0 * density0)

        if random.random() <= probability1:
            classifications[id] = 1
        else:
            classifications[id] = 0

    return classifications