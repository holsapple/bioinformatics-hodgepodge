"""
Final Project - This project simulates genetic drift of the Homo sapiens hemoglobin subunit beta.

Author: Raymond Holsapple
Ron username: rh1054

Last Modification: December 11, 2016
"""


import numpy as np
import matplotlib
import matplotlib.cm as cm
from scipy import stats
matplotlib.use('agg')
from matplotlib import pyplot as plt


def str2seq(dna_seq_str):
    """
    This function converts a string DNA sequence to a numpy array DNA sequence.
    :param dna_seq_str: (string) DNA sequence in string format.
    :return dna_seq_array: (numpy array) DNA sequence in numpy array format.
    """
    dna_seq_array = np.asarray(list(dna_seq_str))
    return dna_seq_array


def seq2str(dna_seq_array):
    """
    This function converts a numpy array DNA sequence to a string DNA sequence.
    :param dna_seq_array: (numpy array of strings) DNA sequence in numpy array format.
    :return dna_seq_str: (string) DNA sequence in string format.
    """
    dna_seq_str = ''
    for n in np.nditer(dna_seq_array):
        dna_seq_str += str(n)
    return dna_seq_str


def sim_mutations():
    """
    This function simulates the mutation process a specified number of generations with a specified number of mutations per generation.
    :return prot_mutation, snps: NumPy arrays of the amino acid mutation scores and the SNP counts after each generation.
    """
    prot_mutation = np.empty([num_gens, num_aa], dtype=int)
    snps = np.empty([3, num_gens], dtype=int)
    for k in range(num_gens):
        cds.mutate(gen_mutations)
        cds.count_snps()
        cds.diff()
        prot_mutation[k] = cds.aa_mut_score
        snps[0, k] = cds.snp
        snps[1, k] = cds.syn_snp
        snps[2, k] = cds.nonsyn_snp
    return prot_mutation, snps


def mc_sim_mutations():
    """
    This function runs a Monte Carlo simulation of the function sim_mutations() and computes an element-wise mean of the resulting mutation score matrices.
    :return: (NumPy arrays) Element-wise means & standard errors of (a) the amino acid mutation scoring matrices and (b) the SNPs.
    """
    mciter_prot_mutation = np.empty([monte_carlo_its, num_gens, num_aa], dtype=int)
    mciter_snps = np.empty([monte_carlo_its, 3, num_gens], dtype=int)
    print()
    print('We have {} Monte Carlo simulation iterations to run.'.format(monte_carlo_its))
    print('-----------------------------------------------------')
    print()
    for k in range(monte_carlo_its):
        print('Iteration {} is running right now.'.format(k + 1))
        cds.__init__(sequence)
        mciter_prot_mutation[k], mciter_snps[k] = sim_mutations()
    prot_mutation = np.mean(mciter_prot_mutation, axis=0)
    pm_stderr = stats.sem(mciter_prot_mutation, axis=0)
    snps = np.mean(mciter_snps, axis=0)
    snps_stderr = stats.sem(mciter_snps, axis=0)
    return prot_mutation, pm_stderr, snps, snps_stderr


def plot_map(label='unk'):
    """
    This function plots a heat map of the amino acid mutation scores.
    :param label: (string) A label for the filename of the plot which will be generated.
    :return: None
    """
    mutation_x = list(range(num_aa))
    mutation_y = list(range(num_gens))
    xx, yy = np.meshgrid(mutation_x, mutation_y)
    plt.figure()
    plt.pcolormesh(xx, yy, protein_mutation, cmap=cm.copper)
    plt.axis([xx.min(), xx.max(), yy.min(), yy.max()])
    plt.colorbar()
    plt.xlabel('amino acid index', fontweight='bold')
    plt.ylabel('generation number', fontweight='bold')
    plt.title('Mutation Heat Map for Homo sapiens\nHemoglobin subunit $\\beta$', fontweight='bold')
    plt.savefig('mutation_map_{}.png'.format(label), dpi=300)


def plot_diff_row(label='unk'):
    """
    This function plots the element-wise amino acid mutation scores for a specified generation of the simulation.
    :param label: (string) A label for the filename of the plot which will be generated.
    :return: None
    """
    row_y = protein_mutation[diff_row]
    row_x = list(range(num_aa))
    if monte_carlo and inside_mc:
        num_sig = num_sd * pm_stderr[diff_row]
        max_err = max(num_sig)
        plt.figure()
        plt.errorbar(row_x, row_y, yerr=num_sig, marker='o', ms=2, color='black', linewidth=1, elinewidth=0.2)
        plt.axis([min(row_x), max(row_x), (min(row_y) - max_err - 0.5), (max(row_y) + max_err + 0.5)])
        plt.xlabel('amino acid index', fontweight='bold')
        plt.ylabel('amino acid mutation score', fontweight='bold')
        plt.title('Amino Acid Mutation Scores for Homo sapiens\nHemoglobin subunit $\\beta$ (Generation {})'.format(diff_row), fontweight='bold')
        plt.savefig('aa_mut_scores_{}.png'.format(label), dpi=300)
        plt.figure()
        plt.errorbar(row_x[0:51], row_y[0:51], yerr=num_sig[0:51], marker='o', ms=2, color='black', linewidth=1, elinewidth=0.2)
        plt.axis([min(row_x[0:51]), max(row_x[0:51]), (min(row_y[0:51]) - max_err - 0.5), (max(row_y[0:51]) + max_err + 0.5)])
        plt.xlabel('amino acid index', fontweight='bold')
        plt.ylabel('amino acid mutation score', fontweight='bold')
        plt.title('Amino Acid Mutation Scores for Homo sapiens\nHemoglobin subunit $\\beta$ (Generation {})'.format(diff_row), fontweight='bold')
        plt.savefig('aa_mut_scores_third_{}.png'.format(label), dpi=300)
    else:
        plt.figure()
        plt.plot(row_x, row_y, marker='o', ms=2, color='black', linewidth=1)
        plt.axis([min(row_x), max(row_x), (min(row_y) - 1), (max(row_y) + 1)])
        plt.xlabel('amino acid index', fontweight='bold')
        plt.ylabel('amino acid mutation score', fontweight='bold')
        plt.title('Amino Acid Mutation Scores for Homo sapiens\nHemoglobin subunit $\\beta$ (Generation {})'.format(diff_row), fontweight='bold')
        plt.savefig('aa_mut_scores_{}.png'.format(label), dpi=300)


def plot_snps(label='unk'):
    """
    This function plots the SNPs (total, synonymous, nonsynonymous) counts following a range of generational mutations.
    :param label: (string) A label for the filename of the plot which will be generated.
    :return: None
    """
    plt.figure()
    plt.plot(snp[0], '-k', linewidth=2, label='Total SNPs')
    plt.plot(snp[1], '-g', linewidth=2, label='Synonymous SNPs')
    plt.plot(snp[2], '-r', linewidth=2, label='Nonsynonymous SNPs')
    plt.xlim(0, (num_gens - 1))
    plt.xlabel('generation index', fontweight='bold')
    plt.ylabel('number of SNPs', fontweight='bold')
    plt.title('SNP Comparison for Homo sapiens\nHemoglobin subunit $\\beta$ ({} Generations)'.format(num_gens), fontweight='bold')
    plt.legend(loc=7, prop={'size': 10}, numpoints=1)
    plt.savefig('snp_comparison_{}.png'.format(label), dpi=300)


"""
This is one bodacious CDS class. Run my code and see how gnarly it operates.
"""
class CDS:
    def __init__(self, seq):
        self.seq = str2seq(seq)
        self.mutated_seq = np.copy(self.seq)
        self.syn_snp = 0
        self.nonsyn_snp = 0
        self.snp = 0
        self.aa_mut_count = np.zeros((int(len(self.seq) / 3),), dtype=np.int)
        self.aa_mut_score = np.zeros((int(len(self.seq) / 3),), dtype=np.int)
        self.mut_ind_list = []

    def prot(self, sequence):
        """
        This function translates a nucleotide sequence into an amino acid sequence.
        :param sequence: (NumPy array) The nucleotide sequence.
        :return: A NumPy array amino acid sequence.
        """
        seq = seq2str(sequence)
        position = 0
        protein = ''
        while position < len(seq):
            codon = seq[position:position + 3]
            protein += codons[codon]
            position += 3
        return str2seq(protein)

    def mutate(self, n):
        """
        This function mutates the nucleotide sequence a specified number of times in each generation.
        :param n: (int) The number of nucleotide mutations to attempt.
        :return: None
        """
        seq_len = len(self.seq)
        self.mut_ind_list = []
        mutation_count = 0
        while mutation_count < n:
            mut_ind = np.random.randint(0, seq_len - 1)
            self.mut_ind_list.append(mut_ind)
            mut_nuc = self.mutated_seq[mut_ind]
            mut_choices = np.asarray(['transition', 'transversion'])
            mut_type = np.random.choice(mut_choices, p=[0.75, 0.25])
            if mut_type == 'transition':
                mutated_nuc = t_ition[mut_nuc]
            else:
                mutated_nuc = np.random.choice(t_version[mut_nuc], p=[0.5, 0.5])
            if mut_ind % 3 == 0:
                new_codon = str(mutated_nuc) + str(self.mutated_seq[mut_ind + 1]) + str(self.mutated_seq[mut_ind + 2])
                if (new_codon != 'TAA') and (new_codon != 'TAG') and (new_codon != 'TGA'):
                    self.mutated_seq[mut_ind] = mutated_nuc
                    mutation_count += 1
            elif mut_ind % 3 == 1:
                new_codon = str(self.mutated_seq[mut_ind - 1]) + str(mutated_nuc) + str(self.mutated_seq[mut_ind + 1])
                if (new_codon != 'TAA') and (new_codon != 'TAG') and (new_codon != 'TGA'):
                    self.mutated_seq[mut_ind] = mutated_nuc
                    mutation_count += 1
            else:
                new_codon = str(self.mutated_seq[mut_ind - 2]) + str(self.mutated_seq[mut_ind - 1]) + str(mutated_nuc)
                if (new_codon != 'TAA') and (new_codon != 'TAG') and (new_codon != 'TGA'):
                    self.mutated_seq[mut_ind] = mutated_nuc
                    mutation_count += 1

    def count_snps(self):
        """
        This function computes the current number of SNPs (total, synonymous, and nonsynonymous) for the mutated nucleotide sequence.
        :return: None
        """
        self.snp = np.sum(self.seq != self.mutated_seq)
        nuc_diffs = np.where(self.seq != self.mutated_seq)
        self.nonsyn_snp = 0
        for mut_ind in nuc_diffs[0]:
            if mut_ind % 3 == 0:
                seq_codon = str(self.seq[mut_ind]) + str(self.seq[mut_ind + 1]) + str(self.seq[mut_ind + 2])
                mutated_codon = str(self.mutated_seq[mut_ind]) + str(self.mutated_seq[mut_ind + 1]) + str(self.mutated_seq[mut_ind + 2])
            elif mut_ind % 3 == 1:
                seq_codon = str(self.seq[mut_ind - 1]) + str(self.seq[mut_ind]) + str(self.seq[mut_ind + 1])
                mutated_codon = str(self.mutated_seq[mut_ind - 1]) + str(self.mutated_seq[mut_ind]) + str(self.mutated_seq[mut_ind + 1])
            else:
                seq_codon = str(self.seq[mut_ind - 2]) + str(self.seq[mut_ind - 1]) + str(self.seq[mut_ind])
                mutated_codon = str(self.mutated_seq[mut_ind - 2]) + str(self.mutated_seq[mut_ind - 1]) + str(self.mutated_seq[mut_ind])
            if codons[seq_codon] != codons[mutated_codon]:
                self.nonsyn_snp += 1
        self.syn_snp = self.snp - self.nonsyn_snp

    def diff(self):
        """
        This function computes the element-wise amino acid mutation scores for the mutated amino acid sequence.
        :return: None
        """
        seq_prot = self.prot(self.seq)
        mutated_seq_prot = self.prot(self.mutated_seq)
        prot_mut_ind = []
        for k in range(len(self.mut_ind_list)):
            prot_mut_ind.append(int(self.mut_ind_list[k] / 3))
        for ind in prot_mut_ind:
            if seq_prot[ind] != mutated_seq_prot[ind]:
                self.aa_mut_count[ind] += 1
                if aa_polarity[seq_prot[ind]] == aa_polarity[mutated_seq_prot[ind]]:
                    self.aa_mut_score[ind] += 1
                else:
                    self.aa_mut_score[ind] += 4

if __name__ == '__main__':
    # get polarity data and create a dictionary #
    aa_polarity = {}
    with open('polarity.csv') as f:
        for line in f:
            if 'nonpolar' in line:
                aa_polarity[line[0]] = 0
            else:
                aa_polarity[line[0]] = 1

    # Homo sapiens Hemoglobin subunit beta CDS #
    heme_cds = 'GTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCAC'

    # codon to amino acid translation dictionary #
    codons = {'TTT': 'F', 'TTC': 'F', 'TTY': 'F',
              'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CTN': 'L',
              'YTN': 'L', 'TTR': 'L',
              'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATH': 'I',
              'ATG': 'M',
              'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GTN': 'V',
              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TCN': 'S',
              'AGT': 'S', 'AGC': 'S', 'WSN': 'S', 'AGY': 'S',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCN': 'P',
              'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ACN': 'T',
              'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GCN': 'A',
              'TAT': 'Y', 'TAC': 'Y', 'TAY': 'Y',
              'TAA': 'X', 'TAG': 'X', 'TGA': 'X', 'TRR': 'X', 'TAR': 'X', 'NNN': 'X',
              'CAT': 'H', 'CAC': 'H', 'CAY': 'H',
              'CAA': 'Q', 'CAG': 'Q', 'CAR': 'Q',
              'AAT': 'N', 'AAC': 'N', 'AAY': 'N',
              'AAA': 'K', 'AAG': 'K', 'AAR': 'K',
              'GAT': 'D', 'GAC': 'D', 'GAY': 'D',
              'GAA': 'E', 'GAG': 'E', 'GAR': 'E',
              'TGT': 'C', 'TGC': 'C', 'TGY': 'C',
              'TGG': 'W',
              'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'MGN': 'R', 'CGN': 'R', 'AGR': 'R',
              'AGA': 'R', 'AGG': 'R',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGN': 'G'}

    # nucleotide transition dictionary #
    t_ition = {'A': 'G',
               'G': 'A',
               'C': 'T',
               'T': 'C'}

    # nucleotide transversion dictionary #
    t_version = {'A': ['C', 'T'],
                 'G': ['C', 'T'],
                 'C': ['A', 'G'],
                 'T': ['A', 'G']}

    # -------------------------------------------------------------------- #
    #                 The main simulation code is below...                 #
    # -------------------------------------------------------------------- #

    sequence = heme_cds
    cds = CDS(sequence)
    num_aa = int(len(sequence) / 3)
    gen_mutations = 5
    num_gens = 101
    diff_row = 300
    num_sd = 3
    monte_carlo_its = 100
    # monte_carlo = False
    monte_carlo = True
    inside_mc = False

    protein_mutation, snp = sim_mutations()
    plot_map('{}gens'.format(num_gens))
    plot_diff_row('gen{}'.format(diff_row))
    plot_snps('{}gens'.format(num_gens))

    if monte_carlo:
        inside_mc = True
        protein_mutation, pm_stderr, snp, snp_stderr = mc_sim_mutations()
        plot_map('{}mciters'.format(monte_carlo_its))
        plot_diff_row('mcgen{}'.format(diff_row))
        plot_snps('{}mciters'.format(monte_carlo_its))

    print()
    print()
    print('DONE')
