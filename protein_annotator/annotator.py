"""
Author: Raymond Holsapple
Ron username: rh1054

This script annotates a genomic sequence similar to what PROKKA does. It takes a fasta file of contigs, finds all open reading frames,
translates the ORFs into proteins, blasts those proteings against the Swiss-Prot database, and then writes an annotated protein fasta file.

Last Update: November 12, 2016
"""

import sys
import re
import numpy as np
from subprocess import Popen, PIPE
from collections import namedtuple
from seq_tools import orf_finder

column_names = 'qseqid pident stitle'  # These are the columns to be printed to standard out when blastp is run.
blast_record = namedtuple('Blast', column_names)


def create_contigs_dictionary(filename):
    """
    This function creates a dictionary of contigs where the keys are the contig headers and the values are the sequences.
    :param filename: (string) Name of the contig fasta file
    :return contigs: (dictionary) Dictionary of contigs
    """
    contigs = {}
    header = ''
    with open(filename, 'r') as file:  # open the contig file
        for line in file:  # iterate over each line in the contig file
            if line[0] == '>':
                header = line.rstrip()
                contigs[header] = ''  # if the line is a header, initialize a new dictionary entry
            else:
                contigs[header] += line.rstrip()  # if the line is a sequence, append it to the value of the key that was just made
    return contigs


def str2seq(dna_seq_str):
    """
    This function converts a string DNA sequence to a numpy array DNA sequence.
    :param dna_seq_str: (string) DNA sequence in string format
    :return dna_seq_array: (numpy array) DNA sequence in numpy array format
    """
    temp_list = []
    for c in dna_seq_str:  # iterate over the characters in the string dna_seq_str
        temp_list.append(c)  # append each character to the list temp_list
    dna_seq_array = np.asarray(temp_list)  # convert the list temp_list to a numpy array
    return dna_seq_array


def filter_short_orfs(dict, length):
    """
    This function filters out open reading frames that are shorter than the parameter 'length'.
    :param dict: (dictionary) dictionary of ORF indices for each contig and its reverse complement
    :param length: (integer) minimum length an ORF must be to remain in the dictionary
    :return orf_ind_dict: (dictionary) dictionary filtered of ORF indices shorter than 'length'
    """
    orf_ind_dict = {}
    for key, value in dict.items():  # iterate over the contig headers and the ORF indices array
        orf_dir_list = []
        for k in range(2):  # each contig has two lists of ORFs: 1) forward direction; 2) reverse complement
            orf_array = dict[key][k]  # orf_array is the full list of ORFs for each direction
            long_orf_list = []
            for orf in orf_array:  # iterate over each ORF in the full list of ORFs
                if (orf[1] - orf[0]) >= length:  # check  to see if the number of nucleotides in the ORF is >= 'length'
                    long_orf_list.append([orf[0], orf[1]])  # if the ORF is long enough, append it to a long ORF list
            orf_dir_list.append(long_orf_list)  # append the list of long ORFs for each direction to an array with two elements: one for each direction
        orf_ind_dict[key] = orf_dir_list  # for each contig, save the bi-directional array of long ORFs to the dictionary
    return orf_ind_dict


def revcomp(seq):
    """
    This function computes the reverse complement of a DNA sequence.
    :param seq: (string) DNA sequence
    :return rc_seq: (string) reverse complement of the DNA sequence
    """
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  # dictionary of nucleotide complements
    rc_seq = ''
    for nuc in seq:  # iterate over the nucleotides in the DNA sequence
        rc_seq += comp_dict[nuc]  # append the complementary nucleotide to the (non-reversed) complement sequence
    rc_seq = rc_seq[::-1]  # reverse the complement sequence
    return rc_seq


def call_process(command):
    """
    This function calls a command line process and returns standard out.
    :param command: (string) the command to be executed
    :return stdout: (string) standard out result of the command
    """
    p = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)  # this line runs the command
    stdout, stderr = p.communicate()  # retrieve the results of the command
    stdout = stdout.decode()  # decode the results of the command
    return stdout


def make_temp_file(sequence):
    """
    This function creates a temporary file which will be used as the query file in the blastp command.
    :param sequence: (string) the string that will be written to the temporary file
    :return: none
    """
    with open('temp.fa', 'w') as file:
        file.write(sequence)


def blast_against_db(sequence):
    """
    This function blasts a protein query sequence against a protein database.
    :param sequence: (string) protein query sequence that will be blasted
    :return best_hits_dict: (dictionary) dictionary of highest percent identity blast records for each protein
    """
    make_temp_file(sequence)  # make the temporary query file
    command = 'blastp -query temp.fa -db /data/blastdb/uniprot_sprot -outfmt "6 {}" -num_threads 20'.format(column_names)  # define the blast command
    result = call_process(command)  # blast the query against the database and save the output in the variable 'result'
    best_hits_dict = {}
    for line in result.splitlines():  # iterate over each line in the result to pick out the ones that have the highest percent identity for each protein
        elements = line.split('\t')  # make a list of each line's contents
        br = blast_record(*elements)  # convert a line's contents into a blast_record named_tuple
        if br.qseqid not in best_hits_dict:
            best_hits_dict[br.qseqid] = br  # if a protein record isn't in the dictionary, add it
        elif float(best_hits_dict[br.qseqid].pident) < float(br.pident):
            best_hits_dict[br.qseqid] = br  # if a protein record's percent identity is higher than the one already saved, replace the lower one in the dictionary
    return best_hits_dict


def annotate_proteins(prot_seqs, prot_labels, prot_ids):
    """
    This function is the heart of the script. It creates the annotated protein fasta file.
    :param prot_seqs: (dictionary) dictionary of protein sequences: protein ID (key), protein sequence (value)
    :param prot_labels: (dictionary) dictionary of protein labels: protein ID (key), protein annotation (value)
    :param prot_ids: (list strings) ordered list of protein IDs
    :return: none
    """
    with open('annotated_proteins.faa', 'w') as file:  # create/open the protein fasta file
        prot_annotation = ''
        for id in prot_ids:  # iterate over the protein IDs which will be the headers of the fasta file
            rem = len(prot_seqs[id]) % 80  # this line is helpful for writing a fasta file with lines that have a maximum of 80 characters
            aa_count = 0  # a counter of amino acid number within a protein sequence
            prot_seq = ''
            for aa in prot_seqs[id]:  # iterate over each amino acid
                aa_count += 1  # add to the aa counter; preparing to append number 'aa_count' amino acid to the fasta protein sequence
                if (aa_count % 80) == 0:
                    prot_seq += aa + '\n'  # if the amino acid number is a multiple of 80, append a newline character after that amino acid
                else:
                    prot_seq += aa  # if the amino acid number isn't a multiple of 80, append it without a newline character after it
            if rem == 0:
                prot_annotation += id + ' ' + prot_labels[id] + '\n' + prot_seq  # if the number of amino acids is a multiple of 80, the sequence will already end with a newline character
            else:
                prot_annotation += id + ' ' + prot_labels[id] + '\n' + prot_seq + '\n'  # if the number of amino acids isn't a multiple of 80, the sequence will need a newline character
        file.write(prot_annotation)  # write the annotated protein fasta file


if __name__ == '__main__':
    filename = sys.argv[1]  # retrieve the contig file from the input
    contigs = create_contigs_dictionary(filename)  # create a contig dictionary: contig header (key), nucleotide sequence (value)

    temp_orf_ind_dict = {}  # a temporary ORF dictionary - includes short ORFs: contig header (key), ORF start/stop indices array (value)
    for key, value in contigs.items():  # iterate over the contig dictionary
        temp_orf_ind_dict[key] = orf_finder(str2seq(value))  # call the orf_finder function for each contig and use the results to bulid the temporary ORF indices dictionary
    orf_ind_dict = filter_short_orfs(temp_orf_ind_dict, 300)  # filter out all ORFs that are less than 300 nucleotides long

    orf_dict = {}  # a nucleotide sequence ORF dictionary: protein ID (key), ORF nucleotide sequence (value)
    count = 1  # a counter to keep track of the number of proteins that will be translated and blasted
    for key, value in orf_ind_dict.items():  # iterate over the ORF indices dictionary
        seqstr = contigs[key]  # retrieve the nucleotide sequence for each contig
        for k in range(2):  # each contig has two lists of ORFs: 1) forward direction; 2) reverse complement
            ind_count = 0  # a counter that counts the number of ORFs for each direction (to be used as an index)
            if k == 1:
                seqstr = revcomp(seqstr)  # if k==1, then the ORFs refer to the reverse complement of the nucleotide sequence
            for orf_ind in orf_ind_dict[key][k]:  # iterate over the ORF indices for each direction
                start, stop = orf_ind_dict[key][k][ind_count]  # identify the start/stop indices for each ORF
                orf_dict['>protein_%s' % count] = seqstr[start:stop]  # build an ORF nucleotide entry in the ORF dictionary using the start/stop indices and the contig sequence
                count += 1  # add 1 to the counter used to build the unique identifier for each protein
                ind_count += 1  # add 1 to the ORF counter for each direction

    # The following dictionary is used to translate codons into amino acids.
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

    protein_dict = {}  # a dictionary of translated protein sequences: protein ID (key), translated protein sequence (value)
    for prot_name, orf in orf_dict.items():  # iterate over the protein nucleotide ORF dictionary
        position = 0  # a counter to keep track of nucleotide position
        protein = ''
        while position < len(orf):  # keep translating codons until you reach a nucleotide position that is >= the length of the ORF
            codon = orf[position:position + 3]  # retrieve three nucleotides to define a codon
            protein += codons[codon]  # translate a codon into an amino acid and append it to the protein sequence
            position += 3  # add three to the nucleotide position and move to the next codon
        protein_dict[prot_name] = protein  # add the translated protein sequence to the protein dictionary

    protein_key_list = sorted(protein_dict)  # sort the keys of the protein dictionary; result: p_1, p_10, ... p_19, p_2, p_20, ... p_29, p_3, etc.

    protein_key_lengths = []  # a clever way to sort the protein keys properly is to use numpy.lexsort and the lengths of the protein IDs, i.e., then p_2 will come before p_10, etc.
    for id in protein_key_list:  # iterate over the protein IDs in 'protein_key_list'
        protein_key_lengths.append(len(id))  # append the length of each protein ID string to the protein key length list
    indices = np.lexsort((protein_key_list, protein_key_lengths))  # get a sorted list of protein ID indices so we can order them appropriately: p_1, ... p_9, p_10, ...

    index_counter = 0  # this will be the index of the protein ID in 'protein_key_list' (first sort... from the protein dictionary)
    file_seq = ''
    sorted_protein_keys = []  # this list will use 'indices' and 'protein_key_list' to properly sort the protein IDs, used to construct an annotation file with properly sorted protein IDs
    for index in indices:  # iterate over each index in 'indices'
        sorted_protein_keys.append(protein_key_list[index])  # append the appropriate protein ID to the sorted protein ID list, result p_1, p_2, p_3, ... p_10, ...
        if (index_counter + 1) == len(indices):  # this 'if' statement constructs the string of protein IDs and sequences which will be blasted against the database
            file_seq += protein_key_list[index] + '\n' + protein_dict[protein_key_list[index]]  # if the final protein sequence is being appended, do not add a newline character to the end of the file
        else:
            file_seq += protein_key_list[index] + '\n' + protein_dict[protein_key_list[index]] + '\n'  # if the protein sequence isn't the final one, it must have a newline character at the end
        index_counter += 1  # add one to the index counter to determine when the last protein is being appended to the string

    best_hits = blast_against_db(file_seq)  # pass the 'file_seq' protein IDs and sequences string to the blast function

    protein_labels = {}  # a protein label dictionary used for annotation
    for prot_id, record in best_hits.items():  # iterate over the protein IDs and blast records in the blast results
        match = re.search(' (.+) OS=', record.stitle)  # retrieve the protein label from the blast results
        if float(record.pident) > 50.0:
            protein_label = match.group(1)  # the protein annotation is capture group 1
        else:
            protein_label = 'hypothetical protein'  # if the protein percent identity is < 50, label it as a 'hypothetical protein'
        protein_labels['>' + prot_id] = protein_label  # add the protein label to the protein label dictionary

    annotate_proteins(protein_dict, protein_labels, sorted_protein_keys)  # pass the sequences, labels, and ordered IDs to the function that writes the annotation file
