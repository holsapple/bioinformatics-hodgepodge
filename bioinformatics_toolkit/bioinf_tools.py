"""
This is a beginner's bioinformatics module.

The module contains these functions: gen_probability()
                                     mkseq()
                                     seq2str()
                                     str2seq()
                                     revcomp()
                                     find_start()
                                     find_stop()
                                     orf_finder()
                                     filter_int()
                                     main

Detailed descriptions of each function and test are located in the corresponding docstrings.

AUTHOR: Raymond Holsapple
LAST UPDATE: 10/12/2016
"""


import numpy as np
import re


def gen_probability(GCcontent=0.5, GCskew=0, ATskew=0):
    """
    This function generates probabilities for the A, T, G, and C nucleotides to be used for generating a DNA sequence in the mkseq function.

    :param GCcontent: (float) The desired GC content for the DNA sequence to be generated; default set to 0.5.
    :param GCskew: (float) The desired GC skewness for the DNA sequence to be generated; default set to zero.
    :param ATskew: (float) The desired AT skewness for the DNA sequence to be generated; default set to zero.
    :return nuc_probs: (numpy array of floats) The calculated probabilities for A, T, G, and C to be used in the mkseq function.
    """

    G = GCcontent * (1 - (1 - GCskew)/2)  # Probability of generating a G when calling np.random.choice() in the mkseq function.
    C = GCcontent - G  # Probability of generating a C when calling np.random.choice() in the mkseq function.
    A = (1 - GCcontent) * (1 - (1 - ATskew)/2)  # Probability of generating an A when calling np.random.choice() in the mkseq function.
    T = 1 - G - C - A  # Probability of generating a T when calling np.random.choice() in the mkseq function.

    nuc_probs = np.array([A, T, G, C])  # Create the numpy array of probabilities.

    return nuc_probs  # Return the nucleotide probabilities that will be used in the mkseq function.


def mkseq(seq_size, GCcontent=0.5, GCskew=0, ATskew=0):
    """
    This function generates a random DNA sequence.

    :param seq_size: (integer) Length of the random DNA sequence to be generated.
    :param GCcontent: (float) The desired GC content for the DNA sequence to be generated; default set to 0.5.
    :param GCskew: (float) The desired GC skewness for the DNA sequence to be generated; default set to zero.
    :param ATskew: (float) The desired AT skewness for the DNA sequence to be generated; default set to zero.
    :return DNA_seq: (numpy array of strings) This is the desired random DNA sequence the function is meant to generate.
    """

    nuc_alpha = np.array(['A', 'T', 'G', 'C'])  # Create the nucleotide alphabet that np.random.choice will use.
    nuc_probs = gen_probability(GCcontent, GCskew, ATskew)  # Call the helper function gen_probability to create the nucleotide probabilities.
    DNA_seq = np.random.choice(nuc_alpha, seq_size, True, nuc_probs)  # Here we use the input parameters to generate the random DNA sequence of desired length.

    return DNA_seq  # Return the random DNA sequence.


def seq2str(DNA_seq_array):
    """
    This function converts a numpy array DNA sequence to a string DNA sequence.

    :param DNA_seq_array: (numpy array of strings) DNA sequence in numpy array format.
    :return DNA_seq_str: (string) DNA sequence in string format.
    """

    DNA_seq_str = ''  # Initialize the DNA sequence string

    for n in np.nditer(DNA_seq_array):  # Here we iterate over the elements in the numpy array DNA_seq_array so we can append them to a string
        DNA_seq_str += str(n)  # Convert the object n to a string and then append it DNA_seq_array

    return DNA_seq_str  # Return the DNA sequence in string format.


def str2seq(DNA_seq_str):
    """
    This function converts a string DNA sequence to a numpy array DNA sequence.

    :param DNA_seq_str: (string) DNA sequnce in string format.
    :return DNA_seq_array: (numpy array of strings) DNA sequence in numpy array format.
    """

    temp_list = []  # Initialize a temporary list that will be converted into a numpy array object

    for c in DNA_seq_str:  # Here we iterate over the characters in the string DNA_seq_str so we can append them to our temporary list
        temp_list.append(c)  # Append each character to the list temp_list

    DNA_seq_array = np.asarray(temp_list)  # Convert the list temp_list to a numpy array using numpy's asarray method

    return DNA_seq_array  # Return the DNA sequence in numpy array format.

def revcomp(DNA_seq):
    """
    This function creates the reverse complement of a DNA sequence.

    :param DNA_seq: (numpy array of strings) This is the DNA sequence for which we need to compute the reverse complement sequence.
    :return DNA_seq_revcomp: (numpy array of strings) This is the reverse complement of the input parameter DNA_seq.
    """

    DNA_seq_str = seq2str(DNA_seq)  # Convert the given DNA numpy array into a string
    DNA_seq_str_rev = DNA_seq_str[::-1]  # Take a negative slice of the DNA sequence string, which will reverse the string
    DNA_seq_str_revcomp = ''  # Initialize a new string which will become the reversed complement of the given DNA sequence string

    for nuc in DNA_seq_str_rev:  # In this for loop we iterate over characters in the string DNA_seq_str_rev so that we can build a complement sequence one element at a time
        if nuc == 'A':
            DNA_seq_str_revcomp += 'T'  # If the nucleotide is an A, change it to a T and append it to the string DNA_seq_str_revcomp
        elif nuc == 'T':
            DNA_seq_str_revcomp += 'A'  # If the nucleotide is a T, change it to an A and append it to the string DNA_seq_str_revcomp
        elif nuc == 'C':
            DNA_seq_str_revcomp += 'G'  # If the nucleotide is a C, change it to a G and append it to the string DNA_seq_str_revcomp
        else:
            DNA_seq_str_revcomp += 'C'  # If the nucleotide is a G, change it to a C and append it to the string DNA_seq_str_revcomp

    DNA_seq_revcomp = str2seq(DNA_seq_str_revcomp)  # Convert the reversed complement string that was just created into a numpy array object

    return DNA_seq_revcomp  # Return the reversed complement numpy array DNA sequence.


def find_start(DNA_seq):
    """
    This function determines the beginning indices for all start codons in the sequence DNA_seq

    :param DNA_seq: (numpy array of strings) This is the DNA sequence we use to search for start codons in.
    :return start_indices: (list of integers) These are the indices for the first nucleotide of a start codon within the sequence DNA_seq.
    """

    start_indices = []  # Initialize the list of beginning indices for the start codons
    A_starts = np.where(DNA_seq == 'A')  # Determine which indices of DNA_seq are an A, since the start codon is ATG
    last_possible_index = len(DNA_seq) - 3  # This is the last possible index where my algorithm could locate the first nucleotide of a 3-nucleotide start codon

    for ind in A_starts[0]:  # Here I iterate through the indices in the list of indices indicating where DNA_seq has an A nucleotide
        if ind <= last_possible_index:  # This line prevents the algorithm from looking for a start codon beginning in either of the last two elements of the DNA sequence
            if DNA_seq[ind + 1] != 'T' or DNA_seq[ind + 2] != 'G':
                continue  # If the A nucleotide isn't immediately followed by a T or doesn't have a G two bases further into the sequence, then it is not a start codon so move on in the loop
            else:
                start_indices.append(ind)  # The negation of the 'if' statement implies that the A is immediately followed by a T which is immediately followed by a G (De Morgan's Laws)

    return start_indices  # Return the list of beginning indices for all start codons in DNA_seq.


def find_stop(DNA_seq):
    """
    This function determines the beginning indices for all stop codons in the sequence DNA_seq

    :param DNA_seq: (numpy array of strings) This is the DNA sequence we use to search for stop codons in.
    :return stop_indices: (list of integers) These are the indices for the first nucleotide of a stop codon within the sequence DNA_seq.
    """

    stop_indices = []  # Initialize the list of beginning indices for the stop codons
    T_stops = np.where(DNA_seq == 'T')  # Determine which indices of DNA_seq are a T, since all stop codons start with a T
    last_possible_index = len(DNA_seq) - 3  # This is the last possible index where my algorithm could locate the first nucleotide of a 3-nucleotide stop codon

    for ind in T_stops[0]:  # Here I iterate through the indices in the list of indices indicating where DNA_seq has an T nucleotide
        if ind <= last_possible_index:  # This line prevents the algorithm from looking for a stop codon beginning in either of the last two elements of the DNA sequence
            final2 = DNA_seq[ind + 1] + DNA_seq[ind + 2]  # For each index in T_stops, I build a 2-letter string consisting of the two nucleotides in DNA_seq that follow that particular T
            if final2 == 'AA' or final2 == 'AG' or final2 == 'GA':
                stop_indices.append(ind)  # If the 2 nucleotides following a T are any of the 3 choices in the 'if' statement above, then the corresponding index points to a stop codon

    return stop_indices  # Return the list of beginning indices for all stop codons in DNA_seq.


def orf_finder(DNA_seq):
    """
    This function takes a numpy array DNA sequnce and finds all of the open reading frames within this sequence.

    Note: A DNA sequence should have 6 possible reading frames: three in the forward direction of DNA_seq (see below) and three for the reverse complement of DNA_seq. The zero-th
          (indicated by 0th) reading frame assumes codons begin with the first nucleotide of the given DNA fragment. The 1st reading frame assumes codons begin with the second
          nucleotide of the given DNA fragment. The 2nd reading frame assumes codons begin with the third nucleotide of the given DNA fragment. The same is true for the reverse
          complement of the given DNA sequence.

    :param DNA_seq: (numpy array of strings) This is the DNA sequence we search for open reading frames (and also its reverse complement).
    :return [orf_for, orf_bak]: (embedded lists of tuples of integers) This function will return two lists: orf_for and orf_bak.
            orf_for is itself a list of three distinct lists: one for each each of the three possible reading frames (0th, 1st, 2nd) of DNA_seq.
            Each of these sub-lists within orf_for is a list of tuples: consisting of a start index (i.e., index of first nucleotide of start codon) and
            stop index (i.e., index of third nucleotide of a corresponding stop codon) for one open reading frame within a specific reading frame.
            orf_bak is just like orf_for, except it applies to the reverse complement of DNA_seq, not DNA_seq itself.
    """

    # -------------------------------------------------------------------------------------------------------------------------------------------------------- #
    # The following two lines were used for debugging the logic of the embedded for loops while writing the function... prior to testing the function in main.
    # This is much simpler than constructing a fake DNA sequence having these start and stop codon indices just to debug the logic in the loops.
    # Because once you have these indices it is simple to construct the lists representing exactly what the function should return.
    # start_indices = [11, 41, 45, 110, 120, 292, 309, 369, 381, 481, 547, 666]
    # stop_indices = [5, 77, 98, 165, 170, 322, 400, 408, 601, 702]
    # -------------------------------------------------------------------------------------------------------------------------------------------------------- #

    start_indices = find_start(DNA_seq)  # Determine the beginning indices of all start codons in DNA_seq
    stop_indices = find_stop(DNA_seq)  # Determine the beginning indices of all stop codons in DNA_seq
    orf_for = [[], [], []]  # Initialize the list (of lists) of all open reading frames (for each of the 3 possible reading frames) in DNA_seq

    for start in start_indices:  # Iterate over all the beginning start codon indices
        for stop in stop_indices:  # Iterate over all the beginning stop codon indices
            if (stop <= start) or (((stop - start) % 3) != 0):  # This 'if' statement determines if the stop index is less than the start index of if the two indices are in different reading frames
                continue  # If either of the conditions mentioned above are true then we continue to the next beginning stop codon index in stop_indices
            read_frame = start % 3  # If neither of the above conditions are true then this might be the stop codon we are looking for, so we compute the reading frame they are in (0, 1, or 2)
            if orf_for[read_frame]:  # If the read frame for these start/stop codons already has an orf in it, then the current start/stop codon pair might not be a valid orf... let's find out
                prev_read_frame_pair = orf_for[read_frame][-1]  # Step 1 of this task: Since the read frame has 'something' in it, retrieve the most recent orf added, which is the only one of concern
                prev_stop = prev_read_frame_pair[1]  # Step 2 of this task: Use this most recent orf and make note of what the beginning index is for the stop codon
                if start > prev_stop:  # Determine if this current start codon comes after the most recent stop codon that was already added to the orf
                    orf_for[read_frame].append((start, stop + 2))  # If condition is true, then this is a new start/stop codon pair that we want to append to the list for the current reading frame
                break  # If you've added a start/stop codon pair to the list, then you are done with this start codon, so break out of the stop codon loop and iterate to the next start index
            else:
                orf_for[read_frame].append((start, stop + 2))  # If the reading frame for the current start/stop codons is empty, then they are definitely a pair that we want to put in the list
                break  # Once you add a start/stop codon pair to the current reading frame codon list, then you have no more need to check other stop codons... move to the next start index

    # --------------------------------------------------------------------------------------------------------------------- #
    # The code below is uncommented because the logic is identical to the logic above in this function. The only difference
    # is that the code below applies to the reverse complement of DNA_seq, whereas above it applies to DNA_seq itself.
    # --------------------------------------------------------------------------------------------------------------------- #
    DNA_seq_revcomp = revcomp(DNA_seq)
    start_indices = find_start(DNA_seq_revcomp)
    stop_indices = find_stop(DNA_seq_revcomp)
    orf_bak = [[], [], []]

    for start in start_indices:
        for stop in stop_indices:
            if (stop <= start) or (((stop - start) % 3) != 0):
                continue
            read_frame = start % 3
            if orf_bak[read_frame]:
                prev_read_frame_pair = orf_bak[read_frame][-1]
                prev_stop = prev_read_frame_pair[1]
                if start > prev_stop:
                    orf_bak[read_frame].append((start, stop + 2))
                break
            else:
                orf_bak[read_frame].append((start, stop + 2))
                break

    return [orf_for, orf_bak]  # Return the lists of open reading frames.


def filter_int(DNA_seq):
    """
    This function takes a DNA sequence and and filters out introns according to the following definition:
    An intron is any subsequence that begins with 'GT', ends with 'AG' and has an 'A somewhere in between.

    Note: The pattern search is not greedy, because we want to match all of the shortest such patterns that
          satisfy the above intron definition.

    :param DNA_seq: (numpy array of strings) This is the DNA sequence we want to filter the introns from.
    :return filtered_DNA_seqs: (list of numpy arrays of strings) The first element of this list is the intron-free version
            of DNA_seq. The second element of this list is the intron-free version of DNA_seq's reverse complement.
    """

    DNA_seq_str = seq2str(DNA_seq)  # Convert the numpy array version of the sequence into a string so we can use the re module
    filtered_DNA_seq_str = re.sub('GT.*?A.*?AG', '', DNA_seq_str)  # This is where we find our regex pattern in the sequence and then cut it out
    filtered_DNA_seq = str2seq(filtered_DNA_seq_str)  # Here I convert the intron-filtered string back into a numpy array

    DNA_seq_revcomp = revcomp(DNA_seq)  # Here we compute the reverse complement of the DNA sequence, so we can remove intons from it as well
    DNA_seq_revcomp_str = seq2str(DNA_seq_revcomp)  # Convert the numpy array version of the reverse complement into a string so we can use the re module
    filtered_DNA_seq_revcomp_str = re.sub('GT.*?A.*?AG', '', DNA_seq_revcomp_str)  # This is where we find our regex pattern in the reverse complement and then cut it out
    filtered_DNA_seq_revcomp = str2seq(filtered_DNA_seq_revcomp_str)  # Here I convert the intron-filtered string back into a numpy array

    return [filtered_DNA_seq, filtered_DNA_seq_revcomp]  # Return the intron-filtered numpy arrays for the original DNA sequence and its reverse complement.


if __name__ == '__main__':
    """
    This is the main function where I test each of the functions defined above. The descriptions for the individual tests are contained in each function's docstring.
    """

    # -------------------- Test 1: test gen_probability --------------------
    """
    This is a basic test of the gen_probability function. The GC & AT skews are set to their default values of zero, and the GC content is set to its default value of 0.5.
    The desired output is [A, T, G, C] = [0.25, 0.25, 0.25, 0.25].
    """
    print()
    print('---------- Test 1: Testing gen_probability ----------')
    print()
    nuc_probs = gen_probability()  # Call gen_probability with the desired parameters.
    print(nuc_probs)
    eps = 1e-10  # Computed values should vary no more than this amount from the true (theoretical) value.
    # This if/else statement determines if the computed A, T, G, and C values lie within the desired variance. If yes, the test is successful. If no, the test fails.
    if (abs(nuc_probs[0] - 0.25) < eps) and (abs(nuc_probs[1] - 0.25) < eps) and (abs(nuc_probs[2] - 0.25) < eps) and (abs(nuc_probs[3] - 0.25) < eps):
        print('Successful Test')
    else:
        print('Unsuccessful Test')
    print()
    print()

    # -------------------- Test 2: test gen_probability --------------------
    """
    In this test of the gen_probability function, the AT skew is set to the default value of zero, but the GC skew is set to 0.05.
    The GC content is set to 0.4, so the desired output is [A, T, G, C] = [0.3, 0.3, 0.21, 0.19].
    """
    print()
    print('---------- Test 2: Testing gen_probability ----------')
    print()
    nuc_probs = gen_probability(0.4, 0.05)  # Call gen_probability with the desired parameters.
    print(nuc_probs)
    eps = 1e-10  # Computed values should vary no more than this amount from the true (theoretical) value.
    # This if/else statement determines if the computed A, T, G, and C values lie within the desired variance. If yes, the test is successful. If no, the test fails.
    if (abs(nuc_probs[0] - 0.3) < eps) and (abs(nuc_probs[1] - 0.3) < eps) and (abs(nuc_probs[2] - 0.21) < eps) and (abs(nuc_probs[3] - 0.19) < eps):
        print('Successful Test')
    else:
        print('Unsuccessful Test')
    print()
    print()

    # -------------------- Test 3: test gen_probability --------------------
    """
    In this test of the gen_probability function, the AT skew is set to 0.1, and the GC skew is set to 0.05.
    The GC content is set to 0.35, so the desired output is [A, T, G, C] = [0.3575, 0.2925, 0.18375, 0.16625].
    """
    print()
    print('---------- Test 3: Testing gen_probability ----------')
    print()
    nuc_probs = gen_probability(0.35, 0.05, 0.1)  # Call gen_probability with the desired parameters.
    print(nuc_probs)
    eps = 1e-10  # Computed values should vary no more than this amount from the true (theoretical) value.
    # This if/else statement determines if the computed A, T, G, and C values lie within the desired variance. If yes, the test is successful. If no, the test fails.
    if (abs(nuc_probs[0] - 0.3575) < eps) and (abs(nuc_probs[1] - 0.2925) < eps) and (abs(nuc_probs[2] - 0.18375) < eps) and (abs(nuc_probs[3] - 0.16625) < eps):
        print('Successful Test')
    else:
        print('Unsuccessful Test')
    print()
    print()

    # -------------------- Test 4: test mkseq --------------------
    """
    In this test of the mkseq function, I will generate a random DNA sequence of length 20, with GC content set to its default value of 0.5.
    I leave the AT & GC skews set to their default values of zero.
    """
    print()
    print('---------- Test 4: Testing mkseq ----------')
    print()
    DNA_seq = mkseq(20)  # Call mkseq to generate the random DNA sequence.
    print(DNA_seq)

    # Now I want to check the GC content, GC skew, and AT skew of the randomly generated sequence and compare it to the desired values of 0.5, 0, 0 (respectively).
    # As the size of the sequence gets larger these computed differences should get closer to zero.
    A = len(np.where(DNA_seq == 'A')[0])  # Number of A's generated.
    T = len(np.where(DNA_seq == 'T')[0])  # Number of T's generated.
    G = len(np.where(DNA_seq == 'G')[0])  # Number of G's generated.
    C = len(np.where(DNA_seq == 'C')[0])  # Number of C's generated.
    GCcontent = (G + C) / (A + T + G + C)   # GC content of the generated sequence
    GCcontent_change = round((abs(0.5 - GCcontent) / 0.5) * 100)  # Compute the percentage deviation of GC content.
    GCskew = round((G - C) / (G + C), 2)  # GC skew of the generated sequence
    ATskew = round((A - T) / (A + T), 2)  # AT skew of the generated sequence

    print('The number of A nucleotides generated is', A)
    print('The number of T nucleotides generated is', T)
    print('The number of G nucleotides generated is', G)
    print('The number of C nucleotides generated is', C)
    print('The GC content of this DNA sequence is', GCcontent, 'which is a', GCcontent_change, '% deviation from the desired GC content.')
    print('The GC skew of this DNA sequence is', GCskew, 'and the desired GC skew was 0.')
    print('The AT skew of this DNA sequence is', ATskew, 'and the desired AT skew was 0.')
    print()
    print()

    # -------------------- Test 5: test mkseq --------------------
    """
    In this test of the mkseq function, I will generate a random DNA sequence of length 100, with GC content set to 0.4.
    The AT skew is set to 0.1 & the GC skew is set to 0.05.
    """
    print()
    print('---------- Test 5: Testing mkseq ----------')
    print()
    DNA_seq = mkseq(100, 0.4, 0.05, 0.1)  # Call mkseq to generate the random DNA sequence.
    print(DNA_seq)

    # Now I want to check the GC content, GC skew, and AT skew of the randomly generated sequence and compare it to the desired values of 0.5, 0, 0 (respectively).
    # As the size of the sequence gets larger these computed differences should get closer to zero.
    A = len(np.where(DNA_seq == 'A')[0])  # Number of A's generated.
    T = len(np.where(DNA_seq == 'T')[0])  # Number of T's generated.
    G = len(np.where(DNA_seq == 'G')[0])  # Number of G's generated.
    C = len(np.where(DNA_seq == 'C')[0])  # Number of C's generated.
    GCcontent = (G + C) / (A + T + G + C)   # GC content of the generated sequence
    GCcontent_change = round((abs(0.4 - GCcontent) / 0.4) * 100)  # Compute the percentage deviation of GC content.
    GCskew = round((G - C) / (G + C), 2)  # GC skew of the generated sequence
    ATskew = round((A - T) / (A + T), 2)  # AT skew of the generated sequence

    print('The number of A nucleotides generated is', A)
    print('The number of T nucleotides generated is', T)
    print('The number of G nucleotides generated is', G)
    print('The number of C nucleotides generated is', C)
    print('The GC content of this DNA sequence is', GCcontent, 'which is a', GCcontent_change, '% deviation from the desired GC content.')
    print('The GC skew of this DNA sequence is', GCskew, 'and the desired GC skew was 0.05.')
    print('The AT skew of this DNA sequence is', ATskew, 'and the desired AT skew was 0.1.')
    print()
    print()

    # -------------------- Test 6: test mkseq --------------------
    """
    In this test of the mkseq function, I will generate a random DNA sequence of length 500, with GC content set to 0.25.
    The AT skew is set to 0.01 & the GC skew is set to 0.1.
    """
    print()
    print('---------- Test 6: Testing mkseq ----------')
    print()
    DNA_seq = mkseq(500, 0.25, 0.1, 0.01)  # Call mkseq to generate the random DNA sequence.
    print(DNA_seq)

    # Now I want to check the GC content, GC skew, and AT skew of the randomly generated sequence and compare it to the desired values of 0.5, 0, 0 (respectively).
    # As the size of the sequence gets larger these computed differences should get closer to zero.
    A = len(np.where(DNA_seq == 'A')[0])  # Number of A's generated.
    T = len(np.where(DNA_seq == 'T')[0])  # Number of T's generated.
    G = len(np.where(DNA_seq == 'G')[0])  # Number of G's generated.
    C = len(np.where(DNA_seq == 'C')[0])  # Number of C's generated.
    GCcontent = (G + C) / (A + T + G + C)  # GC content of the generated sequence
    GCcontent_change = round((abs(0.25 - GCcontent) / 0.25) * 100, 1)  # Compute the percentage deviation of GC content.
    GCskew = round((G - C) / (G + C), 2)  # GC skew of the generated sequence
    ATskew = round((A - T) / (A + T), 2)  # AT skew of the generated sequence

    print('The number of A nucleotides generated is', A)
    print('The number of T nucleotides generated is', T)
    print('The number of G nucleotides generated is', G)
    print('The number of C nucleotides generated is', C)
    print('The GC content of this DNA sequence is', GCcontent, 'which is a', GCcontent_change, '% deviation from the desired GC content.')
    print('The GC skew of this DNA sequence is', GCskew, 'and the desired GC skew was 0.1.')
    print('The AT skew of this DNA sequence is', ATskew, 'and the desired AT skew was 0.01.')
    print()
    print()

    # -------------------- Test 7: test seq2str --------------------
    """
    In this test of the seq2str function, I will generate a random DNA sequence of length 15, with GC content set to the default value of 0.5.
    The AT & GC skews are set to their default values of zero. Then the numpy array DNA sequence is converted to a string and both are printed for comparison.
    """
    print()
    print('---------- Test 7: Testing seq2str ----------')
    print()
    DNA_seq = mkseq(15)  # Call mkseq to generate the random DNA sequence in numpy array format.
    DNA_seq_str = seq2str(DNA_seq)  # Convert the numpy array DNA_seq to a string.
    print('The starting array version is', DNA_seq)
    print('Just to verify, the array type is', type(DNA_seq))
    print('The coverted string version is', DNA_seq_str)
    print('Just to verify, the string type is', type(DNA_seq_str))
    print()
    print()

    # -------------------- Test 8: test str2seq --------------------
    """
    In this test of the str2seq function, I will first generate a random DNA sequence of length 15 in numpy array format since I have already built the tools to do so.
    Then I will use the seq2str function to convert that numpy array into a string. That function has already been defined and tested. Once I have the DNA sequence in
    string format I will pass it to the str2seq function and convert it back into a numpy array. The random DNA sequence will have the GC content set to the default
    value of 0.5, and the AT & GC skews set to their default values of zero. Once all of this is completed, I print the string version of the generated sequence,
    and then the result of passing that string to the str2seq function.
    """
    print()
    print('---------- Test 8: Testing str2seq ----------')
    print()
    DNA_seq = mkseq(15)  # Call mkseq to generate the random DNA sequence in numpy array format.
    DNA_seq_str = seq2str(DNA_seq)  # Convert the numpy array DNA_seq to a string.
    DNA_seq_array = str2seq(DNA_seq_str)  # Convert the string DNA_seq_str back into a numpy array.
    print('The starting string version is', DNA_seq_str)
    print('Just to verify, the string type is', type(DNA_seq_str))
    print('The converted array version is', DNA_seq_array)
    print('Just to verify, the array type is', type(DNA_seq_array))
    print()
    print()

    # -------------------- Test 9: test revcomp --------------------
    """
    In this test of the revcomp function, I will generate a random DNA sequence of length 15, with GC content set to the default value of 0.5.
    The AT & GC skews are set to their default values of zero. Then the numpy array DNA sequence is is passed to the revcomp function where its
    reverse complement is created and returned. I then print both the orignial numpy array DNA sequence and its reverse complement for comparison.
    """
    print()
    print('---------- Test 9: Testing revcomp ----------')
    print()
    DNA_seq = mkseq(15)  # Call mkseq to generate the random DNA sequence in numpy array format.
    DNA_seq_revcomp = revcomp(DNA_seq)  # Create the reverse complement of the numpy array DNA_seq.
    print('The original sequence is', DNA_seq)
    print('The reverse complement is', DNA_seq_revcomp)
    print()
    print()

    # -------------------- Test 10: test find_start --------------------
    """
    In this test of the find_start function, I will generate a test DNA sequence of length 50. I populate the test sequence with zero start codons.
    The test sequence is passed to find_start and then I print the returned list of start codon beginning indices and compare this list to the known
    beginning indices of the start codons in the test sequence. The function should return an empty list [].
    """
    print()
    print('---------- Test 10: Testing find_start ----------')
    print()
    DNA_seq = np.array(['A', 'C', 'A', 'T', 'C', 'C', 'T', 'A', 'A', 'T',  # Here I build the test DNA sequence in numpy array format.
                        'C', 'G', 'C', 'A', 'A', 'T', 'T', 'G', 'T', 'G',  # There are no start codons in this numpy array.
                        'T', 'C', 'T', 'T', 'A', 'G', 'A', 'A', 'G', 'T',
                        'T', 'A', 'T', 'T', 'C', 'C', 'A', 'T', 'C', 'A',
                        'T', 'T', 'G', 'A', 'C', 'T', 'G', 'A', 'C', 'A'])
    start_indices = find_start(DNA_seq)  # Pass DNA_seq to the find_start function to find the beginning indices of the start codons.
    truth_indices = []  # These are the truth values for the desired indices which were determined by construction of the numpy array DNA_seq.
    same_flag = 1  # A flag which will be used to indicate whether or not start_indices and truth_indices are identical: the default is set to true.

    for i in start_indices:  # This 'for' loop iterates over the elements in start_indices to determine if the list is the same as truth_indices.
        if i not in truth_indices:
            same_flag = 0  # If an index from start_indices is not in truth_indices, set same_flag equal to false... indicating that the lists are not identical.

    if same_flag == 1:
        print('Successful Test: find_start returned the correct list')  # If same_flag is true, print a successful test message.
    else:
        print('Unsuccessful Test: find_start returned the incorrect list')  # If same_flag is false, print an unsuccessful test message.
    print('The true start indices are', truth_indices)
    print('The find_start function returned these start indices', start_indices)
    print()
    print()

    # -------------------- Test 11: test find_start --------------------
    """
    In this test of the find_start function, I will generate a test DNA sequence of length 100. I populate the test sequence with multiple start codons.
    The test sequence is passed to find_start and then I print the returned list of start codon beginning indices and compare this list to the known
    beginning indices of the start codons in the test sequence. The function should return [2, 17, 26, 31, 44].
    """
    print()
    print('---------- Test 11: Testing find_start ----------')
    print()
    DNA_seq = np.array(['A', 'C', 'A', 'T', 'G', 'C', 'T', 'A', 'A', 'T',  # Here I build the test DNA sequence in numpy array format.
                        'C', 'G', 'C', 'A', 'A', 'T', 'T', 'A', 'T', 'G',  # The beginning indices of the start codons are:
                        'T', 'C', 'T', 'T', 'A', 'G', 'A', 'T', 'G', 'T',  # 2, 17, 26, 31, 44, 52, 67, 76, 81, and 94.
                        'T', 'A', 'T', 'G', 'C', 'C', 'A', 'T', 'C', 'A',
                        'T', 'T', 'G', 'A', 'A', 'T', 'G', 'A', 'C', 'A',
                        'A', 'C', 'A', 'T', 'G', 'C', 'T', 'A', 'A', 'T',
                        'C', 'G', 'C', 'A', 'A', 'T', 'T', 'A', 'T', 'G',
                        'T', 'C', 'T', 'T', 'A', 'G', 'A', 'T', 'G', 'T',
                        'T', 'A', 'T', 'G', 'C', 'C', 'A', 'T', 'C', 'A',
                        'T', 'T', 'G', 'A', 'A', 'T', 'G', 'A', 'C', 'A'])
    start_indices = find_start(DNA_seq)  # Pass DNA_seq to the find_start function to find the beginning indices of the start codons.
    truth_indices = [2, 17, 26, 31, 44, 52, 67, 76, 81, 94]  # These are the truth values for the desired indices which were determined by construction of the numpy array DNA_seq.
    same_flag = 1  # A flag which will be used to indicate whether or not start_indices and truth_indices are identical: the default is set to true.

    for i in start_indices:  # This 'for' loop iterates over the elements in start_indices to determine if the list is the same as truth_indices.
        if i not in truth_indices:
            same_flag = 0  # If an index from start_indices is not in truth_indices, set same_flag equal to false... indicating that the lists are not identical.

    if same_flag == 1:
        print('Successful Test: find_start returned the correct list')  # If same_flag is true, print a successful test message.
    else:
        print('Unsuccessful Test: find_start returned the incorrect list')  # If same_flag is false, print an unsuccessful test message.
    print('The true start indices are', truth_indices)
    print('The find_start function returned these start indices', start_indices)
    print()
    print()

    # -------------------- Test 12: test find_start --------------------
    """
    In this test of the find_start function, I will generate a random DNA sequence of length 1000, with GC content set to the default value of 0.5.
    The AT & GC skews are set to their default values of zero. The DNA sequence is passed to find_start and then I print the returned list of start
    codon beginning indices. Since this is a randomly generated sequence, I do not know the true locations of the start codons.
    """
    print()
    print('---------- Test 12: Testing find_start ----------')
    print()
    DNA_seq = mkseq(1000)  # Here we generate a random DNA sequence using our mkseq function.
    start_indices = find_start(DNA_seq)  # Pass DNA_seq to the find_start function to find the beginning indices of the start codons.
    start_codons = []  # Initialize a start codon list to verify all returned indices point to a start codon.

    for ind in start_indices:  # Iterate over the indices in start_indices to build the start_codons list.
        codon_str = DNA_seq[ind] + DNA_seq[ind + 1] + DNA_seq[ind + 2]  # For each returned beginning start codon index, create the 3-letter string that it points to.
        start_codons.append(codon_str)  # Append the 3-letter word just created to the start_codons list.

    print('The find_start function returned these start indices', start_indices)  # Here we simply print the list of beginning indices of all start codons from DNA_seq.
    print('The list of codons from the random DNA sequence having beginning indices as indicated in the list above', start_codons)  # All of these should be 'ATG'
    print()
    print()

    # -------------------- Test 13: test find_stop --------------------
    """
    In this test of the find_stop function, I will generate a test DNA sequence of length 50. I populate the test sequence with zero stop codons.
    The test sequence is passed to find_stop and then I print the returned list of stop codon beginning indices and compare this list to the known
    beginning indices of the stop codons in the test sequence. The function should return an empty list [].
    """
    print()
    print('---------- Test 13: Testing find_stop ----------')
    print()
    DNA_seq = np.array(['A', 'C', 'A', 'T', 'C', 'C', 'T', 'T', 'A', 'T',  # Here I build the test DNA sequence in numpy array format.
                        'C', 'G', 'C', 'A', 'A', 'T', 'T', 'G', 'T', 'G',  # There are no start codons in this numpy array.
                        'T', 'C', 'T', 'T', 'A', 'C', 'A', 'A', 'G', 'T',
                        'T', 'A', 'T', 'T', 'C', 'C', 'A', 'T', 'C', 'A',
                        'T', 'T', 'G', 'G', 'C', 'A', 'G', 'A', 'C', 'A'])
    stop_indices = find_stop(DNA_seq)  # Pass DNA_seq to the find_stop function to find the beginning indices of the stop codons.
    truth_indices = []  # These are the truth values for the desired indices which were determined by construction of the numpy array DNA_seq.
    same_flag = 1  # A flag which will be used to indicate whether or not stop_indices and truth_indices are identical: the default is set to true.

    for i in stop_indices:  # This 'for' loop iterates over the elements in stop_indices to determine if the list is the same as truth_indices.
        if i not in truth_indices:
            same_flag = 0  # If an index from stop_indices is not in truth_indices, set same_flag equal to false... indicating that the lists are not identical.

    if same_flag == 1:
        print('Successful Test: find_stop returned the correct list')  # If same_flag is true, print a successful test message.
    else:
        print('Unsuccessful Test: find_stop returned the incorrect list')  # If same_flag is false, print an unsuccessful test message.
    print('The true stop indices are', truth_indices)
    print('The find_stop function returned these stop indices', stop_indices)
    print()
    print()

    # -------------------- Test 14: test find_stop --------------------
    """
    In this test of the find_stop function, I will generate a test DNA sequence of length 100. I populate the test sequence with multiple stop codons.
    The test sequence is passed to find_stop and then I print the returned list of stop codon beginning indices and compare this list to the known
    beginning indices of the stop codons in the test sequence. The function should return [2, 26, 44, 55, 70, 79, 89, 97].
    """
    print()
    print('---------- Test 14: Testing find_stop ----------')
    print()
    DNA_seq = np.array(['A', 'C', 'T', 'A', 'A', 'C', 'T', 'A', 'C', 'T',  # Here I build the test DNA sequence in numpy array format.
                        'C', 'G', 'C', 'A', 'A', 'T', 'T', 'A', 'T', 'G',  # The beginning indices of the stop codons are:
                        'T', 'C', 'T', 'C', 'A', 'G', 'T', 'A', 'G', 'T',  # 2, 26, 44, 55, 70, 79, 89, and 97.
                        'T', 'A', 'T', 'G', 'C', 'C', 'A', 'T', 'C', 'A',
                        'T', 'C', 'G', 'A', 'T', 'G', 'A', 'A', 'C', 'A',
                        'A', 'C', 'A', 'T', 'G', 'T', 'A', 'A', 'A', 'T',
                        'C', 'G', 'C', 'A', 'A', 'T', 'T', 'A', 'C', 'G',
                        'T', 'A', 'G', 'T', 'C', 'G', 'A', 'T', 'G', 'T',
                        'G', 'A', 'T', 'G', 'C', 'C', 'A', 'T', 'C', 'T',
                        'A', 'A', 'G', 'A', 'A', 'T', 'T', 'T', 'A', 'G'])
    stop_indices = find_stop(DNA_seq)  # Pass DNA_seq to the find_stop function to find the beginning indices of the stop codons.
    truth_indices = [2, 26, 44, 55, 70, 79, 89, 97]  # These are the truth values for the desired indices which were determined by construction of the numpy array DNA_seq.
    same_flag = 1  # A flag which will be used to indicate whether or not stop_indices and truth_indices are identical: the default is set to true.

    for i in stop_indices:  # This 'for' loop iterates over the elements in stop_indices to determine if the list is the same as truth_indices.
        if i not in truth_indices:
            same_flag = 0  # If an index from stop_indices is not in truth_indices, set same_flag equal to false... indicating that the lists are not identical.

    if same_flag == 1:
        print('Successful Test: find_stop returned the correct list')  # If same_flag is true, print a successful test message.
    else:
        print('Unsuccessful Test: find_stop returned the incorrect list')  # If same_flag is false, print an unsuccessful test message.
    print('The true stop indices are', truth_indices)
    print('The find_stop function returned these stop indices', stop_indices)
    print()
    print()

    # -------------------- Test 15: test find_stop --------------------
    """
    In this test of the find_stop function, I will generate a random DNA sequence of length 1000, with GC content set to the default value of 0.5.
    The AT & GC skews are set to their default values of zero. The DNA sequence is passed to find_stop and then I print the returned list of stop
    codon beginning indices. Since this is a randomly generated sequence, I do not know the true locations of the stop codons.
    """
    print()
    print('---------- Test 15: Testing find_stop ----------')
    print()
    DNA_seq = mkseq(1000)  # Here we generate a random DNA sequence using our mkseq function.
    stop_indices = find_stop(DNA_seq)  # Pass DNA_seq to the find_stop function to find the beginning indices of the stop codons.
    stop_codons = []  # Initialize a stop codon list to verify all returned indices point to a stop codon.

    for ind in stop_indices:  # Iterate over the indices in stop_indices to build the stop_codons list.
        codon_str = DNA_seq[ind] + DNA_seq[ind + 1] + DNA_seq[ind + 2]  # For each returned beginning stop codon index, create the 3-letter string that it points to.
        stop_codons.append(codon_str)  # Append the 3-letter word just created to the stop_codons list.

    print('The find_stop function returned these stop indices', stop_indices)  # Here we simply print the list of beginning indices of all stop codons from DNA_seq.
    print('The list of codons from the random DNA sequence having beginning indices as indicated in the list above', stop_codons)  # All of these should be either 'TAA' or 'TAG' or 'TGA'
    print()
    print()

    # -------------------- Test 16: test orf_finder --------------------
    """
    In this test of the orf_finder function, I create a DNA sequence having no start or stop codons so that both find_start and find_stop will return empty lists for
    the start and stop indices. The purpose of this test is to illustrate that the orf_finder function can handle such situations without crashing due to some error.
    """
    print()
    print('---------- Test 16: Testing orf_finder ----------')
    print()
    DNA_seq = np.array(['A', 'A', 'C', 'G', 'C', 'C', 'A', 'G', 'A', 'G',  # Here I build the test DNA sequence in numpy array format.
                        'G', 'A', 'G', 'A', 'C', 'A', 'A', 'C', 'C', 'A',
                        'C', 'A', 'A', 'G', 'C', 'A', 'G', 'A', 'C', 'A'])
    [orf_for, orf_bak] = orf_finder(DNA_seq)  # Here we find all the possible open reading frames using the orf_finder function
    DNA_seq_str = seq2str(DNA_seq)  # Convert the DNA sequence to a string
    print('The constructed DNA sequence is', DNA_seq_str)  # Print the DNA sequence
    print('The DNA sequence was constructed to have the following collection of open reading frames [[], [], []]')  # Print the desired output of the function
    print('The list of open reading frames (for each of the three possible reading frames) of the DNA sequence is', orf_for)  # Print the list of open reading frames of DNA_seq
    print('The list of open reading frames (for each of the three possible reading frames) of the reverse complement of the DNA sequence is', orf_bak)  # List of orfs for DNA_seq's reverse complement
    print()
    print()

    # -------------------- Test 17: test orf_finder --------------------
    """
    In this test of the orf_finder function, I create a specific DNA sequence with a known collection of open reading frames. The DNA sequence is passed
    to orf_finder, which will find all possible open reading frames for the DNA sequence (orf_for) and its reverse complement (orf_bak).
    Note: by "open reading frame" I mean a tuple of indices that represent the actual open reading frame of the DNA sequence. The first element of a tuple
    is the beginning index of the open reading frame (i.e., the index of the first nucleotide of the start codon), and the second element of a tuple is the
    last index of the open reading frame (i.e., the index of the third nucleotide of the corresponding stop codon). After both lists of reading frames are
    returned I print each and compare them to the known reading frames. For the sequence constructed below, the open reading frame for the 0th reading frame
    is (18, 47). The open reading frame for the 1st reading frame is (1, 36). The open reading frame for the 2nd reading frame is (14, 43). So the function
    should return [[(18, 47)], [(1, 36)], [(14, 43)]]. Due to complexity, the sequence was not designed to produce any specific collection of open reading
    frames for the reverse complement of the DNA sequence. However, as long as the revcomp function works properly, the process for computing the open reading
    frames of the reverse complement is identical to the one for the constructed sequence. Hence, we can assume it works properly.
    """
    print()
    print('---------- Test 17: Testing orf_finder ----------')
    print()
    DNA_seq = np.array(['A', 'A', 'T', 'G', 'C', 'C', 'T', 'T', 'A', 'T',  # Here I build the test DNA sequence in numpy array format.
                        'C', 'G', 'C', 'A', 'A', 'T', 'G', 'G', 'A', 'T',
                        'G', 'C', 'T', 'T', 'A', 'C', 'A', 'A', 'G', 'T',
                        'T', 'A', 'T', 'T', 'T', 'A', 'A', 'T', 'C', 'A',
                        'T', 'T', 'A', 'G', 'C', 'T', 'G', 'A', 'C', 'A'])
    [orf_for, orf_bak] = orf_finder(DNA_seq)  # Here we find all the possible open reading frames using the orf_finder function
    DNA_seq_str = seq2str(DNA_seq)  # Convert the DNA sequence to a string
    print('The constructed DNA sequence is', DNA_seq_str)  # Print the DNA sequence
    print('The DNA sequence was constructed to have the following collection of open reading frames [[(18, 47)], [(1, 36)], [(14, 43)]]')  # Print the desired output of the function
    print('The list of open reading frames (for each of the three possible reading frames) of the DNA sequence is', orf_for)  # Print the list of open reading frames of DNA_seq
    print('The list of open reading frames (for each of the three possible reading frames) of the reverse complement of the DNA sequence is', orf_bak)  # List of orfs for DNA_seq's reverse complement
    print()
    print()

    # -------------------- Test 18: test orf_finder --------------------
    """
    In this test of the orf_finder function, I create a specific DNA sequence with a known collection of open reading frames. The only difference between this
    test and the previous test is that this one inserts extra start codons (in frame) between existing start/stop codons which are already in frame. The inserted
    start codons do not change the the true open reading frames of this sequence, so the output should be identical to the output from the previous test. The DNA
    sequence is passed to orf_finder, which will find all possible open reading frames for the DNA sequence (orf_for) and its reverse complement (orf_bak).
    Note: by "open reading frame" I mean a tuple of indices that represent the actual open reading frame of the DNA sequence. The first element of a tuple
    is the beginning index of the open reading frame (i.e., the index of the first nucleotide of the start codon), and the second element of a tuple is the
    last index of the open reading frame (i.e., the index of the third nucleotide of the corresponding stop codon). After both lists of reading frames are
    returned I print each and compare them to the known reading frames. For the sequence constructed below, the open reading frame for the 0th reading frame
    is (18, 47). The open reading frame for the 1st reading frame is (1, 36). The open reading frame for the 2nd reading frame is (14, 43). So the function
    should return [[(18, 47)], [(1, 36)], [(14, 43)]]. Due to complexity, the sequence was not designed to produce any specific collection of open reading
    frames for the reverse complement of the DNA sequence. However, as long as the revcomp function works properly, the process for computing the open reading
    frames of the reverse complement is identical to the one for the constructed sequence. Hence, we can assume it works properly.
    """
    print()
    print('---------- Test 18: Testing orf_finder ----------')
    print()
    DNA_seq = np.array(['A', 'A', 'T', 'G', 'C', 'C', 'T', 'A', 'T', 'G',  # Here I build the test DNA sequence in numpy array format.
                        'C', 'G', 'C', 'A', 'A', 'T', 'G', 'G', 'A', 'T',
                        'G', 'C', 'T', 'T', 'A', 'C', 'A', 'T', 'G', 'T',
                        'T', 'A', 'T', 'T', 'T', 'A', 'A', 'T', 'G', 'A',
                        'T', 'T', 'A', 'G', 'C', 'T', 'G', 'A', 'C', 'A'])
    [orf_for, orf_bak] = orf_finder(DNA_seq)  # Here we find all the possible open reading frames using the orf_finder function
    DNA_seq_str = seq2str(DNA_seq)  # Convert the DNA sequence to a string
    print('The constructed DNA sequence is', DNA_seq_str)  # Print the DNA sequence
    print('The DNA sequence was constructed to have the following collection of open reading frames [[(18, 47)], [(1, 36)], [(14, 43)]]')  # Print the desired output of the function
    print('The list of open reading frames (for each of the three possible reading frames) of the DNA sequence is', orf_for)  # Print the list of open reading frames of DNA_seq
    print('The list of open reading frames (for each of the three possible reading frames) of the reverse complement of the DNA sequence is', orf_bak)  # List of orfs for DNA_seq's reverse complement
    print()
    print()

    # -------------------- Test 19: test orf_finder --------------------
    """
    In this test of the orf_finder function, I create a specific DNA sequence with a known collection of open reading frames. The only difference between this
    test and Test 16 is that an extra start/stop codon pair is inserted into each of the 0th and 2nd reading frames. These two frames will each have two open
    reading frames and the 1st reading frame will have only the same one from the previous two tests. The DNA sequence is passed to orf_finder, which will find
    all possible open reading frames for the DNA sequence (orf_for) and its reverse complement (orf_bak). Note: by "open reading frame" I mean a tuple of indices
    that represent the actual open reading frame of the DNA sequence. The first element of a tuple is the beginning index of the open reading frame (i.e., the
    index of the first nucleotide of the start codon), and the second element of a tuple is the last index of the open reading frame (i.e., the index of the third
    nucleotide of the corresponding stop codon). After both lists of reading frames are returned I print each and compare them to the known reading frames. For
    the sequence constructed below, the open reading frame for the 0th reading frame is (18, 47). The open reading frame for the 1st reading frame is (1, 36). The
    open reading frame for the 2nd reading frame is (14, 43). So the function should return [[(18, 26), (36, 47)], [(1, 36)], [(5, 13), (14, 43)]]. Due to complexity, the
    sequence was not designed to produce any specific collection of open reading frames for the reverse complement of the DNA sequence. However, as long as the
    revcomp function works properly, the process for computing the open reading frames of the reverse complement is identical to the one for the constructed sequence.
    Hence, we can assume it works properly.
    """
    print()
    print('---------- Test 19: Testing orf_finder ----------')
    print()
    DNA_seq = np.array(['A', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'A', 'T',  # Here I build the test DNA sequence in numpy array format.
                        'C', 'T', 'A', 'G', 'A', 'T', 'G', 'G', 'A', 'T',
                        'G', 'C', 'T', 'T', 'T', 'A', 'A', 'A', 'G', 'T',
                        'T', 'A', 'T', 'T', 'T', 'A', 'A', 'T', 'G', 'A',
                        'T', 'T', 'A', 'G', 'C', 'T', 'G', 'A', 'C', 'A'])
    [orf_for, orf_bak] = orf_finder(DNA_seq)  # Here we find all the possible open reading frames using the orf_finder function
    DNA_seq_str = seq2str(DNA_seq)  # Convert the DNA sequence to a string
    print('The constructed DNA sequence is', DNA_seq_str)  # Print the DNA sequence
    print('The DNA sequence was constructed to have the following collection of open reading frames [[(18, 26), (36, 47)], [(1, 36)], [(5, 13), (14, 43)]]')  # Print the desired output of the function
    print('The list of open reading frames (for each of the three possible reading frames) of the DNA sequence is', orf_for)  # Print the list of open reading frames of DNA_seq
    print('The list of open reading frames (for each of the three possible reading frames) of the reverse complement of the DNA sequence is', orf_bak)  # List of orfs for DNA_seq's reverse complement
    print()
    print()

    # -------------------- Test 20: test orf_finder --------------------
    """
    In this test of the orf_finder function, I will generate a random DNA sequence of length 1000, with GC content set to the default value of 0.5.
    The AT & GC skews are set to their default values of zero. The DNA sequence is passed to orf_finder, which will find all possible open reading
    frames for the randomly generated DNA sequence (orf_for) and its reverse complement (orf_bak). Note: by "open reading frame" I mean a tuple of
    indices that represent the actual open reading frame of the sequence. The first element of a tuple is the beginning index of the open reading
    frame (i.e., the index of the first nucleotide of the start codon), and the second element of a tuple is the last index of the open reading frame
    (i.e., the index of the third nucleotide of the corresponding stop codon). After both lists of reading frames are returned I print each for inspection.
    """
    print()
    print('---------- Test 20: Testing orf_finder ----------')
    print()
    DNA_seq = mkseq(1000)  # Here we generate a random DNA sequence using our mkseq function.
    [orf_for, orf_bak] = orf_finder(DNA_seq)  # Here we find all the possible open reading frames using the orf_finder function
    print('The list of open reading frames (for each of the three possible reading frames) of the DNA sequence is:')
    print(orf_for)  # Print the list of open reading frames (for each of the 3 possible reading frames) of DNA_seq
    print()
    print('The list of open reading frames (for each of the three possible reading frames) of the reverse complement of the DNA sequence is:')
    print(orf_bak)  # Print the list of open reading frames (for each of the 3 possible reading frames) of DNA_seq's reverse complement
    print()
    print()

    # -------------------- Test 21: test filter_int --------------------
    """
    In this test of the filter_int function, I create a specific DNA sequence of length 30. That sequence is 'AGTGTACCAGTCGCGTATGCTAGTCTCAGG'.
    According to the definition of an intron defined in the function, this sequence has two introns that should be removed:
    1) beginning at index 1, ending at index 9 - GTGTACCAG
    2) beginning at index 14, ending at index 22 - GTATGCTAG
    After removing these introns, the function should return the numpy array version of the string 'ATCGCTCTCAGG'.
    I do not design the test DNA sequence to return and specific intron-filtered reverse complement of DNA_seq. The process is identical, but simply
    performed on the reverse complement of the DNA sequence. Hence, the success of that process depends only on the success of the revcomp function
    as long as we have sufficiently shown that filter_int works properly on the original DNA sequence. After the function returns its two numpy arrays
    (one for DNA_seq and one for its reverse complement) I simply print the results and compare to the truth.
    """
    print()
    print('---------- Test 21: Testing filter_int ----------')
    print()
    DNA_seq = np.array(['A', 'G', 'T', 'G', 'T', 'A', 'C', 'C', 'A', 'G',   # Construct the DNA sequence having the introns described above
                        'T', 'C', 'G', 'C', 'G', 'T', 'A', 'T', 'G', 'C',
                        'T', 'A', 'G', 'T', 'C', 'T', 'C', 'A', 'G', 'G'])
    [filtered_DNA_seq, filtered_DNA_seq_revcomp] = filter_int(DNA_seq)  # Pass DNA_seq to filter_int and retrieve the intron-filtered sequences
    true_filtered_str = 'ATCGCTCTCAGG'  # This is string version of the correctly filtered DNA sequence
    true_filtered_seq = seq2str(true_filtered_str)  # This is numpy array version of the correctly filtered DNA sequence
    print('The original DNA sequence is', seq2str(DNA_seq))
    print('The correct intron-filtered string is', true_filtered_str)
    print('The string version of the sequence returned by filter_int is', seq2str(filtered_DNA_seq))
    print('The correct intron-filtered sequence is', str2seq(true_filtered_seq))
    print('The sequence returned by filter_int is', filtered_DNA_seq)
    print()
    print('The reverse complement of the original DNA sequence is', seq2str(revcomp(DNA_seq)))
    print('The string version of the reverse complement sequence returned by filter_int is', seq2str(filtered_DNA_seq_revcomp))
    print('The reverse complement sequence returned by filter_int is', filtered_DNA_seq_revcomp)
    print()
    print()

    # -------------------- Test 22: test filter_int --------------------
    """
    In this test of the filter_int function, I create a specific DNA sequence of length 30 having no introns to be removed. That sequence is 'AGTGTCCCAGTCGCGTATGCTCGTCTCATG'.
    The function should return the original numpy array. I do not design the test DNA sequence to return and specific intron-filtered reverse complement of DNA_seq. The
    process is identical, but simply performed on the reverse complement of the DNA sequence. Hence, the success of that process depends only on the success of the revcomp
    function as long as we have sufficiently shown that filter_int works properly on the original DNA sequence. After the function returns its two numpy arrays
    (one for DNA_seq and one for its reverse complement) I simply print the results and compare to the truth.
    """
    print()
    print('---------- Test 22: Testing filter_int ----------')
    print()
    DNA_seq = np.array(['A', 'G', 'T', 'G', 'T', 'C', 'C', 'C', 'A', 'G',   # Construct the DNA sequence having no introns to be removed
                        'T', 'C', 'G', 'C', 'G', 'T', 'A', 'T', 'G', 'C',
                        'T', 'C', 'G', 'T', 'C', 'T', 'C', 'A', 'T', 'G'])
    [filtered_DNA_seq, filtered_DNA_seq_revcomp] = filter_int(DNA_seq)  # Pass DNA_seq to filter_int and retrieve the intron-filtered sequences
    print('The original DNA sequence is', seq2str(DNA_seq))
    print('The string version of the sequence returned by filter_int is', seq2str(filtered_DNA_seq))
    print('The sequence returned by filter_int is', filtered_DNA_seq)
    print()
    print('The reverse complement of the original DNA sequence is', seq2str(revcomp(DNA_seq)))
    print('The string version of the reverse complement sequence returned by filter_int is', seq2str(filtered_DNA_seq_revcomp))
    print('The reverse complement sequence returned by filter_int is', filtered_DNA_seq_revcomp)
    print()
    print()

    # -------------------- Test 23: test filter_int --------------------
    """
    In this test of the filter_int function, I generate a random DNA sequence of length 3000, with GC content set to the default value of 0.5.
    The AT & GC skews are set to their default values of zero. After the function returns its two numpy arrays, the original DNA sequence,
    the intron-filtered sequence, and the inron-filtered reverse complement of the original sequence are printed (string versions only).
    """
    print()
    print('---------- Test 23: Testing filter_int ----------')
    print()
    DNA_seq = mkseq(3000)  # Here we generate a random DNA sequence using our mkseq function.
    [filtered_DNA_seq, filtered_DNA_seq_revcomp] = filter_int(DNA_seq)  # Pass DNA_seq to filter_int and retrieve the intron-filtered sequences
    print('The original DNA sequence is', seq2str(DNA_seq))
    print('The string version of the sequence returned by filter_int is', seq2str(filtered_DNA_seq))
    print()
    print('The reverse complement of the original DNA sequence is', seq2str(revcomp(DNA_seq)))
    print('The string version of the reverse complement sequence returned by filter_int is', seq2str(filtered_DNA_seq_revcomp))
    print()
    print()