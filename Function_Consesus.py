def consensus(Alignment_list):
    import numpy as np
    import textwrap
    # First will create a list for each nucleotide following by another two new lists.
    C = []
    G = []
    T = []
    A = []
    All = []
    consensus = []

    # We will take the alignment list from the user and convert it to array the rows will be
    # the number of sequences which length of the list and the columns will be the number of
    # nucleotides in any sequence in the list since they are all the same.
    M = np.array(Alignment_list).reshape(len(Alignment_list), len(Alignment_list[0]))

    # Then we will for loop to loop over evey column using range the length of any sequence.
    for base in range(len(Alignment_list[0])):
        # we will make a counter for every nucleotide.
        A_count = 0
        C_count = 0
        G_count = 0
        T_count = 0
        # Here we will loop in every column changing only the y axis=base every iteration.
        for NUC in M[:, base]:
            # Now if we find any of the nucleotides in the entire column we will add it to the counter.
            if NUC == "A":
                A_count += 1
            elif NUC == "C":
                C_count += 1
            elif NUC == "G":
                G_count += 1
            elif NUC == "T":
                T_count += 1
        # Now we append the values of each counter to its respective list.
        A.append(A_count)
        C.append(C_count)
        G.append(G_count)
        T.append(T_count)

    # After going over each column in the matrix we will append the values in each counter list to
    # a our new list 'p'.
    All.append(A)
    All.append(C)
    All.append(G)
    All.append(T)

    # Then we will change our 'p' list to a matrix.
    # Each column in this matrix  which has 4 rows will represent the
    # count of A,C,G,T respectively at that index.
    All_array = np.array(All).reshape(4, len(A))

    # Meaning if we can identify at which cell the max value is we will get the dominant nucleotide
    # since each value represent the count of that nucleotide at that index.
    for NUC in range(len(A)):
        # So we will go column by column using for loop.
        # if the max value was in the first row so nucleotide 'A' is the dominant.
        # and if it was in the second row it will be 'C' and if it is the third row then it is G
        # finally if it was the fourth row then it must be T.
        # this loop keeps going until we cover all the columns
        # and after we identify the dominant nucleotide we append it our new list consensus.
        if max(All_array[:, NUC]) == All_array[0, NUC]:
            consensus.append("A")
        elif max(All_array[:, NUC]) == All_array[1, NUC]:
            consensus.append("C")
        elif max(All_array[:, NUC]) == All_array[2, NUC]:
            consensus.append("G")
        elif max(All_array[:, NUC]) == All_array[3, NUC]:
            consensus.append("T")
    # Then we combine our list of nucleotides to form our consensus sequence.
    combined_consensus = "".join(consensus)
    # this one is optional : to make our consensus sequence more appealing to the fasta format
    # we used the function textwrap.fill which insert the desired number of elements in each line
    # in our case the line in the fasta file has 70 nucleotide.
    wrapped = textwrap.fill(combined_consensus, 70)
    return wrapped
