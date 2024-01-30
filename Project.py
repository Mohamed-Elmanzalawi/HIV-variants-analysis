# Code Done by Group 1 for CIT656: Programming to Bioinformatics.
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
import Function_Consesus
import Function_CG_Content
import os
from Bio import Entrez
from Bio import Phylo
import matplotlib.pyplot as plt
from Bio import AlignIO
import numpy as np

print('Welcome to our Script')
print('please enter the location of the following files.')

# first will get the folder the user keep all the data in .
# then get the location of Clustalw2 which we will use it for the alignment.
# lastly, the 2 accession lists the user wants to precess.
Download_seq = input('1- Desired folder to keep all the results'
                     ' example:(F:\Files\PycharmProjects\Programming) ')
Alignment_exe = input('2- Clustalw2.exe folder example:(F:\Files\Clustral W) ')
First_Accession_List = input('please enter your first accession list '
                             'example:(MZ767421.1,MZ767468.1) ').split(',')
Second_Accession_List = input('please enter your second accession list').split(',')

# Using the given information we create our paths.
clustalw_exe = "%s\clustalw2.exe" % Alignment_exe
First_Alignment = '%s\Clustal_Alignment.ALN' % Download_seq
second_Alignment = '%s\second_Alignment.ALN' % Download_seq
consensus_Sequences = '%s\consensus_Sequences.fasta' % Download_seq
Third_Alignment = '%s\Main_Alignment.ALN' % Download_seq
Final_Alignment = '%s\Final_Alignment.ALN' % Download_seq
CG_content = '%s\CG_content.txt' % Download_seq
Dissimilar_Positions_txt = '%s\Dissimilar_Positions.txt' % Download_seq

# Then we start downloading our sequences by Entrez.efetch
# command using the lists that the user provided.
First_Sequences = "%s\First_Sequences.fasta" % Download_seq
Entrez.email = 'medo611@hotmail.com'
net_handle = Entrez.efetch(
    db="nucleotide", id=First_Accession_List, rettype="fasta", retmode="text")
out_handle = open(First_Sequences, "w")
out_handle.write(net_handle.read())
out_handle.close()
net_handle.close()
Second_Sequences = "%s\Second_Sequences.fasta" % Download_seq
net_handle_2 = Entrez.efetch(
    db="nucleotide", id=Second_Accession_List, rettype="fasta", retmode="text")
out_handle_2 = open(Second_Sequences, "w")
out_handle_2.write(net_handle_2.read())
out_handle_2.close()
net_handle_2.close()

# Since we are going to compare the consensus sequence of the first list with the
# consensus of the second list will first do MSA for each list separately using clustalw.
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=First_Sequences, outfile=First_Alignment)
clustalw_cline()
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=Second_Sequences, outfile=second_Alignment)
clustalw_cline()

# Then we will capture each sequence in each alignment and put it in a separate list
Alignment_1_list = []
Alignment_2_list = []
for seq_record in SeqIO.parse(First_Alignment, format='clustal'):
    Alignment_1_list.append(seq_record.seq)
for seq_record in SeqIO.parse(second_Alignment, format='clustal'):
    Alignment_2_list.append(seq_record.seq)

# And get the consensus sequence of both alignments using our already made function.
# Show the user both consensus sequences and write it in a file.
consensus_Seq_file_in = open(consensus_Sequences, "w")
consensus_Seq_1 = Function_Consesus.consensus(Alignment_1_list)
consensus_Seq_2 = Function_Consesus.consensus(Alignment_2_list)
consensus_Seq_file_in.write(">The_First_list_Consensus_Sequence" + "\n" + consensus_Seq_1 +
                            '\n\n' + ">The_Second_list_Consensus_Sequence" + "\n" + consensus_Seq_2)
consensus_Seq_file_in.close()

# Now we print our consensus Sequences for the user.
print('This is the consensus sequences :')
for seq_record in SeqIO.parse(consensus_Sequences, format='fasta'):
    print(seq_record.id)
    print(seq_record.seq)

# Here we will make our Third alignment between the two consensus sequences.
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=consensus_Sequences,
                                     outfile=Third_Alignment)
clustalw_cline()

# first we will create two empty lists and open a file to store the dissimilar region data.
Dissimilar_Positions_open = open(Dissimilar_Positions_txt, 'w')
A = []
y = []
# List A will capture the sequences in the clustal alignment using for loop.
for SeqRecord in AlignIO.read(Third_Alignment, 'clustal'):
    A.append(SeqRecord.seq)
# We will create a matrix using the numpy package.
profile = np.array(A)
difference = True
Starting_Gaps = True
# Using for loop we will go over each column in the matrix and give it 4 conditions.
for x in range(len(A[0])):
    # Condition 1 : this condition was made to neglect the first gaps in the alignment
    # so if it saw '-' gaps it should neglect the entire column as long as Staring_Gaps==True
    # which will be false only when the script find the first similar region.
    if '-' in profile[:, x] and Starting_Gaps:
        continue
    # Condition 2 : since sets neglect repetition,if an entire list contain the same value
    # transformed to a set it will have only one value so if an entire column has the same value
    # which is our similar region is converted to a set it will have one value,So when the length==1
    # that means this is a similar region and we neglect it and Starting_Gaps = False so that we
    # include gaps in our dissimilar region and not negect it.
    if len(set(profile[:, x])) == 1:
        Starting_Gaps = False
        difference = False
    # Condition 3 :  index has a dissimilar region,so we will add it to the list Y
    #  and difference will be True.
    if len(set(profile[:, x])) != 1:
        y.append(x)
        difference = True
    # Condition 4: the script will only get the index and print it if the difference is false.
    if not difference:
        # if the y list contains only one element that means it is only one index the dissimilarity exists,
        # So the script will print the following text.
        if len(y) == 1:
            Dissimilar_Positions_open.write('There is a dissimilar region at index %s \n' % y[0])
            # then the script will take that index and compare at in every sequence with the consensus
            # sequence to identify the mutated nucleotide using for loop.
            for SeqRecord in AlignIO.read(Third_Alignment, 'clustal'):
                if SeqRecord[y[0]] != A[0][y[0]]:
                    Dissimilar_Positions_open.write(
                        'Which is %s in %s and in the Main consensus  is %s \n' %
                        (SeqRecord[y[0]], SeqRecord.id, A[0][y[0]]))
            # After the script has printed the dissimilar region we had to clear the list y
            # to prevent it from overlapping with the next dissimilar region.
            y.clear()
        # if the list y has no values which may occur if the first part of the sequence has a similarity,so
        # difference will be false and the script will try to print y and this not a dissimilar region
        # So this condition was made to continue to the next iteration and skip the rest of the code
        # if that happened.
        if len(y) == 0:
            continue
        # However if y had values ,not just one ,but multiple adjacent dissimilar nucleotides this part will
        # start in the script where it will get the first value of y which is the start of the dissimilar region
        # and the last value which is the end and compare it in each sequence with the consensus sequence.
        if len(y) != 1:
            Dissimilar_Positions_open.write(
                'There is a dissimilar region start at index %s and end at index %s \n' % (y[0], y[-1]))
            for SeqRecord in AlignIO.read(Third_Alignment, 'clustal'):
                # While using a range of indexes the last index is not included
                # so we added 1 to the last value to include it.
                if SeqRecord.seq[y[0]:(y[-1] + 1)] != A[0][y[0]:(y[-1] + 1)]:
                    Dissimilar_Positions_open.write('Which is %s in %s and in the Main consensus  is %s \n' %
                                                    (SeqRecord.seq[y[0]:(y[-1] + 1)], SeqRecord.id,
                                                     A[0][y[0]:(y[-1] + 1)]))
            y.clear()

Dissimilar_Positions_open.close()

# Here we will make a file to that has the first consensus sequence
# and the second accession list sequences in a combined file
# which will be used for the phylogenetic tree.
Final_Fasta = '%s\Final_Fasta.fasta' % Download_seq
Final_Fasta_in = open(Final_Fasta, 'w')
Second_Sequences_in = open(Second_Sequences, 'r')
Final_Fasta_in.write(">The_First_CS" + "\n" + consensus_Seq_1
                     + '\n\n' + Second_Sequences_in.read())
Second_Sequences_in.close()
Final_Fasta_in.close()

# Then we will get the CG content in each sequence using for loop and our already made function.
# and write it in a file for the user.
# Also we will make a scatter plot representing these results.
IDs = []
CG_percentage = []
CG_content_open = open(CG_content, 'w')
for Sequence in SeqIO.parse(Final_Fasta, format='fasta'):
    IDs.append(Sequence.id)
    CG_percentage.append(Function_CG_Content.cg_content(Sequence.seq))
    CG_content_open.write('The sequence %s has CG content= %f \n'
                          % (Sequence.id, Function_CG_Content.cg_content(Sequence.seq)))
CG_content_open.close()

# we will create a create a scatter plot for the CG content and save it for the user.
plt.title('The CG Content', fontsize=20)
plt.scatter(x=IDs, y=CG_percentage, marker='o', color='r')

# Here will use the for loop to annotate each point on the scatter plot.
for i, txt in enumerate(CG_percentage):
    plt.annotate(txt, (IDs[i], CG_percentage[i]))
plt.xlabel('IDs')
plt.ylabel('CG Percentage %')
plt.xticks(rotation='vertical')
plt.grid()
plt.savefig("%s\The CG Content scatter plot" % Download_seq)

# Inorder to get the phylogenetic tree we have to make an alignment first with the sequences
# in the combined file which will produce the dnd file for the tree.
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=Final_Fasta, outfile=Final_Alignment)
clustalw_cline()

# Now we will use the dnd file to get the phylogenetic tree.
dnd_file = ("%s\Final_Fasta.dnd" % Download_seq)
tree = Phylo.read(dnd_file, "newick")
fig_2 = plt.figure(figsize=(13, 5), dpi=100)  # create figure & set the size
fig_2.suptitle('The Phylogenetic Tree',  # setting the title
               fontsize=20)
axes = fig_2.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes)
fig_2.savefig("%s\The Phylogenetic Tree" % Download_seq)  # saving it for the user

# Here we removed the unnecessary files file since it's no longer needed.
# to make the script produce only the needed files but this part optional.
os.remove(second_Alignment)
os.remove("%s\First_Sequences.dnd" % Download_seq)
os.remove("%s\Second_Sequences.dnd" % Download_seq)
os.remove('%s\consensus_Sequences.dnd' % Download_seq)
os.remove(Final_Alignment)
os.remove(Final_Fasta)
