from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
import re
## Using SeqIO to read fasta data
file = input("Please input the fasta file: ")
record = SeqIO.read(file, "fasta")
## Using translation table 11 which is for bacterial
table = 11
## Set the min protein length, so won't include amino acid sequence that len < 50 
min_pro_len = 50
## Create a list to store protein sequences
orf = []

for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
    for frame in range(3):
        length = 3 * ((len(record)-frame) // 3) #Multiple of three
        nc = nuc[frame:frame+length].translate(table)
        com = re.compile(r"M([A-Z]*)\*")
        nuc_list = re.findall(com, str(nc))
        for pro in nuc_list:
            if len(pro) >= min_pro_len:
                orf.append(pro)
print("Using this script, it can find ", len(orf), "of open reading frames in target genomics")
### Calculate Moleculer Weight for amino acid sequences
## Create a amino acid weight map
amino_acid_weight = {
    "G": 75.07, "A": 89.09, "V": 117.1, "L": 131.2, "F": 165.2,
    "Y": 181.2, "N": 132.1, "Q": 146.1, "S": 105.09, "T": 119.1,
    "C": 121.2, "I": 131.2, "M": 149.2, "P": 115.1, "W": 204.2,
    "D": 133.1, "E": 147.1, "H": 155.2, "K": 146.2, "R": 174.2
}
### Create empty map for amino acid sequence and its weight
### Create empty list for amino acid sequence which len < 1000
sequence_len = []
sequence_weight = {}
### Using for loop to count every sequences
with open("Weight.txt", "w") as f:
    for i in orf:
    ### initial count to 0
        count = 0
    ### Count every amino acids in sequence
        for j in i:
        ### If the amino acid matchs to the amino acid in map, then takes value and adds to count
            count += amino_acid_weight[j]
    ### After counting, add count to map that the amino acid sequence is key, and weight is value
        sequence_weight[i] = round(count, 2)
    ### Print the sequences and their weights
        f.write(f"Amino acid sequence: {i}\nWeight: {sequence_weight[i]}\n")
        if 600 < len(i) < 900:
            sequence_len.append(i)
#### Blast ########
for i in range(5):
    print("Start Blast")
    i_seq = Seq(sequence_len[i])
    result_handle = NCBIWWW.qblast("blastp", "pdb", i_seq)
    print("End blast")
    blast_records = NCBIXML.read(result_handle)
    print("Write File")
    with open(f"my_blast_{i}.txt", "w") as f:
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                f.write('****Alignment****\n')
                f.write('sequence:' + alignment.title + "\n")
                f.write('length:' + str(alignment.length) + "\n")
                f.write('e value:' + str(hsp.expect) + "\n")
                f.write(hsp.query[0:75] + '...\n')
                f.write(hsp.match[0:75] + '...\n')
                f.write(hsp.sbjct[0:75] + '...\n')
