# Function to read fasta files
def read_fasta(file_name):
    from Bio import SeqIO
    seqsdic = {}
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    for fasta in fasta_sequences:
        seqsdic[fasta.id] = str(fasta.seq)
    return seqsdic


# Writes a dictionary to a fasta file

def write_fasta(dictionary, filename):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when calling the function
    """
    with open(filename, "w") as outfile:
        for key, value in dictionary.items():
            outfile.write(">" + str(key) + "\n" + str(value) + "\n")

    # print("Success! File written")


# Swaps DNA sequencs for proteins

def swap_dna(dnastring):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    protein = []
    end = len(dnastring) - (len(dnastring) % 3) - 1
    for i in range(0, end, 3):
        codon = dnastring[i:i + 3]
        if codon in table:
            aminoacid = table[codon]
            protein.append(aminoacid)
        else:
            protein.append("X")
    return "".join(protein)


# Required function for the frame_id function


def rev_seq(seq):
    trans = []
    for i in seq:
        if i == 'A':
            trans.append('T')
        elif i == 'C':
            trans.append('G')
        elif i == 'G':
            trans.append('C')
        elif i == 'T':
            trans.append('A')
        else:
            trans.append(i)
    trans = ''.join(trans)
    seq_rev = trans[::-1]
    return seq_rev


# Generates the six possible frames per one sequence

def frame_id(seq):
    """
    Usage: frame_id(dictionary['key'])
    frame_id(sequences['>Seq1'])
    six_frames = frame_id(sequences['>Seq1'])
"""
    frames = {'f0': [], 'f+1': [], 'f+2': [], 'r0': [], 'r-1': [], 'r-2': []}
    seq_rev = rev_seq(seq)
    for j in range(0, 3):
        temp = ''.join([seq[j:]])
        temp_rev = ''.join([seq_rev[j:]])
        seq_trans = swap_dna(temp)
        seq_rev_trans = swap_dna(temp_rev)
        if j == 0:
            frames['f0'] = seq_trans
            frames['r0'] = seq_rev_trans
        if j == 1:
            frames['f+1'] = seq_trans
            frames['r-1'] = seq_rev_trans
        if j == 2:
            frames['f+2'] = seq_trans
            frames['r-2'] = seq_rev_trans

    return frames


# Generates all the frames for all the sequences

def gen_frames(dictionary):
    all_dict = {}
    for key, value in dictionary.items():
        all_dict[key] = frame_id(dictionary[key])

    return all_dict


def get_ORF(sixpack_dic):
    for k in sixpack_dic.values():
        for j, q in k.items():
            k[j] = list(q.split('*'))
    return sixpack_dic


#  For command line Usage

import getopt
import sys


def main(argv):
    inputfile = ""
    outputfile = ""
    printprots = False
    try:
        opts, args = getopt.getopt(argv, "hi:o:p", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('dna2proteins.py -i <inputfile> -o <outputfile> -p')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('dna2proteins.py -i <inputfile> -o <outputfile> -p')
            print("-p prints the protein sequences in the terminal")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt == "-p":
            printprots = True

    return inputfile, outputfile, printprots


if __name__ == "__main__":
    inputfile, outputfile, printprots = main(sys.argv[1:])
