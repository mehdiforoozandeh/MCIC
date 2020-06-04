import os

def blastx(queryfile, database, outfile): 
    blast_options = './blastx -outfmt "10 qseqid qlen length sseqid qstart qseq qend sstart send pident evalue bitscore"'
    command = str(blast_options + ' -query ' + queryfile + ' -db ' + database + ' -out ' + outfile)
    os.system(command)


if __name__ == '__main__':
    blastx('sequence.fasta',
           'final_ref.fasta',
           'result_blastx.txt')
