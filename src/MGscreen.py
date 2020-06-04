from blast import blastx
import funcs
from datetime import datetime
import pandas as pd
import numpy as np
import os, shutil
import random, string
from functools import partial
import multiprocessing as mp

def gene_screen(filename, bs_filter):
    # _____RUN_BLAST_____
    blastresname = 'temp/result_blastx_' + ''.join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(7)) + '.txt'
    blastx(filename,
           'final_ref.fasta',
           blastresname)
    matches = []
    input_seqs = funcs.read_fasta(filename)
    sixpack = funcs.gen_frames(input_seqs)
    with open(blastresname, 'r') as file:
        blast_res = file.readlines()
        br = []
        for i in blast_res:
            br.append(i[:-1].split(','))
    del blast_res
    os.remove(blastresname)
    columns = ['qseqid', 'qlen', 'length', 'sseqid', 'qstart', 'qseq',
               'qend', 'sstart', 'send', 'pident', 'evalue', 'bitscore']
    df = pd.DataFrame(br, columns=columns)
    del br
    df = df.drop(['qstart', 'qend', 'sstart', 'send'], axis=1)
    # pruning non-homologous results
    # ref: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/

    to_drop = []
    for i in range(len(df)):
        if float(df.evalue[i]) > 1e-6:
            to_drop.append(i)
        elif float(df.bitscore[i]) < 50:
            to_drop.append(i)

    df = df.drop(to_drop)
    df = df.reset_index()

    # choose the best match
    tlis = {}
    for i in range(len(df)):
        if str(df.qseqid[i]) not in tlis.keys():
            tlis[str(df.qseqid[i])] = [i, float(df.bitscore[i])]
        else:
            if tlis[str(df.qseqid[i])][1] < float(df.bitscore[i]):
                tlis[str(df.qseqid[i])] = [i, float(df.bitscore[i])]

    p = list(np.arange(len(df)))
    q = [bs[0] for bs in tlis.values()]
    drop_list = [item for item in p if item not in q]
    df = df.drop(drop_list)
    df = df.reset_index()
    df = df.drop(['level_0', 'index'], axis=1)
    for i in range(len(df)):
        qid = df.qseqid[i]
        qseq = df.qseq[i].replace('-', '')
        for frame, trans in sixpack[qid].items():
            trans = trans.split('*')
            for contig in trans:
                if qseq[:20] in contig:
                    match_start = contig.find(qseq[:20])
                    bef = contig[:match_start]
                    if 'M' in bef:
                        bef = bef[::-1]
                        methionin = bef.find('M')
                        translated_region = contig[(match_start - methionin - 1):]
                    else:
                        translated_region = contig[contig.find(qseq[:20]):]
                    # seemingly Strict & powerful filters : Bit-score of 400_450 and a e-value(for this db) of 1e-150
                    if float(df.bitscore[i]) > bs_filter:
                        matches.append([qid, translated_region, frame, df.bitscore[i], df.evalue[i]])
    # matches_key: columns=['q_id', 'translation', 'frame', 'bitscore', 'evalue'])
    del df
    return matches


def bigfile_handler(filename, bs_filter):
    NUM_OF_LINES = 200000
    tempor_fol_name = 'temp_' + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(7))
    try:
        os.mkdir(tempor_fol_name)
        with open(filename) as fin:
            fout = open(tempor_fol_name + '/tempgene0.fasta', "w")
            for i, line in enumerate(fin):
                fout.write(line)
                if (i + 1) % NUM_OF_LINES == 0:
                    fout.close()
                    fout = open(tempor_fol_name + '/tempgene%d.fasta' % (i / NUM_OF_LINES + 1), "w")
            fout.close()
        temp_files = os.listdir(tempor_fol_name)
        for i in range(len(temp_files)):
            temp_files[i] = tempor_fol_name + '/' + str(temp_files[i])
        results = []
        # mp.freeze_support()
        pool = mp.Pool()
        N = pool.map(partial(gene_screen, bs_filter=bs_filter), temp_files)
        for i in N:
            results += i
        # for tempor_name in temp_files:
        #     batch = gene_screen(tempor_fol_name + '/' + str(tempor_name), bs_filter)
        #     print('One batch screened & loaded: ' + str(len(batch)))
        #     results += batch
        #     del batch
        while [] in results:
            results.remove([])
        shutil.rmtree(tempor_fol_name)
        print(np.array(results).shape)
        return np.array(results)
    except:
        print("Something Went Wrong!")
        shutil.rmtree(tempor_fol_name)
        exit()
        


if __name__ == '__main__':
    t0 = datetime.now()

    # matches = bigfile_handler('metagene_test.fa',
    #                           bs_filter=400)
    # print(matches)
    # print(matches.shape)
    # dataframe = pd.DataFrame(np.array(matches), columns=['q_id', 'translation', 'trans_len',
    #                                                     'q_len', 'frame', 'bitscore', 'evalue'])
    # print(dataframe.shape)

    matches = bigfile_handler('sequence.fasta',
                              bs_filter=400)
    print(len(matches))
    matches = np.array(matches)
    print(matches.shape)
    dataframe = pd.DataFrame(matches, columns=['q_id', 'translation', 'frame', 'bitscore', 'evalue'])
    print(dataframe.shape)

    # seqsdic = {}
    # for i in range(len(dataframe)):
    #     seqsdic[str(dataframe.q_id[i])] = str(dataframe.translation[i])
    # funcs.write_fasta(seqsdic, 'screened_cellulases.fasta')
    # print('Run Time = ', time.time() - t0)
