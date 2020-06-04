from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, _verify_alphabet
from Bio import Entrez, SeqIO
import numpy as np
import random, string
import pandas as pd
import joblib
import os

def iFeature_seq_extract(seq, desc, K=None):
    advanced_descs = ['CKSAAP', 'CKSAAGP', 'KSCTriad']
    rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(7))
    with open('temp/'+rand+'.fasta', 'w') as temp_seq_file:
        temp_seq_file.write('>'+rand+'\n'+seq)

    try:
        if desc in advanced_descs:
            if K != None:
                os.system('python3 iFeature/codes/{}.py temp/{} {} temp/{}'.format(desc, rand+'.fasta', K, rand+'.tsv'))
            elif K == None:
                os.remove('temp/'+rand+'.fasta')
                raise ValueError('What is the K value ?')
        else:
            os.system('python3 iFeature/iFeature.py --file temp/{} --type {} --out temp/{}'.format(rand+'.fasta', desc, rand+'.tsv'))
        
        outp = pd.read_csv('temp/'+rand+'.tsv', sep='\t', header=0)
        outp = outp.drop(['#'], axis=1)
        # outp = outp.astype('float').dtype
        os.remove('temp/'+rand+'.fasta')
        os.remove('temp/'+rand+'.tsv')
        return outp

    except:
        os.remove('temp/'+rand+'.fasta')
        print('Something went wrong, Please check for corrections')
        exit()

def PFeature_seq_extract(seq):
    rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(7))
    os.system('cd PFeat && python3 PFeats.py Extract_all {} {}'.format(seq, rand))
    res = pd.read_csv('PFeat/temp/'+ rand + '.csv')
    os.remove('PFeat/temp/'+ rand + '.csv')
    return res

def feature_extraction(seq):
    descriptors = ['AAC', 'CKSAAP', 'DDE', 'GAAC', 'CKSAAGP', 'GDPC', 'GTPC','Moran',
                'Geary', 'NMBroto', 'CTDC', 'CTDT', 'CTDD', 'KSCTriad', 'SOCNumber', 'QSOrder', 
                'PAAC', 'APAAC', 'PFeat']

    f_vector = []
    for des in descriptors:
        if des == 'CKSAAP':
            feats = iFeature_seq_extract(seq, des, K=5)
        elif des == 'CKSAAGP':
            feats = iFeature_seq_extract(seq, des, K=5)
        elif des == 'KSCTriad':
            feats = iFeature_seq_extract(seq, des, K=5)
        elif des == 'PFeat':
            feats = PFeature_seq_extract(seq)
        else:
            feats = iFeature_seq_extract(seq, des)

        f_vector.append(feats)
    f_vector = pd.concat(f_vector, axis=1) #as pd.df
    # f_vector = np.array(f_vector) #as np.array
    return f_vector


def read_fasta(file_name):
    
    seqsdic = {}
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    for fasta in fasta_sequences:
        seqsdic[fasta.id] = str(fasta.seq)
    return seqsdic


def obtain_seq_from_ent(ent):
    
    Entrez.email = "mehdiforoozandehsh@gmail.com"

    try:
        handle = Entrez.efetch(db="protein", id=ent, rettype="gp")
        record = SeqIO.read(handle, "gb")
        seq = str(record.seq)

    except:
        
        seq = str(pypro.GetProteinSequence(ent))

    if not _verify_alphabet(Seq(seq, IUPAC.protein)):
        print('Inserted Entry Is Not Valid!')
    return seq


def seq_repair(seq):
    forbidden_chars = ['B','J','O','U','X','Z','@','_','!','#','$','%','^','&','*','(',')','<','>','?','/','|','}','{','~',':','-', '\n', ' ']
    for f in forbidden_chars:
            while f in seq:
                seq = seq.replace(f,'')
    return seq


def ph_f_selection(dataset):
    with open('Models/pHdupidx.txt', 'r') as conf:
        idx = conf.readlines()[0].split(',')
    for i in range(len(idx)):
        idx[i] = int(idx[i])
    idx = np.array(idx)
    dataset = dataset[:, np.sort(idx)]
    scaler = joblib.load('Models/pH_Scaler.sav')
    dataset = scaler.transform(dataset)
    selector = joblib.load('Models/pH_Selector.sav')
    x = selector.transform(dataset)
    return x


def temp_f_selection(dataset):
    with open('Models/tempdupidx.txt', 'r') as conf:
        idx = conf.readlines()[0].split(',')
    for i in range(len(idx)):
        idx[i] = int(idx[i])
    idx = np.array(idx)
    dataset = dataset[:, np.sort(idx)]
    scaler = joblib.load('Models/temp_Scaler.sav')
    dataset = scaler.transform(dataset)
    selector = joblib.load('Models/temp_Selector.sav')
    x = selector.transform(dataset)
    return x


def phpredict(x):
    voting = joblib.load('Models/pH_model.sav')
    elected = voting.predict(x).reshape(1, -1)[0]
    # invkey = {0: 'Extreme Acidic', 1: 'Acidic', 2: 'Neutral', 3: 'Basic'}
    # elected = invkey[elected]
    return elected


def tempredict(x):
    voting = joblib.load('Models/temp_model.sav')
    elected = voting.predict(x).reshape(1, -1)[0]
    # invkey = {0: 'Mesophilic', 1: 'Thermophilic', 2: 'HyperThermophilic'}
    # elected = invkey[elected]
    return elected


# _____TOOL1: Prediction_____
def single_prediction(inp):
    if len(inp) > 20 and _verify_alphabet(Seq(inp, IUPAC.protein)):
        seq = inp
    else:
        seq = obtain_seq_from_ent(inp)
    seq = seq_repair(seq)
    full_features = np.array(feature_extraction(seq))
    ph_x = ph_f_selection(full_features)
    opt_ph = phpredict(ph_x)[0]
    temp_x = temp_f_selection(full_features)
    opt_temp = tempredict(temp_x)[0]
    return opt_temp, opt_ph


def fasta_prediction(file_name):
    seqsdic = read_fasta(file_name)
    result = []
    for seq_id, seq in seqsdic.items():
        seq = seq_repair(seq)
        full_features = np.array(feature_extraction(seq))
        ph_x = ph_f_selection(full_features)
        opt_ph = phpredict(ph_x)[0]
        temp_x = temp_f_selection(full_features)
        opt_temp = tempredict(temp_x)[0]
        result.append([str(seq_id), opt_ph, opt_temp])
        # result[str(seq_id)] = ' pH Dependence Prediction: ' + opt_ph + \
                              # ', Temperature Dependence Prediction: ' + opt_temp
    return result
