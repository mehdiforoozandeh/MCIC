import multiprocessing
multiprocessing.freeze_support()
import MGscreen
import PredNGn
import funcs
import warnings
import sys
import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")


def Cel_Screen_Pred(filename, bs_filter=50, export_file=True, output_name=None):
    file_len = sum(1 for __ in open(filename))
    print('Screening...')
    if file_len > 200000:
        matches = MGscreen.bigfile_handler(filename, bs_filter)
        if len(matches)==0:
           print('No Cellulase Detected!')
           sys.exit()
    else:
        matches = np.array(MGscreen.gene_screen(filename, bs_filter))
        if len(matches)==0:
            print('No Cellulase Detected!')
            sys.exit()
    screened = pd.DataFrame(matches, columns=['q_id', 'translation', 'frame', 'bitscore', 'evalue'])
    print('Screening Results were loaded successfully! ')
    print('Predicting...')
    screened['pHpred'] = pd.Series()
    screened['Temppred'] = pd.Series()

    for i in range(len(screened)):
        seq = str(screened.translation[i])
        seq_id = str(screened.q_id[i])
        seq = PredNGn.seq_repair(seq)
        full_features = np.array(PredNGn.feature_extraction(seq))
        ph_x = PredNGn.ph_f_selection(full_features)
        opt_ph = PredNGn.phpredict(ph_x)
        temp_x = PredNGn.temp_f_selection(full_features)
        opt_temp = PredNGn.tempredict(temp_x)
        screened.pHpred[i] = str(opt_ph)
        screened.Temppred[i] = str(opt_temp)

    print('Prediction Done!')
    if export_file:
        print('Exporting result file ... ')
        if output_name==None:
            opn = filename
            while '/' in opn:
                opn = opn.replace('/','')
            output_name = 'results/MCIC_Results_' + str(opn)
        screened.to_csv(output_name + '.csv')
        print("Your prediction results have been written on file: " + output_name + '.csv')
    else:
        print(screened)


def Cel_Screen(filename, bs_filter=50, export_file=True, output_name=None):
    file_len = sum(1 for __ in open(filename))
    if file_len > 200000:
        matches = MGscreen.bigfile_handler(filename, bs_filter)
        if len(matches)==0:
           print('No Cellulase Detected!')
    else:
        matches = np.array(MGscreen.gene_screen(filename, bs_filter))
        if matches==[]:
            print('No Cellulase Detected!')
    screened = pd.DataFrame(matches, columns=['q_id', 'translation', 'frame', 'bitscore', 'evalue'])
    print('Screening Results were loaded successfully!!! ')
    if export_file:
        if output_name==None:
            opn = filename
            while '/' in opn:
                opn = opn.replace('/','')
            output_name = 'results/Screened_Cellulases_' + str(opn)
        screened.to_csv(output_name + '.csv')
        print('Your screening results have been written on '+output_name+'.csv')
    else:
        print(screened) 


def fasta_pred(filename, export_file=True, output_name=None):
    result = pd.DataFrame(PredNGn.fasta_prediction(filename), columns=['seq id', 'pH dependence', 'temperature dependence'])
    if export_file:
        if output_name==None:
            opn = filename
            while '/' in opn:
                opn = opn.replace('/','')
            output_name = 'results/Prediction_Results_' + str(opn)
        result.to_csv(output_name + '.csv')
        print('Your prediction results have been written on '+output_name+'.csv')
    else:
        print(result)


def single_pred(SeqEnt):
    t_opt, p_opt = PredNGn.single_prediction(SeqEnt)
    return 'pH Dependence Prediction: ' + p_opt + \
           ', Temperature Dependence Prediction: ' + t_opt


if __name__ == '__main__':

    # Pyinstaller fix
    helpme ='''Using a sequence similarity based annotation and an ensemble machine learning approach, MCIC (metagenome cellulase identification and classification) aims to identify and classify cellulolytic enzymes from a given metagenomic data as well as any other amino-acid sequence on the basis of optimum temperature and pH. 


*** USAGE:
(arguments are separated with a space)

-firts argument : MCIC 
-second argument options: [-h] [--help] [csp] [CelScreenPred] [cs] [CelScreen] [fp] [FastaPred] [sp] [SinglePred]
-third argument: [query input file]
-forth+ argument options: [-bs] [--bitscore] [-out] [-noexport]


*** FUNCTIONS: 

1. [csp] or [CelScreenPred]: Accepts nucleotide sequences and screens the cellulolytic enzymes and predicts pH and temperature dependence.
2. [cs] or [CelScreen]: Accepts nucleotide sequences and screens the cellulolytic enzymes. (without prediction)
3. [fp] or [FastaPred]: Accepts a fasta file containing protein sequence of cellulolytic enzymes and predicts pH and temperature dependence.    
4. [sp] or [SinglePred]: Accepts a single protein sequence or entry or accession number of cellulolytic enzymes and predicts pH and temperature dependence.


***HOW TO USE EACH FUNCTION:

1. [csp] or [CelScreenPred]: 

    options: 
    - [-bs]/[--bitscore]: choose your bitscore limit for filteration during the screening process. ==> "default: 50"
    - [-out]: choose your output file's name/address (*.csv). ==> "default: inputname.csv"  
    - [-noexport]: if chosen, the results will get printed on screen and no output file gets exported. "default: False"

    << By default, all results will be written in "results" folder.>>

    usage examples: 
    - MCIC csp input.format
    - MCIC CelScreenPred input.format -bs 100 -out John 
    - MCIC csp input.format --bitscore 500 -noexport
    - MCIC CelScreenPred input.format -out home/John
    - MCIC csp input.format -noexport

2. [cs] or [CelScreen]: 
    
    options: 
    - [-bs]/[--bitscore]: choose your bitscore limit for filteration during the screening process. ==> "default: 50"
    - [-out]: choose your output file's name/address (*.csv). ==> "default: inputname.csv"  
    - [-noexport]: if chosen, the results will get printed on screen and no output file gets exported. "default: False"

    << By default, all results will be written in "results" folder.>>

    usage examples: 
    - MCIC cs input.format
    - MCIC CelScreen input.format -bs 300 -out John 
    - MCIC cs input.format --bitscore 50 -noexport
    - MCIC CelScreen input.format -out home/John
    - MCIC cs input.format -noexport

3. [fp] or [FastaPred]:
    
    options:
    - [-out]: choose your output file's name/address (*.csv). ==> "default: inputname.csv"  
    - [-noexport]: if chosen, the results will get printed on screen and no output file gets exported. "default: False"

    << By default, all results will be written in "results" folder.>>

    usage examples: 
    - MCIC fp input.fasta
    - MCIC FastaPred input.fasta -out John
    - MCIC fp input.fasta -noexport

4. [sp] or [SinglePred]:
    
    This function does not have any options. 
    Two types of inputs are possible (1) amin acid sequence (2) protein entry/accession

    usage examples:
    - MCIC sp MKSCAILAALGCLA....
    - MCIC SinglePred Q7Z9M7

Thanks for reading...'''

    #dictionary = ['-h', '--help', '-bs', '--bitscore', '-out', '-noexport'
    #,'csp','CelScreenPred','cs','CelScreen','fp','FastaPred','sp', 'SinglePred']
    if len(sys.argv)==1: print('MCIC: missing operand \nTry -h or --help for more information.')
    elif len(sys.argv)==2:
        if sys.argv[1]!='-h' and sys.argv[1]!='--help':
            print('MCIC: missing operand \nTry -h or --help for more information.')
        elif sys.argv[1]=='-h' or '--help': print(helpme)
    elif len(sys.argv)>2:    
        # try:    
        if sys.argv[1]== 'csp' or sys.argv[1]=='CelScreenPred':
            bs=50
            export=True
            out=None
            if '-bs' in sys.argv:
                bs= int(sys.argv[sys.argv.index('-bs')+1])
            if '--bitscore' in sys.argv:
                bs= int(sys.argv[sys.argv.index('--bitscore')+1])
            if '-out' in sys.argv:
                out= str(sys.argv[sys.argv.index('-out')+1])
            if '-noexport' in sys.argv:
                export= False
            Cel_Screen_Pred(sys.argv[2], bs_filter=bs, output_name=out, export_file=export)

        elif sys.argv[1]== 'cs' or sys.argv[1]== 'CelScreen':
            bs=50
            export=True
            out=None
            if '-bs' in sys.argv:
                bs= int(sys.argv[(sys.argv.index('-bs'))+1])
            if '--bitscore' in sys.argv:
                bs= int(sys.argv[sys.argv.index('--bitscore')+1])
            if '-out' in sys.argv:
                out= str(sys.argv[sys.argv.index('-out')+1])
            if '-noexport' in sys.argv:
                export= False
            Cel_Screen(sys.argv[2], bs_filter=bs, output_name=out, export_file=export)

        elif sys.argv[1]== 'fp' or sys.argv[1]== 'FastaPred' :
            export=True
            out=None
            if '-out' in sys.argv:
                out= str(sys.argv[sys.argv.index('-out')+1])
            if '-noexport' in sys.argv:
                export= False
            fasta_pred(sys.argv[2], output_name=out, export_file=export)

        elif sys.argv[1]== 'sp' or sys.argv[1]== 'SinglePred':
            print(single_pred(sys.argv[2]))
        else:
            print('problem with the command!')
#         except:
#             print('''Error: 
# Something went wrong!
# Please check your input files.
# You can use -h or --help for more information about the functions, usage, and options.''')
