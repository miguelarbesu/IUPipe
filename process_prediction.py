# Created by @miguelarbesu on 11/05/2018

# Last update: 11/05/2018

# This script is designed to read IUPred output 
# from the predict_effectors.py script, load
# it as grouped Dataframe by name, and then
# calculate and save derived parameters - e.g.
# average disorder fraction by protein, etc.
# Output parameters are saved in an output directory
# as .tsv files to be later loaded and merged with
# the mapping.
# 
# For a complete run please provide:
# input_dir(str), output_dir(str)

import sys
import glob
import json
import pandas as pd

def load_dataset(file):
    '''Function:
    Loads IUPred an output file from the 
    predict_effectors.py script and returns 
    fasta name, prediction type, and
    grouped DataFrame by protein.

    input_dir(str) -> 
    fasta_name(str), pred_type(str),
    datagp (gp(pd.DF)
    '''
    dataset_name = file.split("/")[-1].split(".")[0]
    fasta_name = "_".join(dataset_name.split("_")[0:2])
    pred_type = dataset_name.split("_")[-1]
    datadf = pd.read_csv(filepath_or_buffer=file,
                         sep='\t', 
                         names=['point', 'resno', 'aa',
                                'score', 'prot'])
    datagp = datadf.groupby(by='prot')
    print("Loaded dataset {}, {} prediction type".format(fasta_name, pred_type))
    return(fasta_name, pred_type, datagp)
        #print(datadf)
        #disfrac = datagp.score.apply(lambda x: sum(x>0.4)/len(x))
        
def calc_disfrac(input_dir,
                 output_dir = "./processed_data/",
                 input_regexp = "*.dat",
                 disorder_threshold = 0.4):
    '''Function:
    Calculates disorder fraction for proteins
    in a grouped DataFrame according to a specified
    score threshold. The default value is 0.4 according
    to Fuxreiter et al. 2012?.
    Results are saved as a .tsv in an output folder
    
    input_dir(str), 
    input_regexp(str, opt),
    disorder_threshold(float) -> 
    datagp (gp(pd.DF)
    '''
    input_files = glob.glob(input_dir + input_regexp)
    for file in input_files:
        fasta_name, pred_type, datagp = load_dataset(file)
        disfrac = datagp.score.apply(lambda x: sum(x>disorder_threshold)/
                                     len(x))
        disfrac.to_csv(path = output_dir + '/{}_{}_disfrac.tsv'.format(
                       fasta_name, pred_type), 
                       sep='\t',
                       header = 'disorder_fraction')
        print("Calculated disorder fraction for dataset {}, {} prediction type".format(fasta_name, pred_type))
            
def main():
    calc_disfrac(*sys.argv[1:])
    
if __name__ == "__main__":
    main()
