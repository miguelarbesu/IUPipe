# Created by @miguelarbesu on 17/04/2018

# Last update: 10/05/2018

# This is a simplified, compact version of the 
# original disorder prediction pipeline for HGTs for 
# the nonT3-T3-T4-T6 effectors based in the
# original module developed for the PAIdb.

# Functions can be imported separately or 
# either the whole pipeline run completely
# by excecuting this script.

# For a complete run please provide:
# input_dir(str), iupred_flag(str, "short" or "long")

import os
import sys
import glob
import json
import tempfile
import subprocess
import textwrap
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def map_prot_fasta(input_dir, input_regexp = "*_nr100.fasta"):
    '''Function:
    Reads multirecord .fasta files from an input folder,
    creates direct and inverse dictionaries between protein entries
    and the corresponding .fasta files, and saves them as JSON files
    in an output folder created ad hoc.
    input_dir(string), input_regexp(string, opt) -> 
    fasta2prot(JSON), prot2fasta(JSON)'''
    input_files = glob.glob(input_dir + input_regexp)
    output_dir = "iupred_"+input_dir
    fasta2prot = {}
    print(80*"-"+"\nMapping FASTAS and proteins\n"+80*"-")
    for filename in input_files:
        fastaid = filename.split("/")[-1].split("_")[0]
        print("Processing {}.fasta".format(fastaid))
        try:
            fasta2prot[fastaid] = list(record.id.split('|')[1] for record in SeqIO.parse(filename, format="fasta"))
        except IndexError:
            fasta2prot[fastaid] = list(record.id for record in SeqIO.parse(filename, format="fasta"))
        prot2fasta = {value: key for key in fasta2prot for value in fasta2prot[key]}
        
        with open(output_dir + "fasta2prot_map.json",
                  mode="w+",
                  newline="\n") as handle:
            json.dump(fasta2prot, handle, indent=4, sort_keys=True)
        
        with open(output_dir + "prot2fasta_map.json",
                  mode="w+",
                  newline="\n") as handle:
            json.dump(prot2fasta, handle, indent=4, sort_keys=True)
            
def map_prot_taxon(input_dir, input_regexp = "*_nr100.fasta"):
    '''Function:
    Maps each protein from all FASTA files in an input folder to its corresponding taxon.
    The resulting dictionary is saved as a JSON file in the corresponding output folder.
    input_dir(string) -> 
    prot2taxon(JSON)'''
    input_files = glob.glob(input_dir + input_regexp)
    output_dir = "iupred_"+input_dir
    prot2taxon = {}
    print(80*"-"+"\nMapping taxons and proteins\n"+80*"-")
    for filename in input_files:
        print("Processing {}".format(filename))
        handle = SeqIO.parse(filename, format="fasta")
        for record in handle:
            desc_list = record.description.split(" ")
            taxon = "".join([x.lstrip('OX=') for x in desc_list if x.startswith('OX=')])
            try:
                protid = record.id.split('|')[1]
            except IndexError:
                protid = record.id
            prot2taxon[protid] = taxon
        
            with open(output_dir + "prot2taxon_map.json",
                      mode="w+",
                      newline="\n") as handle:
                json.dump(prot2taxon, handle, indent=4, sort_keys=True)

def run_iupred(input_dir,
               iupred_flag,
               input_regexp = "*_nr100.fasta",
               iupred_path = "./software/iupred/"):
    '''Function: Runs IUPred with the specified flag
    over multi-record FASTA files. An input folder
    and a regular expression are also required.
    '''
    my_env = os.environ.copy()
    my_env["IUPred_PATH"] = iupred_path
    
    input_files = glob.glob(input_dir + input_regexp)
    output_dir = "iupred_"+input_dir
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass
    
    for file in input_files:
        fileref = file.split("/")[-1].split(".")[0]
        resfile = output_dir + fileref + "_IUpred_res_" + iupred_flag + ".dat"
        tempres = tempfile.NamedTemporaryFile(mode='a+')
        indexlist = []
        fasta_handle = SeqIO.parse(file, format='fasta')
        for record in fasta_handle:
            try:
                rec_id = record.id.split('|')[1]
            except IndexError:
                rec_id = record.id
            seq = textwrap.fill(str(record.seq), width=60)
            indexlist.extend([rec_id] * len(record.seq))
            tempseq = tempfile.NamedTemporaryFile(mode='a+')
            print('Predicting record {} in tempfile {}'.format(
                rec_id, tempseq.name))
            tempseq.write('>{}\n{}'.format(rec_id, seq))
            # This is very important, makes it go to the beginning of the file!
            tempseq.seek(0)
            # If stdout is set to none, it will stream prediction to the terminal.
            iupredproc = subprocess.Popen(
                [iupred_path + 'iupred', tempseq.name, iupred_flag], env=my_env, stdout=tempres)
            iupredproc.communicate()
            tempseq.close()
        df = pd.read_csv(filepath_or_buffer=tempres.name,
                         header=None, sep=' ', skipinitialspace=True,
                         comment='#', names=['Position', 'Residue', 'Score'])
        indexseries = pd.Series(data=[x for x in indexlist], name='Protein_ID')
        df['Protein_ID'] = indexseries
        df.to_csv(path_or_buf=resfile, sep="\t", header=False)        

def main():
    map_prot_fasta(sys.argv[1])
    #map_prot_taxon(sys.argv[1])
    #run_iupred(*sys.argv[1:])
    
if __name__ == "__main__":
    main()
      


