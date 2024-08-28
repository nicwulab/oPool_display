from collections import Counter
from rich.progress import track
from Bio import AlignIO
import pandas as pd
import numpy as np
import itertools
import random
import time, os
import multiprocessing

# codon from: https://www.biologicscorp.com/tools/CodonUsage

C_Table = {"UUU":"F","UCU":"S","UAU":"Y","UGU":"C",
"UUC":"F","UCC":"S","UAC":"Y",
"UUA":"L","UCA":"S","UAA":"*","UGA":"*",
"UCG":"S","UGG":"W",
"CUU":"L","CCU":"P","CAU":"H","CGU":"R",
"CCC":"P","CGC":"R",
"CCA":"P","CAA":"Q",
"CCG":"P","CAG":"Q",
"AUU":"I","ACU":"T","AAU":"N","AGU":"S",
"AUC":"I","ACC":"T","AAC":"N","AGC":"S",
"ACA":"T","AAA":"K",
"AUG":"M","ACG":"T","AAG":"K",
"GUU":"V","GCU":"A","GAU":"D","GGU":"G",
"GUC":"V","GCC":"A","GAC":"D","GGC":"G",
"GUA":"V","GCA":"A","GAA":"E","GGA":"G",
"GUG":"V","GCG":"A","GAG":"E","GGG":"G"}

def read_F2D(Fasta_DB):
    with open(Fasta_DB, 'r') as F:
        Seqs = F.read().split('\n')
    Seq_dict = {}
    for i in track(range(int(len(Seqs)/2)), description="Turning Fasta into dict ..."):
        Seq_dict.update({Seqs[i*2].replace(">", ""):Seqs[i*2+1]})
    return Seq_dict

TB_VL = pd.read_csv('result/VL.tsv.gz', sep ='\t')

#TB_VL[TB_VL.sequence_id == id].column

def Dict_update(Seq_dict, TB):
    for id in Seq_dict.keys():
        Seq_dict[id] = {'Nr': Seq_dict[id]}
        Seq_dict[id].update({
            "CDRL1" :TB.CDRL1_AA[TB.Name == id].iloc[0],
            "CDRL2" : TB.CDRL2_AA[TB.Name == id].iloc[0],
            "CDRL3" : TB.CDRL3_AA[TB.Name == id].iloc[0],
            "CDRH1" : TB.CDRH1_AA[TB.Name == id].iloc[0],
            "CDRH2" : TB.CDRH2_AA[TB.Name == id].iloc[0],
            "CDRH3" : TB.CDRH3_AA[TB.Name == id].iloc[0],
            "VL"    : TB.VL_AA[TB.Name == id].iloc[0],
            "VH"    : TB.VH_AA[TB.Name == id].iloc[0],
            "AA"   : TB.VL_AA[TB.Name == id].iloc[0] + "GGGGS" * 3 + TB.VH_AA[TB.Name == id].iloc[0]
        })
    return Seq_dict

def Check_quality(Seq_dict):
    print([[id, "".join([C_Table[Seq_dict[id]['Nr'][i*3:i*3+3]] for i in range(int(len(Seq_dict[id]['Nr'])/3))]), Seq_dict[id]['AA']]  for id in Seq_dict if "".join([C_Table[Seq_dict[id]['Nr'][i*3:i*3+3]] for i in range(int(len(Seq_dict[id]['Nr'])/3))])!=Seq_dict[id]['AA']  ])

def Blast_check():
    TMP_tmp = pd.DataFrame()
    i = 3
    while len(TMP_tmp) != 25 *3:
        i += 1
        os.system(f"blastn -query Primer/{group}  -db blastDB/{group}  -outfmt '6 qacc sacc evalue pident qcovs' -evalue 1e-1 -num_threads 8 -max_hsps 2 -word_size {i} > result/blast/{group}")
        os.system("cat result/blast/" + group + "| awk '{print $1}'| sort| uniq -c| awk '{print $1,$2}'| grep '^1 ' >  result/segs_id/" + group + '.list')
        with open(f"result/segs_id/{group}.list", 'r') as F:
            TMP = F.read().split('\n')[:-1]
        LIST = [i.split(' ')[-1].split(':') for i in TMP]
        TMP_tmp = pd.concat([TMP_tmp, pd.DataFrame(LIST)])
        TMP_tmp[2] = TMP_tmp[2].astype(int)
        TMP_tmp = TMP_tmp.sort_values(list(TMP_tmp.columns))
        TMP_tmp = TMP_tmp[~TMP_tmp.iloc[:,:2].duplicated()]
    TMP_tmp.to_csv(f'result/segs_id/{group}.csv')


# read the table

Group = [i for i in os.listdir('result') if "Re_assembled_" in i]
TB = pd.read_excel('data/TableS1.xlsx', header= 1)
AA_L = 10
os.system('mkdir Primer')
os.system('mkdir blastDB')
os.system('mkdir result/blast')
os.system('mkdir result/segs_id')

for group in Group:
    Fasta_DB = "result/" + group
    Seq_dict = read_F2D(Fasta_DB) # read the fasta
    Seq_dict = Dict_update(Seq_dict, TB) #Update Seq_dict by adding the seq from the table 
    # Check_quality(Seq_dict)
    '''
    for id in Seq_dict.keys():
        i_CDRL3 = 65#Seq_dict[id]["AA"].find(Seq_dict[id]["CDRL3"])
        i_CDRH1 = 130#Seq_dict[id]["AA"].find(Seq_dict[id]["CDRH1"])
        i_CDRH3 = 190#Seq_dict[id]["AA"].find(Seq_dict[id]["CDRH3"])
        Seq_fragment1 = Seq_dict[id]["AA"][         :i_CDRL3 + AA_L]
        Seq_fragment2 = Seq_dict[id]["AA"][i_CDRL3  :i_CDRH1 + AA_L]
        Seq_fragment3 = Seq_dict[id]["AA"][i_CDRH1  :i_CDRH3 + AA_L]
        Seq_fragment4 = Seq_dict[id]["AA"][i_CDRH3  :]
    '''
    # From table to dictionary
    # select the overlap sequence
    F = open("Primer/" + group, 'w')
    F.close()
    F = open("Primer/" + group, 'a')
    for id in Seq_dict.keys():
        try:
            i_CDRL3 = Seq_dict[id]["AA"].find(Seq_dict[id]["CDRL3"]) #65
        except:
            i_CDRL3 = 88
        i_CDRH1 = Seq_dict[id]["AA"].find(Seq_dict[id]["CDRH1"]) #120
        i_CDRH3 = Seq_dict[id]["AA"].find(Seq_dict[id]["CDRH3"]) #180
        for i in range(30):
            Seq_fragment1 = Seq_dict[id]["Nr"][           :(i_CDRL3 + AA_L)*3-20+i]
            Seq_fragment2 = Seq_dict[id]["Nr"][i_CDRL3*3  :(i_CDRH1 + AA_L)*3-20+i]
            Seq_fragment3 = Seq_dict[id]["Nr"][i_CDRH1*3  :(i_CDRH3 + AA_L)*3-20+i]
            #Seq_fragment4 = Seq_dict[id]["Nr"][i_CDRH3*3  :]
            F.write(f">{id}:Overlap1:{i}\n{Seq_fragment1[-30:]}\n")
            F.write(f">{id}:Overlap2:{i}\n{Seq_fragment2[-30:]}\n")
            F.write(f">{id}:Overlap3:{i}\n{Seq_fragment3[-30:]}\n")
    F.close()
    os.system(f"makeblastdb -in result/{group} -dbtype nucl -parse_seqids -out blastDB/{group}")
    Blast_check()



