'''
1. Deleted unpaired antibodies
2. `AnnoG_Merge`
    1. Extract the nucl sequences and run igblast (pyir)
    2. Read the annotation result and retrieve the score of v gene
    3. Merge the table and sorting it based the V gene alignment score
    why should I care about v gene score?
3. `Table_clean`:
    1. Remove the common V gene family: IGHV1-69, IGHV6-1, and IGHV1-18.
    2. Remove the D gene family: IGHD3-9.
    3. Assign the id for unique clonotype (Clone type means the similarities of the antibody. The same clone type means they are similar. We just keep one sequence for each clone type.) 
    4. Exclude the clonotype 17.
    5. Keep the clonotype which are "HA:Ukn" only. (Because if the antibody was identified as the Stem antibody, we are not surprise other antibody from the same clone type are Stem-andtibody, too. So, here exclude them and only kept the "HA:Unk" and they are not similar to Other Stem antibody)
    6. Finally, we only kept one sequence from each clonotype
4. `Re_trans`:
    For somehow, some amino acid sequence-translation doesn't started in a correct position and containing stop codon. So, we try to translating them again when the sequences contain "*". After translation, all sequences end as TVSS or other similar sequences.
5. Antibody intact checking:
    For checking the intact of the amino acid, we aligned them into the Kabat number and checking if the first and the end aa are missing. 
6. In this function, we read the annotation result from the Pyir to retrieve the germlines id and get the germlines sequence from the database. And than, we fill the missing part if the germline sequence is complete.
7. Finally, we save the completed results
'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','-I','--input')     
parser.add_argument('-d','-D','--dlist', nargs='+')     
parser.add_argument('-v','-V','--vlist', nargs='+')
parser.add_argument('-g','-G','--germline')     
parser.add_argument('-o','-U','--output')     

args = parser.parse_args()
INPUT = args.input
D_VLIST = args.vlist
D_DLIST = args.dlist
GM_LOC  = args.germline
OUTPUT = args.output

''' Parameters for debug
INPUT = "data/TableS1.xlsx"
D_VLIST = ["IGHV1-69", "IGHV6-1", "IGHV1-18"]
D_DLIST = ["IGHD3-9"]
GM_LOC = "/home/wenkanl2/miniconda3/envs/Abs/lib/python3.9/site-packages/crowelab_pyir/data/germlines/Ig/human" 
OUTPUT = "result/2024_0228_filtered.csv"
'''

import os, time
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from abnumber import Chain
from rich.progress import track

import warnings
warnings.simplefilter(action='ignore')

def AnnoG_Merge(TB):
    '''
    This function is for:
        1. Extract the nucl sequences and run igblast (pyir)
        2. Read the annotation result and retrieve the score of v gene
        3. Merge the table and sorting it based the V gene alignment score
    PS: make sure the pyir is installed appropriately
    The intermediate file would be stored as "result/VH.a", "result/VL.a", "result/VH.a", 'result/VH.tsv.gz', and 'result/VL.tsv.gz'. With the exist of them, the annotation wouldn't start again. 
    '''
    with open("result/VH.fa", 'w') as F:
        F.write("\n".join([f">{i[0]}\n{i[1]}" for i in TB[['Name', 'VH_nuc']].to_numpy()]) + "\n")
    with open("result/VL.fa", 'w') as F:
        F.write("\n".join([f">{i[0]}\n{i[1]}" for i in TB[['Name', 'VL_nuc']].to_numpy()]) + "\n")
    # Germline annotation
    # read and integrate the result
    try:
        TB_anno = pd.read_csv('result/VH.tsv.gz', sep = '\t')
    except:
        os.system("pyir -m 50 result/VH.fa --outfmt tsv -o VH -s human")
        os.system("pyir -m 50 result/VL.fa --outfmt tsv -o VL -s human")
        os.system("mv VH.tsv.gz VL.tsv.gz result")
        TB_anno = pd.read_csv('result/VH.tsv.gz', sep = '\t')
    TB = pd.merge( TB, TB_anno[['sequence_id', 'v_score']], left_on = 'Name', right_on = 'sequence_id')
    TB.v_score = TB.v_score.fillna(0)
    TB = TB.sort_values("v_score", ascending = True)
    return TB

def Re_trans(TB):
    '''
    Try to translating them again when the sequences contain "*"
    After translation, all sequences end as TVSS or other similar sequences.
    '''
    List = TB.Name[["*" in i for i in TB.VH_AA]].to_list()
    for id in List:
        for i in range(3):
            STOP = "*" in str(Seq(TB.VH_nuc[TB.Name == id].iloc[0][i:]).translate())
            if STOP == False:
                break
        if STOP == False:
            TB.VH_AA[TB.Name == id] = str(Seq(TB.VH_nuc[TB.Name == id].iloc[0][i:]).translate())
    return TB

def Table_clean(TB, D_VLIST, D_DLIST):
    '''
    This function is for:
        1. Remove the common V gene family: IGHV1-69, IGHV6-1, and IGHV1-18.
        2. Remove the D gene family: IGHD3-9.
        3. Assign the id for unique clonotype (Clone type means the similarities of the antibody. The same clone type means they are similar. We just keep one sequence for each clone type.) 
        4. Exclude the clonotype 17.
        5. Keep the clonotype which are "HA:Ukn" only. (Because if the antibody was identified as the Stem antibody, we are not surprise other antibody from the same clone type are Stem-andtibody, too. So, here exclude them and only kept the "HA:Unk" and they are not similar to Other Stem antibody)
        6. Finally, we only kept one sequence from each clonotype
    '''
    # Filter Based on Heavy V Gene
    V_counts = [sum([ii in str(i) for ii in D_VLIST ]) for i in TB.Heavy_V_gene]
    TB = TB[np.array(V_counts)==0]
    # Exclude Specific D Gene
    D_counts = [sum([ii in str(i) for ii in D_DLIST ]) for i in TB.Heavy_D_gene]
    TB = TB[np.array(D_counts)==0]
    #TB = TB[["IGHD3-9" not in str(i)  for i in TB.Heavy_D_gene]]

    # Assign Unique ID to Missing Clonotypes
    TB.clonotype =  TB.clonotype.astype(str)
    TB.clonotype[TB.clonotype == 'nan'] = [f"uniq_{i}" for i in range(len(TB.clonotype[TB.clonotype == 'nan']))]

    # Filter Clonotypes Binding to Specific Antigen
    Ctype_lst = TB.clonotype.unique()[[list(set(TB['Binds to'][TB.clonotype == ctype].to_list())) == ['HA:Unk'] for ctype in TB.clonotype.unique()]]
    # Exclude Specific Clonotype
    TB = TB[TB.clonotype != "17.0"]
    TB = TB[TB.clonotype.isin(Ctype_lst)]

    # Remove Duplicate Clonotypes and save
    TB = TB[~TB.clonotype.duplicated()]
    #TB = TB[["X" not in i for i in  TB.VL_AA]]
    #TB = TB[["X" not in i for i in  TB.VH_AA]]
    #TB = TB[["*" not in i for i in  TB.VL_AA]]
    #TB = TB[["*" not in i for i in  TB.VH_AA]]
    return TB

def Kabt_Annot(TB, column = 'VH_AA'):
    '''
    kabat numbering the amino acid to find the missing sequences
    '''
    Kabat = {}
    Out = []
    #for i in range(len(TB[column])):
    for i in track(range(len(TB[column])), description=f"Kabat numbering {column}..."):
        try:
            Kabat.update({TB.Name.iloc[i]:dict(Chain(TB[column].iloc[i], scheme='kabat'))})
        except:
            Out += [TB.Name.iloc[i]]
    KABAT = pd.DataFrame(Kabat).T
    KABAT.columns = [str(i) for i in KABAT.columns]
    return KABAT, Out
    
def Seq_Complete(TB, chain, Lst1, Lst2):
    '''
    In this function, we read the annotation result from the Pyir to retrieve the germlines id and get the germlines sequence from the database.
    '''
    TB_anno = pd.read_csv(f'result/{chain}.tsv.gz', sep = '\t')
    for id in Lst1:
        #id = Lst1[0]
        tmp = TB_anno[TB_anno.sequence_id == id]
        v_call = tmp.v_call.iloc[0].split(',')[0]
        v_start = int(tmp.v_germline_start.iloc[0]//3)
        if v_start ==0 and tmp.v_germline_start.iloc[0] > 1:
            v_start += 1
        Filled_seq = str(Seq(Ref_V[v_call]).translate())[:v_start] + TB[TB.Name == id][f'{chain}_AA'].iloc[0]
        TB[f'{chain}_AA'][TB.Name == id] = Filled_seq
    for id in Lst2:
        tmp = TB_anno[TB_anno.sequence_id == id]
        j_call = tmp.j_call.iloc[0].split(',')[0]
        j_germ = tmp.j_germline_alignment_aa.iloc[0]
        J_germ = False
        for i in range(3):
            seq_g = str(Seq(Ref_J[j_call][i:]).translate())
            if tmp.j_germline_alignment_aa.iloc[0] in seq_g:
                J_germ = seq_g
                break
        Tail = J_germ[J_germ.find(j_germ)+len(j_germ):]
        TB[f'{chain}_AA'][TB.Name == id] += Tail
    return TB

def read_Ref(GM_LOC, G):
    file = f"{GM_LOC}/human_gl_{G}"
    with open(file, 'r') as F:
        Ref_V = {i.split('\n')[0]:"".join(i.split("\n")[1:]) for i in F.read().split(">")}
    #if G == 'V':
        #with open("/raid/home/wenkanl2/miniconda3/envs/Abs/lib/python3.9/site-packages/crowelab_pyir/data/germlines/Ig/human_VH/human_gl_V.fasta", 'r') as F:
        #    for i in F.read().split(">"):
        #        Ref_V.update({i.split('\n')[0]:"".join(i.split("\n")[1:])})
        #with open("/raid/home/wenkanl2/miniconda3/envs/Abs/lib/python3.9/site-packages/crowelab_pyir/data/germlines/Ig/human/human_gl_V.fasta", 'r') as F:
        #    for i in F.read().split(">"):
        #        Ref_V.update({i.split('\n')[0]:"".join(i.split("\n")[1:])})
    return Ref_V

# Load Data

print(time.ctime(), "Loading the data")
TB = pd.read_excel(INPUT, header = 1)
TB = TB[~TB.VH_AA.isna()]
TB = TB[~TB.VL_AA.isna()]

# Germline annotation 
print(time.ctime(), "Germline annotation")
TB = AnnoG_Merge(TB)
# Clean 
print(time.ctime(), "Sequences Cleaning")
TB = Table_clean(TB, D_VLIST, D_DLIST)
# Translation again
TB = Re_trans(TB)

######################################
## Kabat Annotation and head fill up
######################################


print(time.ctime(), "Kabat numbering and sequences completing")
Ref_V = read_Ref(GM_LOC, "V")
Ref_J = read_Ref(GM_LOC, "J")

# Kabat and missing list
KABAT, Out = Kabt_Annot(TB, column = 'VH_AA')
Lst1 = KABAT.index[KABAT.H1.isna()].to_list()
Lst2 = KABAT.index[KABAT.H113.isna()].to_list()
TB = Seq_Complete(TB, "VH", Lst1, Lst2)

KABAT, Out = Kabt_Annot(TB, column = 'VL_AA')
Lst1 = KABAT.index[KABAT.L1.isna()].to_list()
Lst2 = KABAT.index[KABAT.L107.isna()].to_list()
TB = Seq_Complete(TB, "VL", Lst1, Lst2)

#TB.VL_AA[TB.Name == Lst2[0]].iloc[0]
# fill the head 
print(time.ctime(), "Salving the result")
TB.to_csv(OUTPUT)

print(time.ctime(), f"Done, the result is at {OUTPUT}")
