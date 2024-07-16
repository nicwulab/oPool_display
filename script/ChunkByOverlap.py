import pandas as pd
import numpy as np
import os
from Bio.Seq import Seq

def Fa2Dict(File):
    with open(File, 'r') as F:
        tmp = F.read().split('\n')[:-1]
    Seqs = {tmp[i*2].replace('>', ''):tmp[i*2 + 1] for i in range(len(tmp)//2)}
    return Seqs

def Segment_Seq(seq):
    seq_prim = [i for i in Primer_Slc if seq == i.split(':')[0]]
    seq_prim.sort()
    prim1 = Primer[seq_prim[0]]
    prim2 = Primer[seq_prim[1]]
    prim3 = Primer[seq_prim[2]]
    Range1 = range(Seqs[seq].find(prim1) + len(prim1))
    Range2 = range(Seqs[seq].find(prim1),
                   Seqs[seq].find(prim2) + len(prim2))
    Range3 = range(Seqs[seq].find(prim2),
                   Seqs[seq].find(prim3) + len(prim3))
    Range4 = range(Seqs[seq].find(prim3),
                   len(Seqs[seq]))
    Seg1 = "".join(np.array(list(Seqs[seq]))[Range1])
    Seg2 = "".join(np.array(list(Seqs[seq]))[Range2])
    Seg3 = "".join(np.array(list(Seqs[seq]))[Range3])
    Seg4 = "".join(np.array(list(Seqs[seq]))[Range4])
    return prim1.replace("U", "T"), prim2.replace("U", "T"), prim3.replace("U", "T"), Seg1.replace("U", "T"), Seg2.replace("U", "T"), Seg3.replace("U", "T"), Seg4.replace("U", "T")

Lst_Seq = os.listdir('Primer')

Result = []
for group in Lst_Seq:
    Seqs = Fa2Dict("result/"+group)
    Seqs_ID = list(Seqs.keys())
    Primer = Fa2Dict("Primer/"+group)
    Primer_TB = pd.read_csv("result/segs_id/"+ group + ".csv", index_col = 0)
    Primer_Slc = [":".join([str(ii) for ii in i]) for i in Primer_TB.to_numpy()]
    for seq in Seqs:
        prim1, prim2, prim3, Seg1, Seg2, Seg3, Seg4 = Segment_Seq(seq)
        Result += [{'Name':seq, 'group': group, 'Overlap1':prim1, 'Overlap2': prim2,
                    'Overlap3': prim3, 'Seg1': Seg1, 'Seg2': Seg2, 'Seg3': Seg3, 'Seg4': Seg4}]

TB = pd.read_csv('data/TableS1.csv', header = 1)
TB_seg = pd.DataFrame(Result)
TB = TB[TB.Name.isin(TB_seg.Name)]

# for storing the result
pd.merge(TB, TB_seg).to_csv('result/20240305_Fragmented.csv')

# list for the request
up_seqeunce   = "ACTAAAGGAGTATCCACCATG" #upstream sequence,   21bp
down_seqeunce = "TCCGGAGGATCCGATTACAAG" #downstream seqeunce, 21bp
TB_seg.group  = [i.replace('Re_assembled_','').replace('.fa','').replace('.', '') for i in TB_seg.group]

lst_all = []
for i in range(len(TB_seg)):
    lst_all += [[TB_seg.group.iloc[i], up_seqeunce + TB_seg.Seg1.iloc[i]]]
    lst_all += [[TB_seg.group.iloc[i], str(Seq(TB_seg.Seg2.iloc[i]).reverse_complement())]]
    lst_all += [[TB_seg.group.iloc[i], TB_seg.Seg3.iloc[i]]]
    lst_all += [[TB_seg.group.iloc[i], str(Seq(TB_seg.Seg4.iloc[i]).reverse_complement())]]


TB_req = pd.DataFrame(lst_all, columns = ['lib', 'seq'])

Seq_id = TB_seg.Name[TB_seg.group == '0_Sim_6'].to_list()

#Seqs = Fa2Dict("result/Re_assembled_0_Sim_.6.fa")

from Bio.Data import CodonTable
#from Bio.Alphabet import IUPAC
import random

def amino_acid_to_dna(amino_acid_sequence):
    dna_sequence = ""
    for amino_acid in amino_acid_sequence:
        # Find the codons that can encode the amino acid
        codons = [codon for codon, aa in standard_codon_table.forward_table.items() if aa == amino_acid]
        # If the amino acid is stop, add stop codons
        if amino_acid == "*":
            codons += standard_codon_table.stop_codons
        # Randomly select one of the codons
        selected_codon = random.choice(codons) if codons else None
        if selected_codon:
            dna_sequence += selected_codon
        else:
            print(f"No codon found for {amino_acid}")
    return dna_sequence

standard_codon_table = CodonTable.unambiguous_dna_by_id[1]

Random_seq = {id:amino_acid_to_dna(Seq(Seqs[id]).translate()) for id in Seqs.keys()}
Rand_all = []
for id in Random_seq.keys():
    Ran1 = range(len(TB_seg.Seg1[TB_seg.Name == id].iloc[0]))
    Ran2 = range(Seqs[id].replace("U", 'T').find(TB_seg.Seg2[TB_seg.Name == id].iloc[0]),
                 Seqs[id].replace("U", 'T').find(TB_seg.Seg2[TB_seg.Name == id].iloc[0]) + 
                 len(TB_seg.Seg2[TB_seg.Name == id].iloc[0]))
    Ran3 = range(Seqs[id].replace("U", 'T').find(TB_seg.Seg3[TB_seg.Name == id].iloc[0]),
                 Seqs[id].replace("U", 'T').find(TB_seg.Seg3[TB_seg.Name == id].iloc[0]) + 
                 len(TB_seg.Seg3[TB_seg.Name == id].iloc[0]))
    Ran4 = range(Seqs[id].replace("U", 'T').find(TB_seg.Seg4[TB_seg.Name == id].iloc[0]),
                 Seqs[id].replace("U", 'T').find(TB_seg.Seg4[TB_seg.Name == id].iloc[0]) +
                 len(TB_seg.Seg4[TB_seg.Name == id].iloc[0]))
    Rand_all += [
        up_seqeunce + "".join(np.array(list(Random_seq[id]))[Ran1]),
        str(Seq("".join(np.array(list(Random_seq[id]))[Ran2])).reverse_complement()),
        "".join(np.array(list(Random_seq[id]))[Ran3]),
        str(Seq("".join(np.array(list(Random_seq[id]))[Ran4]) + down_seqeunce).reverse_complement())]

TB_ex = pd.DataFrame([['Random', i] for i in Rand_all], columns  = ['lib', 'seq'])

TB_req2 = pd.concat([TB_req, TB_ex])

TB_req2.to_csv('result/test_request.csv')





