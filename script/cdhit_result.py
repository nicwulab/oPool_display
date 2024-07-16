'''
This script is for reconnected the antibody sequences and create a potential PCR assembly overlap pool
'''

import argparse, os

parser = argparse.ArgumentParser()
parser.add_argument('-i','-I','--input')     
parser.add_argument('-n','-N','--negative')     
parser.add_argument('-nn','-NN','--nnegative', type = int)
parser.add_argument('-ng','-NG','--ngroup', type = int)     
parser.add_argument('-gs','-GS','--groupsize', type = int)     

args = parser.parse_args()
Fasta_DB     = args.input
NEGTIVE   = args.negative
N_neg = args.nnegative 
N_group      = args.ngroup
N_Seq = args.groupsize 


'''
Fasta_DB = "result/TableS1_filtered.fa"
NEGTIVE = 'result/random_neg.csv'
N_Seq = 25
N_group = 11
N_neg = 2
'''

import pandas as pd
import numpy as np
from collections import Counter

def cdhit2DB(File):
    # Initialize lists to store data
    sequence_names = []
    cluster_ids = []
    # Open the .clstr file and parse it
    with open(File, 'r') as file:
        current_cluster_id = None
        for line in file:
            if line.startswith('>Cluster '):
                current_cluster_id = line.split()[1]
            else:
                # Extracting sequence name
                parts = line.split()
                seq_name = parts[2].strip('>').split('...')[0]
                sequence_names.append(seq_name)
                cluster_ids.append(current_cluster_id)
    # Create a DataFrame
    df = pd.DataFrame({'Name': sequence_names, 
                       'ID': cluster_ids, 
                       "Sampling": [int(i.split(':')[-1].split('-')[1]) for i in sequence_names],
                       "Trunk": [int(i.split(':')[-1].split('-')[0]) for i in sequence_names] , 
                       "Name_Trunk": [i.split(':')[0]+":" + i.split(':')[-1].split('-')[0] for i in sequence_names]})
    return df

def CdBind(file):
    File = 'cdhit/' + file
    N_sim = file.split("0.")[-1].split('.')[0]
    df = cdhit2DB(File)
    df.columns  = ['Name', 'Sim_.' + N_sim, 'Sampling', 'Trunk', 'Name_Trunk']
    df = df.sort_values('Name')
    df = df.reset_index()
    return df


def Seq_grep(df_all, ID, N_Seq, Prefer_list, Neg_list):
    DF = df_all[~df_all[['Name_Trunk', ID ]].duplicated()]
    # Display the DataFrame
    Name_Trunk = df_all.Name_Trunk.unique()
    Name_count = dict(Counter(DF.Name_Trunk))
    Name_count = sorted(Name_count.items(), key=lambda x: x[1])
    Result_name = []
    Result_dispose = []
    while len(Name_count) != 0:
        # find the sequences from the prefer_list first
        if len(Prefer_list) != 0:
            name = Prefer_list[0] + ":0"
            Prefer_list.remove(Prefer_list[0])
        else:
            DF = DF[~DF.ID.isin(Neg_list)]
            Name_count = dict(Counter(DF.Name_Trunk))
            Name_count = sorted(Name_count.items(), key=lambda x: x[1])
            name = Name_count[0][0]
        name_cluster = [i[0] for i in Name_count if name.split(":")[0] == i[0].split(":")[0]]
        name_clusterTB = list(set([i  for i in Name_Trunk if name.split(":")[0] == i.split(":")[0]]))
        if len(name_cluster) == len(name_clusterTB):
            name_TB = DF[DF.Name_Trunk.isin(name_cluster)]
            # Counter the Cluster Number
            seq1 = name_TB.Name[~name_TB.Name_Trunk.duplicated()].to_list()
            if len(seq1) == len(name_clusterTB):
                Result_name += [seq1]
                DF = DF[~DF[ID].isin(name_TB[ID][~name_TB.Name_Trunk.duplicated()])]
                DF = DF[~DF.Name_Trunk.isin(name_cluster)]
                Name_count = dict(Counter(DF.Name_Trunk))
                Name_count = sorted(Name_count.items(), key=lambda x: x[1])
            else:
                Result_dispose += [name]
                DF = DF[DF.Name_Trunk != name]
                Name_count = dict(Counter(DF.Name_Trunk))
                Name_count = sorted(Name_count.items(), key=lambda x: x[1])
        else:
            Result_dispose += [name]
            DF = DF[~DF.Name_Trunk.isin(name_cluster)]
            Name_count = dict(Counter(DF.Name_Trunk))
            Name_count = sorted(Name_count.items(), key=lambda x: x[1])
        if len(Result_name) == N_Seq:
            #Result_name = np.array(Result_name)
            return Result_name
    return Result_name

def write_seq(df, Result_name, i, ID, Seq_dict):
    seqs = [id[0].split(":")[0] for id in Result_name]
    df = df[~df.ID.isin(seqs)]
    with open(f"result/Re_assembled_{i}_{ID}.fa", 'w') as F:
        F.close()
    with open(f"result/Re_assembled_{i}_{ID}.fa", 'a') as F:
        for seq_name in Result_name:
            seq_name.sort()
            F.write(">"+seq_name[0].split(':')[0] + "\n")
            F.write(''.join([Seq_dict[i] for i in seq_name]) + "\n")
    return df

def main():
    # combind all cd-hit cluster results
    df_lst = []
    List = [i for i in os.listdir('cdhit') if '.clstr' in i]
    for i in List:
        df_lst += [CdBind(i)]

    df_all = pd.concat([df_lst[0]] + [i.iloc[:,2] for i in df_lst[1:]] , axis = 1)
    df_all['ID'] = [i.split(":")[0] for i in df_all.Name]

    # read the negative list
    TB_neg = pd.read_csv(NEGTIVE)
    Neg_list = TB_neg.Name.to_list()[:26]

    ## Assembly

    import itertools
    from rich.progress import track

    with open(Fasta_DB, 'r') as F:
        Seqs = F.read().split('\n')

    #Name_list = list(itertools.chain.from_iterable(Result_name))
    Seq_dict = {}
    for i in track(range(int(len(Seqs)/2)), description="Turning Fasta into dict ..."):
        Seq_dict.update({Seqs[i*2].replace(">", ""):Seqs[i*2+1]})

    id = .6
    ID = f"Sim_.{str(id)[2:]}"
    i = -1
    while i < (N_group-1):
        i += 1
        Prefer_list = Neg_list[i*N_neg:(i+1)*N_neg]
        print(i, ID, Prefer_list)
        Result_name = Seq_grep(df_all, ID, N_Seq, Prefer_list, Neg_list)
        print(len(Result_name))
        if len(Result_name) == N_Seq:
            df_all = write_seq(df_all, Result_name, i, ID, Seq_dict)
        else:
            id += .05
            ID = f"Sim_.{str(id)[2:4].replace('0', '')}"
            i -= 1

    write_seq(df_all, Result_name, i, ID, Seq_dict)

if __name__ == "__main__":
    main()
