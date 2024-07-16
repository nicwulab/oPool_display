'''
This function is designed for assign random codon to the amino acid to reducing the similarities among genes. 
'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','-I','--input')     
parser.add_argument('-p','-P','--pool', type = int)     
parser.add_argument('-n','-N','--negative')     
parser.add_argument('-o','-U','--output')     

args = parser.parse_args()
INPUT     = args.input
Pool_size = args.pool
NEGTIVE   = args.negative
OUTPUT    = args.output

'''
INPUT  = 'result/2024_0228_filtered.csv'
Pool_size = 2000000
NEGTIVE = 'result/random_neg.csv'
OUTPUT = 'result/TableS1_filtered.fa'
'''

import itertools, random, time, signal, multiprocessing
from multiprocessing import Pool, cpu_count 
from rich.progress import Progress, track
from collections import Counter
from Bio import AlignIO
import pandas as pd
import numpy as np




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

def VL_update(Primer_L, Primer_LAA, i, Extend = False):
    # add the length of the primer and re-select the condom
    len_n = len(Primer_L[i])/3
    if Extend:
        len_n +=1
    NS = TB.VL_AA.iloc[i].find(Primer_LAA[i])
    AA_new = TB.VL_AA.iloc[i][NS:int(NS+len_n)]
    primer_new =  "".join([random.choice(C_Table_R[i]) for i in AA_new])
    Primer_L[i] = primer_new
    Primer_LAA[i] = AA_new
    return Primer_L, Primer_LAA

def VH_update(Primer_H, Primer_HAA, i, Extend = False):
    # add the length of the primer and re-select the condom
    len_n = len(Primer_H[i])/3
    if Extend:
        len_n +=1
    NS = TB.VH_AA.iloc[i].find(Primer_HAA[i])
    AA_new = TB.VH_AA.iloc[i][NS:int(NS+len_n)]
    primer_new =  "".join([random.choice(C_Table_R[i]) for i in AA_new])
    Primer_H[i] = primer_new
    Primer_HAA[i] = AA_new
    return Primer_H, Primer_HAA

def TB_clean(TB):
    TB = TB[~TB.VL_nuc.isna()]
    TB = TB[~TB.VH_nuc.isna()]
    '''
    TB = TB[~TB.CDRL1_AA.isna()]
    TB = TB[~TB.CDRL2_AA.isna()]
    TB = TB[~TB.CDRL3_AA.isna()]
    TB = TB[~TB.CDRH1_AA.isna()]
    TB = TB[~TB.CDRH2_AA.isna()]
    TB = TB[~TB.CDRH3_AA.isna()]
    '''
    TB = TB[["X" not in i for i in  TB.VL_AA]]
    TB = TB[["X" not in i for i in  TB.VH_AA]]
    TB = TB[["*" not in i for i in  TB.VL_AA]]
    TB = TB[["*" not in i for i in  TB.VH_AA]]
    return TB

def AA2DNA_random(Seq):
    return "".join([random.choice(C_Table_R[i]) for i in Seq])

def list_multy(numbers):
    product = 1
    for number in numbers:
        product *= number
    return product

def Sim_Check(Primer_H, Thre = 0):
    MAX = len(Primer_H)
    Sim_lst = []
    duplicates = []
    for i in range(MAX-1):
        for ii in range(i+1, MAX):
            tmp = []
            for s1,s2 in zip(Primer_H[i],Primer_H[ii]):
                tmp +=[s1==s2]
            R_dic = dict(Counter(tmp))
            try:
                if R_dic[False]/len(tmp) <= Thre:
                    duplicates +=[i]
                    Sim_lst += [[i, ii]]
            except:
                duplicates +=[i]
                Sim_lst += [[i, ii]]
    duplicates = list(set(duplicates))
    print("similar seqs:", len(duplicates))
    return duplicates, Sim_lst

def Sim_Check_single(i, primer_new, Primer_base, Thre = 0.15):
    duplicates = []
    #Sim_lst = []
    for index in range(len(Primer_base)):
        tmp = []
        for s1,s2 in zip(Primer_base[index],primer_new):
            tmp +=[s1==s2]
        R_dic = dict(Counter(tmp))
        try:
            if R_dic[False]/len(tmp) <= Thre:
                duplicates +=[index]
                #Sim_lst += [R_dic[False]/len(tmp)]
        except:
            duplicates +=[index]
            #Sim_lst += [0]
    duplicates = list(set(duplicates))
    #print("similar seqs:", len(duplicates))
    return duplicates

def VH_update2(Primer_H, Primer_HAA, i, Thre = .15, Extend = False):
    # add the length of the primer and re-select the condom
    len_n = len(Primer_H[i])/3
    if Extend:
        len_n +=1
    NS = TB.VH_AA.iloc[i].find(Primer_HAA[i])
    AA_new = TB.VH_AA.iloc[i][NS:int(NS+len_n)]
    primer_new = Primer_H[i]
    dup_lst = Sim_Check_single(i, primer_new, Primer_H, Thre)
    primer_all = list(itertools.product(*[C_Table_R[i] for i in AA_new]))
    primer_all = ["".join(i) for i in primer_all]
    Primer_H = Primer_H[:i] + Primer_H[i+1:]
    Primer_base = [i for i in Primer_H if i != '']
    #######################
    ## Thanks for GPT4
    #######################
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    # List to store the results of the parallel execution
    results = []
    # Distribute the execution of Sim_Check_single across multiple processes
    for primer_new in primer_all:
        result = pool.apply_async(process_task, args=(primer_new, i, Primer_base, Thre))
        results.append(result)
    # Iterate through the results to find the first valid primer
    primer = None
    for result in results:
        primer_new = result.get()
        if primer_new is not None:
            primer = primer_new
            break
    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()
    if primer == None:
        raise ValueError
    else:
        Primer_H[i] = primer
        Primer_HAA[i] = AA_new
        return Primer_H, Primer_HAA

def Seq_update2(Primer_H, Primer_HAA, i, Thre = .15):
    # add the length of the primer and re-select the condom
    primer_new = Primer_H[i]
    primer_all = list(itertools.product(*[C_Table_R[i] for i in Primer_HAA[i]]))
    primer_all = ["".join(i) for i in primer_all]
    Primer_H = Primer_H[:i] + Primer_H[i+1:]
    Primer_base = [i for i in Primer_H if i != '']
    #######################
    ## Thanks for GPT4
    #######################
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    # List to store the results of the parallel execution
    results = []
    # Distribute the execution of Sim_Check_single across multiple processes
    for primer_new in primer_all:
        result = pool.apply_async(process_task, args=(primer_new, i, Primer_base, Thre))
        results.append(result)
    # Iterate through the results to find the first valid primer
    primer = None
    for result in results:
        primer_new = result.get()
        if primer_new is not None:
            primer = primer_new
            break
    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()
    if primer == None:
        raise ValueError
    else:
        Primer_H[i] = primer
        return Primer_H

# Function to be executed in parallel
def process_task(primer_new, i, Primer_base, Thre):
    if primer_new not in Primer_base:
        dup_lst = Sim_Check_single(i, primer_new, Primer_base, Thre)
        if len(dup_lst) == 0:
            return primer_new
    return None

def Primer_check_update(Primer_L, Primer_LAA, Thre, TXT):
    A = time.time()
    Thre = 0.40
    Failed_listL = []
    duplicates, Sim_lst = Sim_Check(Primer_L, Thre)
    numbers = range(len(duplicates))
    for N in track(numbers, description=f"Processing {TXT} ..."):
    #for i in [ii for ii in Rank_L_ord if ii in duplicates]:
        i = [ii for ii in Rank_L_ord if ii in duplicates][N]
        try:
            Primer_L = Seq_update2(Primer_L, Primer_LAA, i, Thre)
        except:
            Primer_L[i] = ""
            Failed_listL += [i]
    print(time.time() - A)
    return Primer_L

def Add_FR_region(TB):
    TB['FRL1_AA'] = ''
    TB['FRL2_AA'] = ''
    TB['FRL3_AA'] = ''
    TB['FRL4_AA'] = ''
    TB['FRH1_AA'] = ''
    TB['FRH2_AA'] = ''
    TB['FRH3_AA'] = ''
    TB['FRH4_AA'] = ''

    for i in range(len(TB)):
        CDRL1_NS = TB.VL_AA.iloc[i].find(TB.CDRL1_AA.iloc[i])
        CDRL1_NE = CDRL1_NS + len(TB.CDRL1_AA.iloc[i])
        CDRL2_NS = TB.VL_AA.iloc[i].find(TB.CDRL2_AA.iloc[i])
        CDRL2_NE = CDRL2_NS + len(TB.CDRL2_AA.iloc[i])
        CDRL3_NS = TB.VL_AA.iloc[i].find(TB.CDRL3_AA.iloc[i])
        CDRL3_NE = CDRL3_NS + len(TB.CDRL3_AA.iloc[i])

        CDRH1_NS = TB.VH_AA.iloc[i].find(TB.CDRH1_AA.iloc[i])
        CDRH1_NE = CDRH1_NS + len(TB.CDRH1_AA.iloc[i])
        CDRH2_NS = TB.VH_AA.iloc[i].find(TB.CDRH2_AA.iloc[i])
        CDRH2_NE = CDRH2_NS + len(TB.CDRH2_AA.iloc[i])
        CDRH3_NS = TB.VH_AA.iloc[i].find(TB.CDRH3_AA.iloc[i])
        CDRH3_NE = CDRH3_NS + len(TB.CDRH3_AA.iloc[i])

        TB.FRL1_AA.iloc[i] = TB.VL_AA.iloc[i][        :CDRL1_NS]
        TB.FRL2_AA.iloc[i] = TB.VL_AA.iloc[i][CDRL1_NE:CDRL2_NS]
        TB.FRL3_AA.iloc[i] = TB.VL_AA.iloc[i][CDRL2_NE:CDRL3_NS]
        TB.FRL4_AA.iloc[i] = TB.VL_AA.iloc[i][CDRL3_NE:]

        TB.FRH1_AA.iloc[i] = TB.VH_AA.iloc[i][        :CDRH1_NS]
        TB.FRH2_AA.iloc[i] = TB.VH_AA.iloc[i][CDRH1_NE:CDRH2_NS]
        TB.FRH3_AA.iloc[i] = TB.VH_AA.iloc[i][CDRH2_NE:CDRH3_NS]
        TB.FRH4_AA.iloc[i] = TB.VH_AA.iloc[i][CDRH3_NE:]

    return TB

def init_worker():
    # Ignore KeyboardInterrupt in the worker processes
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def process_seq(index):
    Seq = TB.VL_AA.iloc[index] + "GGGGS" * 3 + TB.VH_AA.iloc[index]
    temp_result = [AA2DNA_random(Seq) for _ in range(R)]
    local_result2 = []
    for ii, seq in enumerate(temp_result):
        Trunk = len(seq) // 99 + (len(seq) % 99 != 0)
        for N in range(Trunk):
            fragment = seq[N * 99:(N + 1) * 99]
            local_result2.append([f"{TB.Name.iloc[index]}:{N}-{ii}", fragment])
    return local_result2

# Create a pool of workers

print(time.ctime(), "Loading the data")
TB = pd.read_csv(INPUT, header= 0, index_col = 0)
# add negative control
tmp = pd.read_csv(NEGTIVE, index_col = 0)
TB = pd.concat([TB, tmp])
TB = TB_clean(TB)

C_Table_R = {}
for i in C_Table.values():
    C_Table_R.update({i:[]})
for i in C_Table.keys():
    C_Table_R[C_Table[i]] += [i]

C_Table_RC = {}
for i in C_Table_R:
    if i != '*':
        C_Table_RC.update({i:len(C_Table_R[i])})

print(time.ctime(), "Star the sampling...")
'''
R = int(Pool_size/len(TB))
Result = []
for i in track(range(len(TB)), description=f"Assign the triplet codon..."):
    Seq = TB.VL_AA.iloc[i] + "GGGGS"*3 + TB.VH_AA.iloc[i]
    Result +=[[AA2DNA_random(Seq) for i in range(R)]]

Result = np.array(Result)
Result2 = []
for i in track(range(len(Result)), description= ""):
    for ii in range(len(Result[i])):
        seq = Result[i,ii]
        Trunk = len(seq)//99
        if len(seq)%99 != 0:
            Trunk += 1
        for N in range(Trunk):
            Result2 += [[f"{TB.Name.iloc[i]}:{N}-{ii}", seq[N*99:(N + 1)*99]]]

print(time.ctime(), f"Done, the result is at {OUTPUT}")
Result2 = []
for i in track(range(len(TB)), description=f"Assigning the triplet codon..."):
    Seq = TB.VL_AA.iloc[i] + "GGGGS"*3 + TB.VH_AA.iloc[i]
    temp_result = [AA2DNA_random(Seq) for _ in range(R)]
    for ii in range(len(temp_result)):
        seq = temp_result[ii]
        Trunk = len(seq) // 99
        if len(seq) % 99 != 0:
            Trunk += 1
        for N in range(Trunk):
            fragment = seq[N*99:(N+1)*99]
            Result2.append([f"{TB.Name.iloc[i]}:{N}-{ii}", fragment])

'''
R = int(Pool_size / len(TB))

with Progress() as progress:
    num_processes = int(cpu_count() * .9)
    task1 = progress.add_task("[green]Processing...", total=len(TB))
    # Create a Pool of processes and ensure the init_worker function is used to ignore KeyboardInterrupt
    with Pool(processes=num_processes, initializer=init_worker) as pool:
        results = []
        # Iterate over results from pool.imap or pool.map
        for i, result in enumerate(pool.imap_unordered(process_seq, range(len(TB)))):
            results.extend(result)
            # Update the progress bar with each iteration
            progress.update(task1, advance=1)

print(time.ctime(), "Saving the result...")
#Result2 = np.array([item for sublist in results if sublist != 0 for item in sublist])
Result2 = np.array(results)
Seq_All = "\n".join(["\n".join([">"+i[0],i[1]]) for i in results]) + "\n" 

with open(OUTPUT, 'w') as F:
        F.write(Seq_All)
print(time.ctime(), f"Done, the result is at {OUTPUT}")
