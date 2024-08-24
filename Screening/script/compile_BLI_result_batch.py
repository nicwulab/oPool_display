import os
import glob

def ReadingData(filename, trial, datatype):
    DataDict = {}
    with open(filename, 'r') as infile:
        printing = 0
        for line in infile.readlines():
            if printing == 1:
                line = line.rstrip().rsplit("\t")
                if len(line) < 3: continue
                if trial == 'Time1': time = line[0]; sign = line[1]; fit = line[2]
                if trial == 'Time2': time = line[3]; sign = line[4]; fit = line[5]
                if trial == 'Time3': time = line[6]; sign = line[7]; fit = line[8]
                if trial == 'Time4': time = line[9]; sign = line[10]; fit = line[11]
                if trial == 'Time5': time = line[12]; sign = line[13]; fit = line[14]
                DataDict[time] = fit if datatype == 'fit' else sign
            if 'Time1' == line[0:5]: printing = 1
    return DataDict

def CompileData(AllData, outfile):
    with open(outfile, 'w') as outfile:
        SampleIDs = AllData.keys()
        outfile.write("\t".join(['Time', 'Signal', 'SampleID']) + "\n")
        for ID in SampleIDs:
            for time in AllData[ID].keys():
                outfile.write("\t".join([time, AllData[ID][time], ID]) + "\n")

def wrapper(In_Folder, Out_Folder, Exp, trial, datatype):
    filenames = glob.glob(In_Folder + '/' + Exp + '_*.txt')
    outfile = Out_Folder + '/' + Exp + '_' + datatype + '_All.compile'
    AllData = {}
    for filename in filenames:
        print("\tReading file: %s" % filename)
        ID = os.path.basename(filename).replace('.txt', '')
        AllData[ID] = ReadingData(filename, trial, datatype)
    CompileData(AllData, outfile)
    print("Compiled Data into %s" % outfile)

def main():
    In_Folder = 'data/BLI_data'
    Out_Folder = 'result'
    
    targets = ['H3-stem', 'H1-stem']
    H3_stem_binders = ['01.h.02-scFv', '2F02-scFv', '2F04-scFv', '31.a.55-scFv', '56.j.01-scFv', 
                       'AG2-G02-scFv', '01.ad.01-scFv','3C06-scFv','3F02-scFv',
                       '01.h.02-Fab', '2F02-Fab', '2F04-Fab', '31.a.55-Fab', '56.j.01-Fab', 
                       '01.ad.01-Fab', '3C06-Fab', '3F02-Fab']
    H1_stem_binders = ['16.ND.92-scFv', '2F01-scFv', '31.a.55-scFv', '3C06-scFv','3F02-scFv',
                       '16.ND.92-Fab', '2F01-Fab', '31.a.55-Fab', '3C06-Fab', '3F02-Fab']
    
    # Store binders in a dictionary
    binders_dict = {
        'H3-stem': H3_stem_binders,
        'H1-stem': H1_stem_binders
    }

    datatypes = ['fit', 'sign']
    Trials = ['Time1']

    for target in targets:
        binders = binders_dict[target]
        
        for binder in binders:
            Exp = binder + '_' + target
            for Trial in Trials:
                for datatype in datatypes:
                    print(f"Processing {Exp}")
                    wrapper(In_Folder, Out_Folder, Exp, Trial, datatype)

if __name__ == "__main__":
    main()