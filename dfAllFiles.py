≈#!/usr/bin/python3
#!/usr/bin/env python3

import pandas as pd
import os 

#
sourcedir = 'tumor_files'
dl = os.listdir('tumor_files')
pw= "tumor_files/"
print(pw)
dfm = pd.DataFrame(columns = ['Sequence', 'Length'])
for Filesindl in dl:
    frames2 = []
    step2 = 4
    #print(type(TumorFile))
    #next(TumorFile)
    Filesindl = pw + Filesindl
    #print(Filesindl)
    with open(Filesindl) as TumorFile:
        next(TumorFile)
        for lineno2, line2 in enumerate(TumorFile):
            tempframe2= []
            if lineno2 % step2 == 0:
                line2=line2.strip()
                end2= len(line2)
                # print(end)
                #print(TumorFile)
                tempframe2.append(line2)
                tempframe2.append(end2)
                frames2.append(tempframe2)
                #print(frames2)
    dfcurr = pd.DataFrame(frames2, columns= ['Sequence', 'Length'])
    # print(dfcurr)
    dfm=dfm.append(dfcurr, ignore_index = True)
print(dfm)
    
dfm.to_csv('sequence.csv', header= 'TRUE')
