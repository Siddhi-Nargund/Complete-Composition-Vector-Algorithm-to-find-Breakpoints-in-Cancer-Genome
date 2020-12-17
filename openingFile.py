#!/usr/bin/python3
#!/usr/bin/env python3
from __future__ import division
import itertools
import pandas as pd
from pandas import DataFrame

frames2 = []
step2 = 4
with open("/home/snargund/CCV/tumor_files/sequenceaa") as TumorFile:
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
            #print(frames)
df2 = pd.DataFrame(frames2, columns= ['Sequence', 'Length'])
#print(df2)
df2.to_csv('sequence.csv', header= 'TRUE')
