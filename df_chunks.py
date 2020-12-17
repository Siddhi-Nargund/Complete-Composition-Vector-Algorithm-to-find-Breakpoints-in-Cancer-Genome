import pandas as pd
#import psutil
from collections import defaultdict
import pickle

md = defaultdict()
kmerDf = pd.DataFrame()
def slidingWindow(sequenceStr, SeqLength):
    i=int(0)
    j=int(8)
    #print(sequenceStr)
    retList = []
    while j<=SeqLength and i<=(SeqLength-8):
        s=sequenceStr[i:j]
        i += 1
        j += 1
        if s in md.keys():
            md[s] += 1
        else:
            md[s] = 1
TotalLength = 0
for chunked_df in pd.read_csv("sequence.csv", chunksize=10000):
    # print(chunked_df)
    # print(type(chunked_df))
    for each in range(10000):
        print(chunked_df.iloc[each]['Sequence'])
        tempLen = chunked_df['Length']
        # print(type(chunked_df['Sequence']))
        slidingWindow(chunked_df.iloc[each]['Sequence'],chunked_df.iloc[each]['Length'])
        currSum = chunked_df['Length'].sum(axis=0)
        TotalLength += currSum
LminusKPlusOne = TotalLength - 9
TotalLengthString = "Total Length is: " +str(LminusKPlusOne)

file= open ('Total_Length.txt', "w")
file.write(TotalLengthString)
file.close()

#print(md)
with open('8mers_count.csv', 'w') as f:
    for key in md.keys():
        f.write("%s,%s\n"%(key,md[key]))

file_to_write = open("kmers.pickle", "wb")

# pickle.dump(md, file_to_write)

# dataf=pd.DataFrame(list(md.items()), columns= ['8mer','Count'])
# dataf.to_pickle('kmer.pkl')


# print("Final Dataframe: ")
# #print(kmerDf.head())
# # dups = kmerDf.pivot_table(index = ['8MER'], aggfunc ='size')
# # print(dups)
# kmerDf.to_pickle('kmers.pkl')
