#!/usr/bin/python3
#!/usr/bin/env python3
from __future__ import division
import itertools
import pandas as pd
from pandas import DataFrame
import numpy as np
pd.set_option("display.max_rows", None, "display.max_columns", None)
NormalDFData = []
def slidingWindow(sequenceStr, SeqLength, types, datafrs):
    i=int(0)
    j=int(8)
    #print(sequenceStr)
    temp = SeqLength
    finObyE = 0.0
    retList = []
    finSeqOByE = []
    while j<=SeqLength and i<=(SeqLength-8):
        tempSeqObyE = []
        s=sequenceStr[i:j]
        # print('**')
        # print(s)
        #retList.append(s)
        i += 1
        j += 1
        currKmers = GeneratingKmers(s)
        freqAll = CalculateFrequency(currKmers, datafrs)
        currIndividualProb = CalculateIndividualProbability(freqAll, datafrs)
        ExpectedProb = CalculateExpectedProbability(currIndividualProb)
        observedProb = CalculateObservedProbability(s, datafrs)
        oByE = CalculateObyEProbabilities(observedProb,ExpectedProb)
        tempSeqObyE.append(s)
        tempSeqObyE.append(oByE)
        # print('Done')
        types.append(tempSeqObyE)
        s = ''

    # tempVarName = 'df'+types
    tempVarName = pd.DataFrame(finSeqOByE, columns= ['KMER', 'oByE'])
    # print(tempVarName)
    # dfn = pd.concat([dfn,tempVarName])
    return finSeqOByE

def GeneratingKmers(sequenceStr):
    # print(sequenceStr)
    alphaOne= sequenceStr[:7]
    returnSolutionGeneratingKmers = {}
    # print("Alpha 1 is : "+str(alphaOne))
    returnSolutionGeneratingKmers['alphaOne'] = alphaOne
    alphaTwo=sequenceStr[1:]
    # print("Alpha 2 is : "+str(alphaTwo))
    returnSolutionGeneratingKmers['alphaTwo'] = alphaTwo
    denominator= sequenceStr[1:7]
    # print("Denominator is "+str(denominator))
    returnSolutionGeneratingKmers['denominator'] = denominator
    return returnSolutionGeneratingKmers

def CalculateFrequency(genKmer, datafrs):
    alphaOne = genKmer['alphaOne']
    alphaTwo = genKmer['alphaTwo']
    denominator = genKmer['denominator']
    calculatingFreqeuncy = {}
    fAlphaOne = datafrs.Sequence.str.contains(alphaOne).sum()
    calculatingFreqeuncy['frequencyAlphaOne'] = fAlphaOne
    #print("Frequency of Alpha 1 :"+str(fAlphaOne))
    fAlphaTwo =datafrs.Sequence.str.contains(alphaTwo).sum()
    calculatingFreqeuncy['frequencyAlphaTwo'] = fAlphaTwo
    fDenominator = datafrs.Sequence.str.contains(denominator).sum()
    calculatingFreqeuncy['frequencyDenominator'] = fDenominator
    #print("Frequency of denominator :"+str(fDenominator))
    return calculatingFreqeuncy

def CalculateIndividualProbability(freqAll, datafrs):
    #Need to change while making function
    LminusKplusOne = datafrs['Length'].sum(axis=0) - 9
    frequencyAlphaOne= freqAll['frequencyAlphaOne']
    frequencyAlphaTwo= freqAll['frequencyAlphaTwo']
    frequencyDenominator= freqAll['frequencyDenominator']
    calculateIndividualProbability= {}
    pAlphaOne= float((frequencyAlphaOne)/(LminusKplusOne))
    calculateIndividualProbability['ProbabilityOfAlphaOne']=pAlphaOne
    #print("Probalitity of alpha 1 is :"+str(pAlphaOne))
    pAlphaTwo=(frequencyAlphaTwo)/(LminusKplusOne)
    calculateIndividualProbability['ProbabilityOfAlphaTwo']=pAlphaTwo
    #print("Probalitity of alpha 2 is :"+str(pAlphaTwo))
    pDenominator=(frequencyDenominator)/(LminusKplusOne)
    calculateIndividualProbability['ProbabilityOfDenominator']=pDenominator
    #print("Probalitity of denominator is :"+str(pDenominator))
    return calculateIndividualProbability

def CalculateExpectedProbability(returnDictionaryFromCalculatingProbability):
    ProbAlphaOne= returnDictionaryFromCalculatingProbability['ProbabilityOfAlphaOne']
    ProbAlphaTwo= returnDictionaryFromCalculatingProbability['ProbabilityOfAlphaTwo']
    ProbDenominator= returnDictionaryFromCalculatingProbability['ProbabilityOfDenominator']
    if ProbDenominator==0:
        return 0
    else:
        E= ProbAlphaOne*ProbAlphaTwo/ProbDenominator
        #print(E)
        return E

#Calculating Observed Probability
def CalculateObservedProbability(eightMerSequence, datafrs):
    TotalLength = datafrs['Length'].sum(axis=0)
    ObservedFrequency= datafrs.Sequence.str.contains(eightMerSequence).sum()
    ObservedProb=ObservedFrequency/TotalLength
    #print("Observed Probability is :" +str(ObservedProb))
    return ObservedProb

#Calculating Observed by Expected Probabilities
def CalculateObyEProbabilities (ObservedProbability, ExpectedProbability):
    Result=ObservedProbability/ExpectedProbability
    return Result

frames = []
step = 4
completeData = ''
with open("gem_normal.R1.fastq") as NormalFile:
    next(NormalFile)
    for lineno, line in enumerate(NormalFile):
        tempframe= []
        if lineno % step == 0:
            line=line.strip()
            end= len(line)
            # print(end)Z
            #print(normalFile)
            tempframe.append(line)
            tempframe.append(end)
            frames.append(tempframe)
            #print(frames)
df = pd.DataFrame(frames, columns= ['Sequence', 'Length'])
# print(df)
normal = []
df['oByE'] = df.apply(lambda row: slidingWindow(row['Sequence'], row['Length'], normal, df),axis=1)
normalDF = pd.DataFrame(normal, columns = ['KMER', 'oByE'])
#print(normalDF)

# Tumor File Operations
'''
frames2 = []
step2 = 4
with open("gem_tumor.R1.fastq") as TumorFile:
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

TotalLength2 = df2['Length'].sum(axis=0)

LminusKplusOne2= TotalLength2 - 9
Tumor = []
df2['oByE'] = df2.apply(lambda row: slidingWindow(row['Sequence'], row['Length'], Tumor, df2),axis=1)
TumorDF = pd.DataFrame(Tumor, columns = ['KMER', 'oByE'])
#print(TumorDF)

matchedDF = pd.merge(normalDF, TumorDF, on='KMER', how='inner', indicator=True)
matchedDF = matchedDF.rename(columns={'oByE_x': 'Normal', 'oByE_y': 'Tumor'})
matchedDF = matchedDF.drop_duplicates().reset_index()
del matchedDF["_merge"]
matchedDF[matchedDF['Normal'] != matchedDF['Tumor']]
matchedDF = matchedDF.drop_duplicates().reset_index()
del matchedDF["_merge"]
matchedDF[matchedDF['Normal'] != matchedDF['Tumor']]
#print(matchedDF)

matchedDF.to_csv('result.csv', header= 'TRUE')

'''
