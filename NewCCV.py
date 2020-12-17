#!/usr/bin/python3
#!/usr/bin/env python3
from __future__ import division
import itertools
import pandas as pd
from pandas import DataFrame
#for normal
fileOpen= open("/home/snargund/CCV/gem_normal_files/gem_normal.R1.fa")
NormalFile=fileOpen.readlines()
# NormalFile = [line.rstrip() for line in open("/home/snargund/CCV/gem_normal_files/gem_normal.R1.fa")]
print(NormalFile)
#Converting file to open from list to string 
fileString = ' '.join([str(elem) for elem in NormalFile]) 
#print(fileString)  
TotalLength=(len(fileString))
#print("Lenght of the sequence is :"+str(TotalLength))


#s is global variable and sequenceSubString is local for the function
def GeneratingKmers(sequenceSubString):
    alphaOne= sequenceSubString[:7]
    returnSolutionGeneratingKmers = {}
    print("Alpha 1 is : "+str(alphaOne))
    returnSolutionGeneratingKmers['alphaOne'] = alphaOne
    alphaTwo=sequenceSubString[1:]
    print("Alpha 2 is : "+str(alphaTwo))
    returnSolutionGeneratingKmers['alphaTwo'] = alphaTwo
    denominator= sequenceSubString[1:7]
    print("Denominator is "+str(denominator))
    returnSolutionGeneratingKmers['denominator'] = denominator
    return returnSolutionGeneratingKmers

#Calculating frequency
def CalculateFrequency(returnDictionaryFromGeneratingKmers):
    alphaOne = returnDictionaryFromGeneratingKmers['alphaOne']
    alphaTwo = returnDictionaryFromGeneratingKmers['alphaTwo']
    denominator = returnDictionaryFromGeneratingKmers['denominator']
    calculatingFreqeuncy = {}
    fAlphaOne = (fileString.count(alphaOne))
    calculatingFreqeuncy['frequencyAlphaOne'] = fAlphaOne
    #print("Frequency of Alpha 1 :"+str(fAlphaOne))
    fAlphaTwo =(fileString.count(alphaTwo))
    calculatingFreqeuncy['frequencyAlphaTwo'] = fAlphaTwo
    #print("Frequency of Alpha 2 :"+str(fAlphaTwo))
    fDenominator = (fileString.count(denominator))
    calculatingFreqeuncy['frequencyDenominator'] = fDenominator
    #print("Frequency of denominator :"+str(fDenominator))
    return calculatingFreqeuncy

#Calculating probability
def CalculateIndividualProbability(returnDictionaryFromCalculatingFrequency):
    LminusKplusOne= TotalLength - 9     #Need to change while making function
    #print(LminusKplusOne)
    frequencyAlphaOne= returnDictionaryFromCalculatingFrequency['frequencyAlphaOne']
    frequencyAlphaTwo= returnDictionaryFromCalculatingFrequency['frequencyAlphaTwo']
    frequencyDenominator= returnDictionaryFromCalculatingFrequency['frequencyDenominator']
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

#Calculating Expected Probalilty by Markov model
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
def CalculateObservedProbability(eightMerSequence):
    ObservedFrequency= (fileString.count(s))
    ObservedProb=ObservedFrequency/TotalLength
    #print("Observed Probability is :" +str(ObservedProb))
    return ObservedProb

#Calculating Observed by Expected Probabilities
def CalculateObyEProbabilities (ObservedProbability, ExpectedProbability):
    Result=ObservedProbability/ExpectedProbability
    return Result

NormalDictArray = []
i=int(0)
j=int(8)
end=TotalLength
for seq in NormalFile:
    while j<=end and i<=(end-8):
        NormalDict= {}
        s=seq[i:j]
        NormalDict['8mer']=s
        #print(s)
        solDictionary = GeneratingKmers(s)
        #print(solDictionary)
        freqDictionary = CalculateFrequency(solDictionary)
        #print(freqDictionary)
        individualProbDictionary= CalculateIndividualProbability(freqDictionary)
        #print(individualProbDictionary)
        ObservedProbability=CalculateObservedProbability(s)
        #print(ObservedProbability)
        ExpectedProbability=CalculateExpectedProbability(individualProbDictionary)
        #print(ExpectedProbability)

        NormalFinalResult=CalculateObyEProbabilities(ObservedProbability,ExpectedProbability)
        #print("Observed Probability Divided by Expected Probability is:")
        #print(NormalFinalResult)
        NormalDict['O/E Prob']= NormalFinalResult
        NormalDictArray.append(NormalDict)
        #print(NormalDict)
        i+=1
        j+=1

#For tumor 

#print("*****TUMOR******")
fileOpen= open("/home/snargund/CCV/gem_normal_files/gem_tumor.R1.fa")
TumorFile=fileOpen.readlines()
# TumorFile = [line.rstrip() for line in open("sample_gem_tumor.R1.fa")]
print(TumorFile)
#Converting file to open from list to string 
fileString = ' '.join([str(elem) for elem in TumorFile]) 
#print(fileString)  
TotalLength=(len(fileString))
#print("Lenght of the sequence is :"+str(TotalLength))

TumorDictArray=[]
i=int(0)
j=int(8)
end=TotalLength
for seq in TumorFile:
    while j<=end and i<=(end-8):
        TumorDict= {}
        s=seq[i:j]
        TumorDict['8mer']=s
        #print(s)
        solDictionary = GeneratingKmers(s)
        #print(solDictionary)
        freqDictionary = CalculateFrequency(solDictionary)
        #print(freqDictionary)
        individualProbDictionary= CalculateIndividualProbability(freqDictionary)
        #print(individualProbDictionary)
        ObservedProbability=CalculateObservedProbability(s)
        #print(ObservedProbability)
        ExpectedProbability=CalculateExpectedProbability(individualProbDictionary)
        #print(ExpectedProbability)
        TumorFinalResult=CalculateObyEProbabilities(ObservedProbability,ExpectedProbability)
        #print("Observed Probability Divided by Expected Probability is:")
        #print(TumorFinalResult)
        TumorDict['O/E Prob']= TumorFinalResult
        TumorDictArray.append(TumorDict)
        #print(TumorDictArray)
    
        i+=1
        j+=1
    
#print(TumorDictArray)
#print(NormalDictArray)

#Temp to access the 8mer data since it is list of dicts
TempNormalDictArray = []
for b in NormalDictArray:
    TempNormalDictArray.append(b['8mer'])

df= pd.DataFrame()
#curr= current element in tumor dict
#index() gives index, try except for values which don;t have index
temp=[]
frames= []
for a in range(len(TumorDictArray)):
    curr = TumorDictArray[a]['8mer']
    
    try:
        finalDict= {}
        
        CurrIndex = TempNormalDictArray.index(curr)
        #print("Found "+str(curr)+" in Normal Dictionary at indices:"+str(CurrIndex))
        #print("Data in Tumor Dictionary is "+str(TumorDictArray[a]))
        #print("Data in Normal Dictionary at "+str(CurrIndex)+" is "+str(NormalDictArray[CurrIndex]))
        #print("Data in Normal Dictionary is "+str(NormalDictArray[CurrIndex]))
        temp=(TumorDictArray[a], NormalDictArray[CurrIndex])
        finalDict['Tumor']=TumorDictArray[a]
        finalDict['Normal']=NormalDictArray[CurrIndex]
        #print(finalDict)
        df=pd.DataFrame(finalDict)
        print(df)
        frames.append(df)

        
    except ValueError as ve:
        continue

result= pd.concat(frames)
result.to_csv('result.csv', header= 'TRUE')

