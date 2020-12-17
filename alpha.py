
import pandas as pd


kmerDf = pd.read_csv("8mers_count.csv")
kmerDf.columns = ['8mer', 'count']
kmerDf['Alpha1'] = kmerDf['8mer'].str.slice(0,7)
kmerDf['Alpha2'] = kmerDf['8mer'].str.slice(1,8)
kmerDf['Denominator'] = kmerDf['8mer'].str.slice(1,7)
for each in kmerDf['Alpha1']:
    #print(each)
    kmerDf['Alpha1_Count']= kmerDf.groupby('Alpha1')['Alpha1'].transform('count') 
    kmerDf['Alpha2_Count']= kmerDf.groupby('Alpha2')['Alpha2'].transform('count') 
    kmerDf['Denominator_Count']= kmerDf.groupby('Denominator')['Denominator'].transform('count')
    kmerDf['P_Alpha1']= kmerDf['Alpha1_Count'].div(441)
    kmerDf['P_Alpha2']= kmerDf['Alpha2_Count'].div(441)
    kmerDf['P_Denominator']= kmerDf['Denominator_Count'].div(441)
    kmerDf['Exp_P']= kmerDf['P_Alpha1']*kmerDf['P_Alpha2']/kmerDf['P_Denominator']
    kmerDf['Obs_P']=kmerDf['count'].div(441)
    kmerDf['Obs/Exp']=kmerDf['Obs_P']/kmerDf['Exp_P']
print(kmerDf)
kmerDf.to_csv("ObsExp.csv")
