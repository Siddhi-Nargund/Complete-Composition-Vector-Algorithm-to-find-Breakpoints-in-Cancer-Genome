from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
count = 0
total_len = 0
df = pd.DataFrame(columns=['Sequence', 'Length'])
with open("/home/snargund/CCV/gem_normal.R1.fastq") as in_handle:
    for title, seq, qual in FastqGeneralIterator(in_handle):
        df.loc[count] = [seq, len(seq)]
        count += 1
        total_len += len(seq)
        print(seq)
print("Reading file")
#print("%i length %i" % (count, total_len))
print(df)
df.to_csv('output_biopython.csv', header= 'TRUE')
