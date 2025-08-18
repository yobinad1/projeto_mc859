import pandas as pd
import gzip

with gzip.open("/home/daniel/Documents/Unicamp/IC/MC859.2/disgenet/download/all_gene_disease_associations.txt.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t")
    
df.info()
#print(df.head())