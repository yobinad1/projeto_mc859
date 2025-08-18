import pandas as pd
import gzip

with gzip.open("/home/daniel/Documents/Unicamp/IC/MC859.2/disgenet/download/all_gene_disease_associations.txt.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t")
    
gene_consistency = df.groupby("geneSymbol")["geneId"].nunique()

genes_inconsistentes = gene_consistency[gene_consistency > 1]

if not genes_inconsistentes.empty:
    print("Existem genes inconsistentes:")
    for gene in genes_inconsistentes.index:
        gene_ids = df[df["geneSymbol"] == gene]["geneId"].unique()
        print(f"{gene}: {list(gene_ids)}")
else:
    print("Todos os geneSymbols estão associados a apenas um geneId.")
    
num_genes = df["geneSymbol"].nunique()
print(f"Número de genes únicos: {num_genes}")