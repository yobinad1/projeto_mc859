import pandas as pd
import gzip

with gzip.open("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/disgenet/download/all_gene_disease_associations.txt.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t")

name_to_ids_count = df.groupby("diseaseName")["diseaseId"].nunique()

nomes_com_varios_ids = name_to_ids_count[name_to_ids_count > 1]

resultado = (
    df[df["diseaseName"].isin(nomes_com_varios_ids.index)]
    .groupby("diseaseName")["diseaseId"]
    .unique()
    .reset_index()
)

print(f"Total de nomes que aparecem com mais de um ID: {len(resultado)}")
print(resultado)

num_unknown = (df["diseaseName"] == "Unknown").sum()
print(f"Número de doenças 'Unknown': {num_unknown}")

num_doencas = df["diseaseName"].nunique()
print(f"Número de doenças únicas: {num_doencas}")