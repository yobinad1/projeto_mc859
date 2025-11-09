import pandas as pd
import gzip

with gzip.open("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/disgenet/download/all_gene_disease_associations.txt.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t")

## Tratamento dos dados
# Selecionar colunas relevantes
usable_columns = ["geneId", "geneSymbol", "geneName", "diseaseId", "diseaseName", "score"]
df = df[usable_columns]

# Remover entradas com "Unknown" em diseaseName
df = df[df["diseaseName"] != "Unknown"]

# Unificando geneId de SFPQ que estava duplicado
df.loc[df["geneSymbol"] == "SFPQ", "geneId"] = df[df["geneSymbol"] == "SFPQ"]["geneId"].min()

# Unificar os ids de doenças que aparecem com mais de um id
doencas_duplicadas = {
    "Anaphylaxis": "umls:C0002792",
    "Appendicitis": "umls:C0003615",
    "Death": "umls:C0011065",
    "Ebstein's anomaly": "umls:C3665605",
    "Goldenhar Syndrome": "umls:C0265240",
    "Hyperkinetic conduct disorder": "umls:C0339004",
    "Hyperlipoproteinemia Type II": "umls:C0020445",
    "Hyperprolinemia": "umls:C0268529",
    "Hypertriglyceridemia": "umls:C0020557",
    "Mastitis": "umls:C0024894",
    "Mastocytoma": "umls:C0024897",
    "Meningioma": "umls:C0025286",
    "Mycetoma": "umls:C0024449",
    "Myofibroma": "umls:C1266121",
    "Nuclear cataract": "umls:C0392557",
    "islet cell tumor": "umls:C0242363"
}

# Substitui o diseaseId de todas as linhas pelo ID padronizado
df["diseaseId"] = df.apply(
    lambda row: doencas_duplicadas[row["diseaseName"]] if row["diseaseName"] in doencas_duplicadas else row["diseaseId"],
    axis=1
)

# Filtra pelo ID desejado
result = df[df["diseaseId"] == "umls:C0011853"]

# Mostra os nomes únicos encontrados
print(result["diseaseName"].unique())
