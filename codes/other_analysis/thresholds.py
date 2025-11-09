import pandas as pd
import gzip
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from collections import defaultdict
import numpy as np
import os

output_dir = "/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/output_images"

# Carregar os dados
with gzip.open("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/disgenet/download/all_gene_disease_associations.txt.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t")

# Selecionar colunas relevantes
usable_columns = ["geneId", "geneSymbol", "geneName", "diseaseId", "diseaseName", "score"]
df = df[usable_columns]

# Remover entradas com "Unknown" em diseaseName
df = df[df["diseaseName"] != "Unknown"]

# Unificar geneId duplicado de SFPQ
df.loc[df["geneSymbol"] == "SFPQ", "geneId"] = df[df["geneSymbol"] == "SFPQ"]["geneId"].min()

# Dicionário de doenças duplicadas
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

# Padronizar os diseaseId
df["diseaseId"] = df.apply(
    lambda row: doencas_duplicadas[row["diseaseName"]] if row["diseaseName"] in doencas_duplicadas else row["diseaseId"],
    axis=1
)

# Criar grafo bipartido doença-gene
B = nx.Graph()
B.add_nodes_from(df["geneId"].unique(), bipartite="gene")
B.add_nodes_from(df["diseaseId"].unique(), bipartite="disease")
edges = [(row["geneId"], row["diseaseId"], {"weight": row["score"]}) for idx, row in df.iterrows()]
B.add_edges_from(edges)

# Agrupar genes -> lista de doenças associadas
gene_to_doencas = defaultdict(list)
for gene, doenca, data in B.edges(data=True):
    if B.nodes[gene]["bipartite"] == "gene":
        gene_to_doencas[gene].append((doenca, data["weight"]))
    else:
        gene_to_doencas[doenca].append((gene, data["weight"]))

# Listas para armazenar resultados
thresholds = np.arange(0.01, 0.26, 0.01).round(2)
num_nodes = []
num_edges = []

# Loop pelos thresholds
for threshold in thresholds:
    G_disease = nx.Graph()
    for doencas in gene_to_doencas.values():
        for (d1, w1), (d2, w2) in combinations(doencas, 2):
            weight = min(w1, w2)
            if weight >= threshold:
                if G_disease.has_edge(d1, d2):
                    G_disease[d1][d2]['weight'] = min(G_disease[d1][d2]['weight'], weight)
                else:
                    G_disease.add_edge(d1, d2, weight=weight)

    num_nodes.append(G_disease.number_of_nodes())
    num_edges.append(G_disease.number_of_edges())
    print(f"Threshold={threshold:.2f} | Nós={G_disease.number_of_nodes()} | Arestas={G_disease.number_of_edges()}")

# Plotar gráficos
plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(thresholds, num_nodes, marker="o")
plt.title("Número de nós por threshold")
plt.xlabel("Threshold")
plt.ylabel("Número de nós")

plt.subplot(1,2,2)
plt.plot(thresholds, num_edges, marker="o", color="orange")
plt.title("Número de arestas por threshold")
plt.xlabel("Threshold")
plt.ylabel("Número de arestas")

plt.tight_layout()
file_path = os.path.join(output_dir, 'distribuicao_graus_log.png')
plt.savefig(file_path)
print(f"\nGráfico salvo em: '{file_path}'")