import pandas as pd
import gzip
import networkx as nx
import numpy as np
import os
from itertools import combinations
from collections import defaultdict

# Carregar os dados
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

## Criação do grafo bipartido doença-gene
# Inicializa grafo bipartido
B = nx.Graph()

# Adiciona nós de genes e doenças
B.add_nodes_from(df["geneId"].unique(), bipartite="gene")
B.add_nodes_from(df["diseaseId"].unique(), bipartite="disease")

# Adiciona arestas com peso
edges = [(row["geneId"], row["diseaseId"], {"weight": row["score"]}) for idx, row in df.iterrows()]
B.add_edges_from(edges)

## Projeção do grafo para doença-doença
# Seleção dos nós de doenças
doenca_nodes = [n for n, d in B.nodes(data=True) if d["bipartite"] == "disease"]

gene_to_doencas = defaultdict(list)

for gene, doenca, data in B.edges(data=True):
    if B.nodes[gene]["bipartite"] == "gene":
        gene_to_doencas[gene].append((doenca, data["weight"]))
    else:
        gene_to_doencas[doenca].append((gene, data["weight"]))

G_disease = nx.Graph()

for doencas in gene_to_doencas.values():
    for (d1, w1), (d2, w2) in combinations(doencas, 2):
        weight = min(w1, w2)  # mantém o mínimo entre os scores
        if G_disease.has_edge(d1, d2):
            # Mantém o menor mínimo entre todos os genes compartilhados
            G_disease[d1][d2]['weight'] = min(G_disease[d1][d2]['weight'], weight)
        else:
            G_disease.add_edge(d1, d2, weight=weight)
                        

# Inicio da implementação do TMFG
doencas = list(G_disease.nodes())
index_map = {d: i for i, d in enumerate(doencas)}
N = len(doencas)

# Inicializa matriz de pesos
W = np.zeros((N, N))

# Preenche a matriz com os pesos
for u, v, data in G_disease.edges(data=True):
    i, j = index_map[u], index_map[v]
    W[i, j] = data['weight']
    W[j, i] = data['weight']  # simétrico
    
def tmfg(W):
    N = W.shape[0]
    nodes = list(range(N))
    
    # Step 1: Escolher 4 nós iniciais (maior soma de pesos)
    degree = W.sum(axis=1)
    init_nodes = np.argsort(degree)[-4:]
    
    G = nx.Graph()
    G.add_nodes_from(nodes)
    
    # Conectar tetraedro inicial
    for i in range(4):
        for j in range(i+1, 4):
            G.add_edge(init_nodes[i], init_nodes[j], weight=W[init_nodes[i], init_nodes[j]])
    
    # Lista de faces (triângulos)
    faces = [tuple(sorted(tri)) for tri in [
        (init_nodes[0], init_nodes[1], init_nodes[2]),
        (init_nodes[0], init_nodes[1], init_nodes[3]),
        (init_nodes[0], init_nodes[2], init_nodes[3]),
        (init_nodes[1], init_nodes[2], init_nodes[3])
    ]]
    
    remaining_nodes = [n for n in nodes if n not in init_nodes]
    
    # Step 2: Adicionar nós restantes
    for u in remaining_nodes:
        max_weight = -np.inf
        best_face = None
        for f in faces:
            weight_sum = W[u, f[0]] + W[u, f[1]] + W[u, f[2]]
            if weight_sum > max_weight:
                max_weight = weight_sum
                best_face = f
        # Adicionar arestas ao novo nó
        for v in best_face:
            G.add_edge(u, v, weight=W[u, v])
        # Atualizar faces
        faces.remove(best_face)
        faces.extend([
            tuple(sorted([u, best_face[0], best_face[1]])),
            tuple(sorted([u, best_face[0], best_face[2]])),
            tuple(sorted([u, best_face[1], best_face[2]]))
        ])
        
    return G

G_tmfg = tmfg(W)

# Renomear nós para os nomes reais das doenças
mapping = {i: d for i, d in enumerate(doencas)}
G_tmfg = nx.relabel_nodes(G_tmfg, mapping)

# Salvar o grafo TMFG em GraphML
output_dir = 'graph'
os.makedirs(output_dir, exist_ok=True)
full_path = os.path.join(output_dir, 'tmfg_cache.graphml')
nx.write_graphml(G_tmfg, full_path)

