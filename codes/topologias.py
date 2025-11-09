import networkx as nx
import matplotlib.pyplot as plt
import os
import pandas as pd
import gzip

with gzip.open("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/disgenet/download/all_gene_disease_associations.txt.gz", "rt") as f:
    df = pd.read_csv(f, sep="\t")

mapping = dict(zip(df["diseaseId"], df["diseaseName"]))
G_tmfg = nx.read_graphml("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/graph/tmfg_cache.graphml")

# Diretório de saída para salvar resultados
output_dir = os.path.join(os.path.dirname(__file__), '../output_images/topologias')
os.makedirs(output_dir, exist_ok=True)

def calcular_centralidades(G):
	print("\n--- Centralidades ---")
	grau = nx.degree_centrality(G)
	betweenness = nx.betweenness_centrality(G)
	closeness = nx.closeness_centrality(G)
	eigenvector = nx.eigenvector_centrality(G, max_iter=1000)
	
	print(f"Nó com maior centralidade de grau: {max(grau, key=grau.get)}")
	print(f"Nó com maior betweenness: {max(betweenness, key=betweenness.get)}")
	print(f"Nó com maior closeness: {max(closeness, key=closeness.get)}")
	print(f"Nó com maior eigenvector: {max(eigenvector, key=eigenvector.get)}")
	return grau, betweenness, closeness, eigenvector

def calcular_clusterizacao(G):
	print("\n--- Clusterização ---")
	clustering = nx.clustering(G)
	avg_clustering = nx.average_clustering(G)
	print(f"Coeficiente de clusterização médio: {avg_clustering:.4f}")
	return clustering, avg_clustering

def calcular_assortatividade(G):
	print("\n--- Assortatividade ---")
	assort = nx.degree_assortativity_coefficient(G)
	print(f"Coeficiente de assortatividade por grau: {assort:.4f}")
	return assort

def detectar_comunidades(G):
	print("\n--- Comunidades (Greedy Modularity) ---")
	from networkx.algorithms.community import greedy_modularity_communities
	communities = list(greedy_modularity_communities(G))
	print(f"Número de comunidades detectadas: {len(communities)}")
	sizes = [len(c) for c in communities]
	print(f"Tamanho das comunidades: {sizes}")
	return communities

def triadic_closure(G):
	print("\n--- Triadic Closure ---")
	triangles = nx.triangles(G)
	total_triangles = sum(triangles.values()) // 3
	transitivity = nx.transitivity(G)
	print(f"Número total de triângulos: {total_triangles}")
	print(f"Transitividade (razão de fechamento triádico): {transitivity:.4f}")
	return total_triangles, transitivity

def plot_top_centralidades(centralidade, nome, mapping, top_n=10):
    # Ordena valores
    top_nodes = sorted(centralidade.items(), key=lambda x: x[1], reverse=True)[:top_n]
    nodes, values = zip(*top_nodes)
    
    # Converte UMLS para nome real (se tiver mapping)
    node_labels = [mapping.get(n, n) for n in nodes]
    
    plt.figure(figsize=(10, 6))
    plt.barh(node_labels[::-1], values[::-1], color="tab:blue", alpha=0.7)
    plt.title(f"Top {top_n} - Centralidade de {nome}")
    plt.xlabel("Valor da centralidade")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"top_{nome.lower()}.png"), dpi=300)
    plt.close()
    print(f"Figura salva em {os.path.join(output_dir, f'top_{nome.lower()}.png')}")

#  --- Centralidades ---
grau, betweenness, closeness, eigenvector = calcular_centralidades(G_tmfg)

plot_top_centralidades(grau, "Grau", mapping)
plot_top_centralidades(betweenness, "Betweenness", mapping)
plot_top_centralidades(closeness, "Closeness", mapping)
plot_top_centralidades(eigenvector, "Eigenvector", mapping)

# --- Clusterização ---
clustering, avg_clustering = calcular_clusterizacao(G_tmfg)

# --- Análise de clusterização ---
# 1. Histograma dos coeficientes individuais
plt.figure(figsize=(8,5))
plt.hist(list(clustering.values()), bins=30, color='tab:blue', alpha=0.7)
plt.xlabel('Coeficiente de clusterização (nó)')
plt.ylabel('Frequência')
plt.title('Distribuição dos coeficientes de clusterização dos nós')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'hist_clusterizacao.png'), dpi=300)
plt.close()
print(f"\nFigura salva em {os.path.join(output_dir, 'hist_clusterizacao.png')}")

# 2. Comparação com rede aleatória
G_rand = nx.gnm_random_graph(G_tmfg.number_of_nodes(), G_tmfg.number_of_edges(), seed=42)
rand_clust = nx.average_clustering(G_rand)
print(f"\nCoeficiente de clusterização médio (rede aleatória): {rand_clust:.4f}")

# 3. Relação clusterização vs grau
graus = dict(G_tmfg.degree())
clust_vals = [clustering[n] for n in G_tmfg.nodes()]
grau_vals = [graus[n] for n in G_tmfg.nodes()]
plt.figure(figsize=(8,5))
plt.scatter(grau_vals, clust_vals, alpha=0.5)
plt.xlabel('Grau do nó')
plt.ylabel('Coeficiente de clusterização')
plt.title('Relação entre grau e clusterização dos nós')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'grau_vs_clusterizacao.png'), dpi=300)
plt.close()
print(f"\nFigura salva em {os.path.join(output_dir, 'grau_vs_clusterizacao.png')}")

# 4. Nós com maior e menor clusterização
top_clust = sorted(clustering.items(), key=lambda x: x[1], reverse=True)[:10]
low_clust = sorted(clustering.items(), key=lambda x: x[1])[:10]

print("\n--- Análise detalhada de clusterização ---")
print("Top 10 nós com maior clusterização:")
for n, v in top_clust:
	print(f"{mapping.get(n, n)}: {v:.3f}")
 
print("\n")
print("Top 10 nós com menor clusterização:")
for n, v in low_clust:
	print(f"{mapping.get(n, n)}: {v:.3f}")


# --- Assortatividade ---
assort = calcular_assortatividade(G_tmfg)

avg_neighbor_degree = nx.average_neighbor_degree(G_tmfg)

plt.scatter(list(dict(G_tmfg.degree()).values()), list(avg_neighbor_degree.values()), alpha=0.5)
plt.xlabel("Grau do nó")
plt.ylabel("Grau médio dos vizinhos")
plt.title("Relação grau x grau médio dos vizinhos")
plt.savefig(os.path.join(output_dir, "grau_vs_grau_medio_vizinhos.png"), dpi=300)
plt.close()


# --- Comunidades ---
communities = detectar_comunidades(G_tmfg)

# --- Triadic Closure ---
total_triangles, transitivity = triadic_closure(G_tmfg)

# --- Análise opcional de triadic closure ---
# 1. Comparação com rede aleatória
triangles_rand = nx.triangles(G_rand)
total_tri_rand = sum(triangles_rand.values()) // 3
transitivity_rand = nx.transitivity(G_rand)
print(f"\nTriadic closure (rede aleatória):")
print(f"Número total de triângulos: {total_tri_rand}")
print(f"Transitividade (rede aleatória): {transitivity_rand:.4f}")

# 2. Top 10 nós que participam de mais triângulos
triangles = nx.triangles(G_tmfg)
top_tri = sorted(triangles.items(), key=lambda x: x[1], reverse=True)[:10]
print("\nTop 10 nós que participam de mais triângulos:")
for n, v in top_tri:
	print(f"{mapping.get(n, n)}: {v}")

