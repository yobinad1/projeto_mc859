import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os

output_dir = "/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/output_images"
os.makedirs(output_dir, exist_ok=True) 

G_tmfg = nx.read_graphml("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/graph/tmfg_cache.graphml")

print("\n--- Informações do Grafo TMFG ---")
is_planar, cert = nx.check_planarity(G_tmfg)
print("Planar?", is_planar)

num_vertices = G_tmfg.number_of_nodes()
num_arestas = G_tmfg.number_of_edges()

grau_medio = (2 * num_arestas) / num_vertices

print(f"Número de vértices: {num_vertices}")
print(f"Número de arestas: {num_arestas}")
print(f"Grau médio dos vértices: {grau_medio:.2f}")

# Número de componentes conexas
num_componentes = nx.number_connected_components(G_tmfg)
print(f"Número de componentes conexas: {num_componentes}")

# Distribuição de graus
degrees_dict = dict(G_tmfg.degree())
degrees = list(degrees_dict.values())


# Plotar a distribuição de graus
plt.figure(figsize=(12, 7))
plt.hist(degrees, bins=50, color='skyblue', edgecolor='black')
plt.title('Distribuição de Graus da Rede de Doenças (TMFG)')
plt.xlabel('Grau (Número de Conexões)')
plt.ylabel('Frequência (Número de Doenças)')
plt.grid(axis='y', alpha=0.75)

max_degree = max(degrees)
ticks = np.arange(0, max_degree + 50, 50) 
plt.xticks(ticks, rotation=0)

file_path = os.path.join(output_dir, 'distribuicao_graus.png')
plt.savefig(file_path)
print(f"\nGráfico salvo em: '{file_path}'")


# Plotar a distribuição de graus em escala logarítmica
plt.figure(figsize=(12, 7))
plt.hist(degrees, bins=50, color='skyblue', edgecolor='black')
plt.title('Distribuição de Graus da Rede de Doenças (Escala Log)')
plt.xlabel('Grau (Número de Conexões)')
plt.ylabel('Frequência (Número de Doenças) - Escala Log')
plt.grid(axis='y', alpha=0.75)
plt.yscale('log')

max_degree = max(degrees)
ticks = np.arange(0, max_degree + 50, 50) 
plt.xticks(ticks, rotation=0)

file_path = os.path.join(output_dir, 'distribuicao_graus_log.png')
plt.savefig(file_path)
print(f"\nGráfico salvo em: '{file_path}'")

degree_counts = np.bincount(degrees)  # índice = grau, valor = frequência
degree_values = np.nonzero(degree_counts)[0]  # apenas graus existentes
frequency = degree_counts[degree_values]

plt.figure(figsize=(12, 7))
plt.loglog(degree_values, frequency, marker='o', linestyle='None', color='red', markersize=2)
plt.title('Distribuição de Graus da Rede de Doenças (Log-Log)')
plt.xlabel('Grau (Número de Conexões)')
plt.ylabel('Frequência (Número de Doenças)')
plt.grid(True, which="both", ls="--", alpha=0.6)

file_path = os.path.join(output_dir, 'distribuicao_graus_loglog.png')
plt.savefig(file_path)
print(f"\nGráfico salvo em: '{file_path}'")