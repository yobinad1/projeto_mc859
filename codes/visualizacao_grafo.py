import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os
import random
from networkx.algorithms.community import greedy_modularity_communities
from pyvis.network import Network

seed = 42
random.seed(seed)
np.random.seed(seed)

output_dir = os.path.join(os.path.dirname(__file__), '../output_images/visualizacao')
os.makedirs(output_dir, exist_ok=True)

# Carrega grafo
G_tmfg = nx.read_graphml("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/graph/tmfg_cache.graphml")

# Detecta comunidades
communities = list(greedy_modularity_communities(G_tmfg))

# Ordena comunidades por tamanho (maior -> menor)
communities_sorted = sorted(communities, key=len, reverse=True)

# Cria índice de comunidade para cada nó
comm_index = {}
for rank, comm in enumerate(communities_sorted):
    for node in comm:
        comm_index[node] = rank

# Cria grafo reduzido
G_comm = nx.Graph()
for rank, comm in enumerate(communities_sorted):
    G_comm.add_node(rank, size=len(comm))

for u, v in G_tmfg.edges():
    cu, cv = comm_index[u], comm_index[v]
    if cu != cv:
        if G_comm.has_edge(cu, cv):
            G_comm[cu][cv]['weight'] += 1
        else:
            G_comm.add_edge(cu, cv, weight=1)

# Visualização interativa
net = Network(notebook=False, height="1100px", width="100%", bgcolor="#ffffff", font_color="black")

max_w = max(d['weight'] for _,_,d in G_comm.edges(data=True))

for rank, comm in enumerate(communities_sorted):
    size = 10 + 5 * math.log(len(comm))  # escala de tamanho
    net.add_node(
        rank,
        label=f"{rank+1}",  # número da comunidade
        size=size,
        title=f"Comunidade {rank+1}, Tamanho: {len(comm)} doenças",
        )
    

for u, v, d in G_comm.edges(data=True):
    net.add_edge(u, v, value=d['weight']/max_w*10)

# Layout mais espaçado
net.repulsion(node_distance=300, spring_length=150, spring_strength=0.01, damping=0.9)

net.write_html(os.path.join(output_dir, "comunidades_supernos.html"))

n = len(G_comm.nodes)
adj = np.zeros((n, n))

for u, v, d in G_comm.edges(data=True):
    adj[u, v] = d['weight']
    adj[v, u] = d['weight']

plt.figure(figsize=(10, 8))
sns.heatmap(adj, cmap='YlGnBu', square=True, cbar_kws={'label': 'Peso da ligação'})
plt.title('Heatmap da matriz de adjacência entre comunidades')
plt.xlabel('Comunidade')
plt.ylabel('Comunidade')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "heatmap_comunidades.png"), dpi=300)
plt.close()

# Heatmap do grafo original, utilizando amostragem de arestas
frac = 0.15  # 15% das arestas
num_edges = int(frac * G_tmfg.number_of_edges())
edges_sample = random.sample(list(G_tmfg.edges()), num_edges)
subG = nx.Graph()
subG.add_edges_from(edges_sample)

if not nx.is_connected(subG):
    subG = subG.subgraph(max(nx.connected_components(subG), key=len)).copy()

# Métrica: clusterização
metric = nx.clustering(subG)
values = np.array(list(metric.values()))
norm = (values - values.min()) / (values.max() - values.min() + 1e-9)
node_colors = [plt.cm.plasma(norm[i]) for i in range(len(subG.nodes()))]

layouts = {
    "Kamada-Kawai": nx.kamada_kawai_layout(subG),
    "Fruchterman-Reingold": nx.spring_layout(subG, seed=seed, k=0.15, iterations=100)
}

for name, pos in layouts.items():
    plt.figure(figsize=(12, 9))
    nx.draw_networkx_nodes(subG, pos, node_color=node_colors, node_size=18, alpha=0.8)
    nx.draw_networkx_edges(subG, pos, alpha=0.08, width=0.5)
    plt.title(f"Subgrafo (15% das arestas) - Layout: {name}")
    plt.axis('off')
    sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=plt.Normalize(vmin=values.min(), vmax=values.max()))
    sm.set_array([])
    plt.colorbar(sm, label='Coeficiente de clusterização')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"subgrafo_heatmap_clusterizacao_{name.replace(' ', '_').lower()}.png"), dpi=300)
    plt.close()