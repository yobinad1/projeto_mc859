import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms.community.quality import modularity

G_tmfg = nx.read_graphml("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/graph/tmfg_cache.graphml")
G_tmfg_dbht = nx.read_graphml("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/graph/tmfg_dbht_cache.graphml")


# --- TMFG puro ---
comunidades_tmfg = list(greedy_modularity_communities(G_tmfg))
mod_tmfg = modularity(G_tmfg, comunidades_tmfg)
print(f"TMFG: {len(comunidades_tmfg)} comunidades, modularidade = {mod_tmfg:.4f}")

# --- TMFG+DBHT ---
# Carregar clusters DBHT salvos em CSV
disease_clusters_dbht = {}
with open("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/graph/dbht_clusters.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        disease_clusters_dbht[row['diseaseId']] = int(row['cluster'])

# Agrupar nós por cluster
clusters_dbht = {}
for node, cid in disease_clusters_dbht.items():
    clusters_dbht.setdefault(cid, set()).add(node)
comunidades_dbht = list(clusters_dbht.values())

mod_dbht = modularity(G_tmfg_dbht, comunidades_dbht)
print(f"TMFG+DBHT: {len(comunidades_dbht)} comunidades, modularidade = {mod_dbht:.4f}")