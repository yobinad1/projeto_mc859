import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

def link_prediction_experiment_auc(
    G, 
    frac_remove=0.1, 
    method="adamic_adar", 
    num_negative=None, 
    seed=42
):
    """
    Experimento de Link Prediction com avaliação via ROC-AUC e PR-AUC.
    - Remove uma fração das arestas reais.
    - Gera candidatos positivos (arestas removidas) e negativos (pares não conectados).
    - Avalia com Adamic-Adar, RA, Jaccard ou Preferential Attachment.
    """

    random.seed(seed)
    G = G.copy()
    edges = list(G.edges())
    num_remove = int(frac_remove * len(edges))
    removed_edges = random.sample(edges, num_remove)

    # Remove arestas para o teste
    G.remove_edges_from(removed_edges)

    # --- Conjunto positivo ---
    positives = [(u, v) for u, v in removed_edges]

    # --- Conjunto negativo (amostra aleatória de não-arestas) ---
    non_edges = list(nx.non_edges(G))
    if num_negative is None:
        num_negative = len(positives)
    negatives = random.sample(non_edges, num_negative)

    # --- Escolhe o preditor ---
    if method == "adamic_adar":
        predictor = nx.adamic_adar_index(G, positives + negatives)
    elif method == "resource_allocation":
        predictor = nx.resource_allocation_index(G, positives + negatives)
    elif method == "jaccard":
        predictor = nx.jaccard_coefficient(G, positives + negatives)
    elif method == "preferential_attachment":
        predictor = nx.preferential_attachment(G, positives + negatives)
    else:
        raise ValueError("Método desconhecido")

    # --- Extrai scores ---
    y_true = [1] * len(positives) + [0] * len(negatives)
    y_score = []
    for _, _, p in predictor:
        y_score.append(p)

    # --- Métricas ---
    roc_auc = roc_auc_score(y_true, y_score)
    pr_auc = average_precision_score(y_true, y_score)

    return {
        "method": method,
        "frac_remove": frac_remove,
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "y_true": y_true,
        "y_score": y_score
    }

output_dir = "output_images/curvas_link_pred"
os.makedirs(output_dir, exist_ok=True)

seed = 42
random.seed(seed)
np.random.seed(seed)

G_tmfg = nx.read_graphml("/home/daniel/Documents/Unicamp/IC/MC859/projeto_mc859/graph/tmfg_cache.graphml")

methods = ["adamic_adar", "resource_allocation", "jaccard", "preferential_attachment"]

results = []
for method in methods:
    scores = []
    y_true_plot = None
    y_score_plot = None
    roc_auc_plot = None
    pr_auc_plot = None
    for seed in range(5):  # 5 rodadas com seeds diferentes
        res = link_prediction_experiment_auc(G_tmfg, frac_remove=0.1, method=method, seed=seed)
        scores.append((res["roc_auc"], res["pr_auc"]))
        # Salva os dados da primeira rodada para plotar
        if seed == 0:
            y_true_plot = res["y_true"]
            y_score_plot = res["y_score"]
            roc_auc_plot = res["roc_auc"]
            pr_auc_plot = res["pr_auc"]
    mean_roc = np.mean([s[0] for s in scores])
    mean_pr = np.mean([s[1] for s in scores])
    results.append((method, mean_roc, mean_pr))

    # Plot das curvas para a primeira seed
    fpr, tpr, _ = roc_curve(y_true_plot, y_score_plot)
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (area = {roc_auc_plot:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve - {method}')
    plt.legend(loc='lower right')
    plt.savefig(os.path.join(output_dir, f"roc_curve_{method}.png"))
    plt.close()

    precision, recall, _ = precision_recall_curve(y_true_plot, y_score_plot)
    plt.figure()
    plt.plot(recall, precision, label=f'PR curve (area = {pr_auc_plot:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve - {method}')
    plt.legend(loc='lower left')
    plt.savefig(os.path.join(output_dir, f"pr_curve_{method}.png"))
    plt.close()

print("\n=== Resultados Médios (5 rodadas, 10% removido) ===")
for method, roc, pr in results:
    print(f"{method:25s} | ROC-AUC: {roc:.4f} | PR-AUC: {pr:.4f}")
