# -*- coding: utf-8 -*-
# ============================================================
# Université de Batna 2
# TP M1 Bioinformatique - Automatisation & Visualisation GEO
# ============================================================

import GEOparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
from datetime import datetime

# =======================
# CONFIGURATION DU TP
# =======================
geo_id = "GSE11121"  # Jeu de données stable et léger
BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
LOG_FILE = os.path.join(BASE_DIR, "logs", "tp_geo_log.txt")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# =======================
# FONCTION LOG
# =======================
def log(msg):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    message = f"[{timestamp}] {msg}"
    print(message)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(message + "\n")

# =======================
# ÉTAPE 1 — Téléchargement GEO
# =======================
log(f"🔽 Téléchargement du dataset {geo_id} depuis NCBI GEO...")
gse = GEOparse.get_GEO(geo=geo_id, destdir=DATA_DIR)
log("✅ Téléchargement terminé.")

# =======================
# ÉTAPE 2 — Métadonnées
# =======================
title = gse.metadata.get("title", ["Titre non disponible"])[0]
log(f"Titre de l'étude : {title}")
log(f"Nombre d'échantillons : {len(gse.gsms)}")

# =======================
# ÉTAPE 3 — Extraction des données d'expression
# =======================
log("📊 Extraction et fusion des données d’expression...")

tables = []
for gsm_name, gsm in gse.gsms.items():
    df = gsm.table.copy()

    # Détection automatique de la colonne d'identifiant
    id_col = None
    for candidate in ["ID_REF", "ID", "Gene", "GENE_SYMBOL", "SPOT_ID"]:
        if candidate in df.columns:
            id_col = candidate
            break

    if id_col is None or "VALUE" not in df.columns:
        log(f"⚠️ Échantillon ignoré (colonnes manquantes) : {gsm_name}")
        continue

    df = df[[id_col, "VALUE"]].set_index(id_col)
    df.rename(columns={"VALUE": gsm_name}, inplace=True)
    tables.append(df)

if not tables:
    raise ValueError("Aucune table valide trouvée — vérifiez le GEO ID.")

expression_data = pd.concat(tables, axis=1).dropna()
expression_data["mean_expression"] = expression_data.mean(axis=1)
log(f"✅ Fusion réussie : {expression_data.shape[0]} gènes, {expression_data.shape[1]-1} échantillons.")

# =======================
# ÉTAPE 4 — Statistiques
# =======================
log("📈 Calcul des statistiques descriptives...")
stats = expression_data["mean_expression"].describe()
stats.to_csv(os.path.join(RESULTS_DIR, "summary_stats.txt"), sep="\t")
log("✅ Statistiques enregistrées dans results/summary_stats.txt")

# =======================
# ÉTAPE 5 — Histogramme
# =======================
plt.figure(figsize=(10, 5))
sns.histplot(expression_data["mean_expression"], bins=50, kde=True, color="#1f77b4")
plt.title("Distribution des niveaux moyens d’expression génique")
plt.xlabel("Expression moyenne (log2)")
plt.ylabel("Fréquence")
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "hist_expression.png"))
plt.close()
log("📊 Histogramme sauvegardé dans results/hist_expression.png")

# =======================
# ÉTAPE 6 — Heatmap + Clustering
# =======================
log("🧠 Génération d'une heatmap avec clustering...")

subset = expression_data.drop(columns=["mean_expression"]).sample(n=min(50, len(expression_data)), random_state=42)
sns.clustermap(
    subset,
    cmap="coolwarm",
    metric="euclidean",
    method="average",
    figsize=(10, 10)
)
plt.savefig(os.path.join(RESULTS_DIR, "heatmap_clustering.png"))
plt.close()
log("🔥 Heatmap sauvegardée dans results/heatmap_clustering.png")

# =======================
# ÉTAPE 7 — Exportation finale
# =======================
expression_data.to_csv(os.path.join(RESULTS_DIR, "expression_cleaned.csv"))
log("✅ Données exportées dans results/expression_cleaned.csv")
log("🎉 TP exécuté avec succès ! Consultez le dossier results/.")
