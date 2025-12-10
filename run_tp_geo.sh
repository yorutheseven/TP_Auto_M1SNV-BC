#!/bin/bash
# ===============================================
# TP_Auto_M1SNV - Script d'exécution automatique
# Université de Batna 2
# ===============================================

echo "🚀 Lancement du TP_Auto_M1SNV"
sleep 1

echo "📦 Installation des bibliothèques Python..."
pip install --quiet GEOparse pandas matplotlib seaborn biopython scipy

echo "▶️ Exécution du script principal..."
PYTHONIOENCODING=utf-8 python tp_geo_main.py

echo "✅ TP terminé. Consultez le dossier results/ pour les graphiques et fichiers."
