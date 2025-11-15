#!/bin/bash
# === export.sh ===
# Envoie les changements locaux vers GitHub

set -e  # Arrête le script en cas d'erreur

echo "🔍 Vérification du dépôt..."
if [ ! -d ".git" ]; then
  echo "❌ Ce dossier n'est pas un dépôt Git."
  exit 1
fi

echo "📦 Ajout des fichiers modifiés..."
git add .

# Message automatique (ou interactif)
if [ -n "$1" ]; then
  msg="$1"
else
  msg="Auto-commit from $(hostname) on $(date '+%Y-%m-%d %H:%M')"
fi

echo "📝 Commit : $msg"
git commit -m "$msg" || echo "ℹ️ Rien à committer (aucun changement)."

echo "Récupération des changements distants..."
git pull --rebase origin main || {
  echo "⚠️  Conflit détecté : résolvez les conflits puis relancez export.sh"
  exit 1
}

echo "🚀 Envoi vers GitHub..."
git push origin main && echo "✅ Export réussi."

