#!/bin/bash
# === import.sh ===
# Récupère la dernière version depuis GitHub
# Usage :
#   ./import.sh        → mise à jour normale avec sauvegarde
#   ./import.sh -f     → forçage (ignore complètement les modifications locales)

set -e

# --- Détection du dépôt ---
if [ ! -d ".git" ]; then
  echo "❌ Ce dossier n'est pas un dépôt Git."
  exit 1
fi

# --- Détection automatique de la branche principale ---
main_branch=$(git remote show origin | sed -n '/HEAD branch/s/.*: //p')
[ -z "$main_branch" ] && main_branch="main"

# --- Option force ---
if [ "$1" = "-f" ]; then
  echo "Mode FORCÉ : les modifications locales seront écrasées."
  sleep 1
  git fetch origin
  git reset --hard origin/$main_branch
  git clean -fd
  echo "Import forcé terminé (toutes les modifications locales ont été supprimées)."
  exit 0
fi

# --- Mode normal ---
echo "🧹 Vérification des modifications locales..."
if [ -n "$(git status --porcelain)" ]; then
  backup_branch="backup_$(date '+%Y%m%d_%H%M%S')"
  echo "💾 Création d'une branche de sauvegarde : $backup_branch"
  git checkout -b "$backup_branch"
  echo "✅ Sauvegarde effectuée."
  git checkout "$main_branch"
fi

echo "🔄 Mise à jour depuis GitHub..."
git fetch origin
git pull --rebase origin "$main_branch" || {
  echo "Conflit détecté : corrigez les conflits manuellement ou exécutez :"
  echo "   ./import.sh -f  pour écraser les modifications locales."
  exit 1
}

echo "Import terminé avec succès."
