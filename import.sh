#!/bin/bash
# Script d'importation (synchronisation depuis GitHub)
# Usage :
#   ./import.sh       → mise à jour normale avec rebase (préserve tes modifs locales)
#   ./import.sh -f    → forçage total (remplace TOUT par la version du dépôt distant)

# Couleurs pour la lisibilité
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}🔍 Vérification du dépôt Git...${NC}"
if [ ! -d .git ]; then
  echo -e "${RED}❌ Ce répertoire n'est pas un dépôt Git !${NC}"
  exit 1
fi

# Si -f est passé → reset complet
if [ "$1" == "-f" ]; then
  echo -e "${YELLOW}⚠️  Import forcé : toutes les modifications locales seront perdues${NC}"
  git fetch origin main || { echo -e "${RED}❌ Échec du fetch${NC}"; exit 1; }
  git reset --hard origin/main || { echo -e "${RED}❌ Échec du reset${NC}"; exit 1; }
  git clean -fdx || { echo -e "${RED}❌ Échec du nettoyage${NC}"; exit 1; }
  echo -e "${GREEN}✅ Répertoire local entièrement synchronisé avec le dépôt distant.${NC}"

else
  # Importation normale
  echo -e "${YELLOW}📦 Mise à jour du dépôt (rebase)...${NC}"
  git pull --rebase origin main || {
    echo -e "${RED}⚠️  Conflit détecté !${NC}"
    echo "➡️  Corrigez les conflits manuellement, ou relancez avec :"
    echo "   ./import.sh -f   pour forcer la version distante."
    exit 1
  }
  echo -e "${GREEN}✅ Mise à jour terminée sans conflit.${NC}"
fi

echo "Import terminé avec succès."
