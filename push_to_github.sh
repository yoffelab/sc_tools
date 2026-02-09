#!/bin/bash
# One-time script to push sc_tools to GitHub. Run from project root:
#   chmod +x push_to_github.sh && ./push_to_github.sh
set -e
cd "$(dirname "$0")"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  git init
  git config core.hooksPath .githooks
fi

git remote remove origin 2>/dev/null || true
git remote add origin https://github.com/yoffelab/sc_tools.git

git add -A
git status

echo ""
echo "If any file over 1MB is staged, the pre-commit hook will block the commit."
echo "Unstage with: git reset HEAD -- <file>"
echo ""

git commit -m "Initial commit: sc_tools scientific analysis tool for spatial and single-cell multiomics"
git branch -M main
git push -u origin main

echo "Done."
