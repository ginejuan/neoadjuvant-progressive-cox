# Cómo subir este repositorio a GitHub

Esta carpeta ya está lista. Ejecuta lo siguiente **en tu ordenador** (necesitas
Git instalado y una cuenta de GitHub).

## Opción A — con GitHub CLI (`gh`), la más rápida
```bash
cd "ruta/a/github_repo_neoadjuvant-progressive-cox"
git init
git add .
git commit -m "Initial commit: reproducible analysis and de-identified data"
gh repo create neoadjuvant-progressive-cox --public --source=. --remote=origin --push
```

## Opción B — sin gh (crear el repo en la web)
1. En https://github.com/new crea un repositorio llamado
   `neoadjuvant-progressive-cox` (público), **sin** README/LICENSE (ya los tienes).
2. En una terminal:
```bash
cd "ruta/a/github_repo_neoadjuvant-progressive-cox"
git init
git add .
git commit -m "Initial commit: reproducible analysis and de-identified data"
git branch -M main
git remote add origin https://github.com/USUARIO/neoadjuvant-progressive-cox.git
git push -u origin main
```
(Sustituye `USUARIO` por tu usuario de GitHub.)

## Recomendado: DOI citable con Zenodo
Conecta el repo en https://zenodo.org (login con GitHub) → activa el repositorio →
crea un *release* en GitHub. Zenodo generará un DOI que puedes citar en el
apartado *Code availability* del manuscrito.
