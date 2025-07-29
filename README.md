# Safelab_ Interface pour l’analyse juridique des molécules soumises à autorisation spécifique

 **Safelab** est une interface interactive qui permet d’afficher les obligations juridiques associées à des substances chimiques réglementées **(stupéfiants, psychotropes, précurseurs)**.
 L’utilisateur peut :
  - faire une recherche par **nom, numéro CAS ou SMILES**,
  - dessiner une structure moléculaire personnalisée,
  - rechercher automatiquement une substance à partir du SMILES généré par le dessin.

## Installation

**Cloner le dépôt et se placer dans le dossier du projet :**
```text
git clone https://github.com/ton-utilisateur/legal_GUI.git

cd legal_GUI
```
**Créer l’environnement Conda à partir du fichier environment.yml :**
```text
conda env create -f environment.yml
conda activate rdkit-env
L’environnement utilise conda-forge pour installer RDKit et les autres dépendances.
```
**Kekule.js : éditeur de molécules**

Le fichier HTML **assets/kekule_editor.html** contient une version locale du code source de Kekule.js (déjà intégré).
Aucune installation supplémentaire via npm n’est requise.
Installation aussi possible manuellement via npm :
```text
npm install kekule
```
**Lancer l'application Dash :**
```text
python app.py
```

## Fonctionnalités principales

° Recherche de substances par nom, numéro CAS ou SMILES
° Dessin de structures moléculaires via un éditeur intégré (Kekule.js), avec génération automatique du SMILES
° Affichage des informations juridiques liées à chaque substance :

   - Classification dans les annexes belges (stupéfiants, psychotropes)

   - Catégorie de précurseur (règlements UE 273/2004 et 111/2005)

   - Messages juridiques dynamiques avec liens vers les sources officielles

° Visualisation de la structure moléculaire à partir du SMILES avec RDKit

° Interface responsive avec navigation claire entre la page de recherche et la base de données
° Filtrage interactif et tri alphabétique (A à Z ou Z à A) dans la base de données.
° Téléchargement de la base au format CSV pour une analyse externe
## Données utilisées

L’application repose sur une **base de données CSV** (`data/donnees.csv`) structurée avec les colonnes suivantes :

- **`name`** : nom de la substance
- **`cas`** : numéro CAS
- **`smiles` / `smiles_canonique`** : représentation chimique
- **`belgian_annex`** : annexe belge applicable (I, II, III, IV)
- **`EU Annex/Category 273/2004` / `EU Annex/Category 111/2005`** : catégorie de précurseur selon les règlements européens.
  
La base de données utilisée n’est pas encore exhaustive : certains stupéfiants, psychotropes et précurseurs peuvent manquer
Cette base est automatiquement chargée au démarrage par le fichier `utils/data_loader.py`.

## Structure du projet

```text
legal_GUI/
├── app.py # Point d'entrée principal de l'application Dash
├── environnement.yml # Fichier pour créer l'environnement Conda
├── assets/ # Fichiers statiques : CSS, icônes, JS, Kekule
│   ├── kekule.js_master/ # Code source de Kekule.js (éditeur moléculaire)
│   ├── kekule_editor.html # Intégration personnalisée de l'éditeur Kekule
│   ├── clientside.js # Fonctions JavaScript côté client pour Dash ,Permet de récupérer le SMILES dessiné dans l’éditeur Kekule
│   ├── style.css # Style personnalisé de l'interface
│   ├── fontawesome.css # Icônes Font Awesome
│   └── favicon.ico # Icône de l'application
├── callbacks/ # Fichiers contenant les callbacks Dash
│   └── substance_callbacks.py # Callbacks pour la recherche et les résultats
├── data/
│   └── donnees.csv # Base de données des substances
├── layout/ # Composants de mise en page de l'interface
│   ├── main_layout.py # Agencement principal de l'application
│   ├── search_page.py # Page de recherche de substances
│   └── databse_page.py # Page d'affichage de la base de données
├── utils/ # Fonctions utilitaires
│   ├── data_loader.py # Chargement de la base de données
│   ├── draw_rdkit.py # Visualisation moléculaire avec RDKit
│   ├── Narc_psy-procedure.py # Logique juridique : stupéfiants/psychotropes
│   └── Prec_procedure.py # Logique juridique : précurseurs
```


## Technologies utilisées
 **VS Code**
 
  **Python 3.9**

  **Dash et Dash Bootstrap Components**

  **RDKit – visualisation moléculaire**

  **Kekule.js – éditeur moléculaire HTML/JS**

  **Pandas, Plotly, HTML/CSS, Font Awesome**

## Remarque
Dans la section**Draw**, il est impératif d’utiliser un **SMILES canonique** pour garantir une détection correcte des substances et un affichage fiable des résultats juridiques.

Les conditions de détection et d’affichage dans l’interface dépendent fortement :

- du **respect exact des noms de colonnes**, y compris les **espaces et majuscules**, dans le fichier CSV  
  (exemples : `Name`, `belgian_annex`, `EU Annex/Category 273/2004`, `EU Annex/Category 111/2005`, etc.),
- de la **façon dont les valeurs sont saisies** dans ces colonnes (ex. : `Annexe I`, `1`, `I` peuvent produire des résultats différents).
  
De plus,**la détection automatique ne prend pas encore en compte la catégorie 4 des précurseurs ni l’annexe IV de l’arrêté royal belge.**

Lors du chargement, les noms de colonnes sont automatiquement transformés en minuscules, les espaces sont remplacés par des tirets bas (`_`), et les espaces superflus sont supprimés.

Il est donc essentiel de **respecter strictement la structure** du fichier `donnees.csv` pour garantir le bon fonctionnement de l’application.

 En cas de modification de la base de données (ajout de colonnes, changement de nom, nouveaux formats),  
**le code Python devra être ajusté en conséquence**, notamment les conditions de détection dans les fonctions d’analyse juridique.

