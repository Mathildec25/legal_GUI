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

- Recherche de substances par nom, numéro CAS ou SMILES

- Dessin de structures moléculaires via un éditeur intégré (Kekule.js), avec génération automatique du SMILES

- Affichage des informations juridiques liées à chaque substance :

   - Classification dans les annexes belges (stupéfiants, psychotropes)

   - Catégorie de précurseur (règlements UE 273/2004 et 111/2005)

   - Messages juridiques dynamiques avec liens vers les sources officielles

- Visualisation de la structure moléculaire à partir du SMILES avec RDKit

- Interface responsive avec navigation claire entre la page de recherche et la base de données

- Filtrage interactif et tri alphabétique (A à Z ou Z à A) dans la base de données.

- Téléchargement de la base au format CSV pour une analyse externe

## Guide d’utilisation
 ### Page de recherche (Search)
Dès l’ouverture de l’application, vous arrivez sur la page de recherche, qui permet plusieurs actions :

**1. Recherche par nom, numéro CAS ou SMILES**
   - Utilisez la zone de saisie centrale (champ déroulant) pour taper :

     - le nom de la molécule (ex. : cannabidiol),

     - son numéro CAS (ex. : 13956-29-1),

     - ou sa représentation SMILES.

  - Cliquez sur le bouton  Search pour afficher les informations juridiques et chimiques associées à la substance recherchée.

**2. Dessiner une molécule manuellement (Draw)**
   - Cliquez sur le bouton Draw pour ouvrir une fenêtre modale (popup).

   - Sur la gauche, un éditeur intégré (Kekule.js) vous permet de dessiner votre molécule à la main.

   - Une fois le dessin terminé, cliquez sur Get SMILES : le SMILES généré sera automatiquement inséré dans la barre en bas.

   - Cliquez ensuite sur Search pour faire apparaître les informations juridiques et celles issues de la base, en bas de la fenêtre (il suffit de scroller).


   Vous pouvez également coller directement un SMILES canonique dans la barre, cliquer sur Draw pour voir le dessin à droite, puis cliquer sur Search.

   **Attention : toujours utiliser un SMILES canonique pour garantir l’exactitude de la détection.**

   - Pour fermer la fenêtre de dessin, cliquez sur le bouton Close.

  ### Navigation entre les pages
Sur la gauche, un menu latéral (sidebar) permet de passer d'une page à l'autre :

  - Search : page par défaut (zone de recherche, éditeur, résultats)

  - Database : visualisation de la base de données complète


  ### Page "Database"
  La page Full Database permet de :

   - voir toutes les substances listées dans la base donnees.csv,
   - trier les colonnes (de A → Z),
   - filtrer dynamiquement les résultats,
   - exporter la base au format CSV via le bouton Export.
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
├── app.py # point d'entrée principal, initialise l’application Dash et enregistre les callbacks.
├── environnement.yml # Fichier pour créer l'environnement Conda
├── assets/ # Fichiers statiques : CSS, icônes, JS, Kekule
│   ├── kekule.js_master/ # Code source de Kekule.js (éditeur moléculaire)
│   ├── kekule_editor.html # Intégration personnalisée de l'éditeur Kekule
│   ├── clientside.js # Fonctions JavaScript côté client pour Dash (récupération du SMILES dessiné dans l’éditeur Kekule)
│   ├── style.css # Style personnalisé de l'interface
│   ├── fontawesome.css # Icônes Font Awesome
│   └── favicon.ico # Icône de l'application
├── callbacks/ # Fichiers contenant les callbacks Dash(Fichiers de gestion des interactions utilisateur)
│   └── substance_callbacks.py # Callbacks Dash pour la recherche, le dessin, l’affichage des résultats

├── data/
│   └── donnees.csv # Base de données des substances
├── layout/ # Composants de mise en page de l'interface
│   ├── main_layout.py # Structure générale (navigation, sidebar, contenu)
│   ├── search_page.py # Page de recherche de substances
│   └── databse_page.py #  Page d’affichage interactif de la base de données
├── utils/ # Fonctions utilitaires
│   ├── data_loader.py # Chargement et nettoyage de la base (normalisation des noms de colonnes)
│   ├── draw_rdkit.py # Dessin de la molécule à partir d’un SMILES (RDKit)
│   ├── Narc_psy-procedure.py # Affichage dynamique des obligations pour les stupéfiants / psychotropes
│   └── Prec_procedure.py # Affichage dynamique des obligations pour les précurseurs
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

