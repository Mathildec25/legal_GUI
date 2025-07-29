# Safelab – Interface for the Legal Analysis of Regulated Chemical Substances

 **Safelab** is an interactive interface designed to display the legal obligations associated with regulated chemical substances **(narcotics, psychotropics, and precursors)**.
 
The user can:
- search by name, CAS number, or SMILES,

- draw a custom molecular structure,

- automatically search for a substance using the SMILES generated from the drawing.

## Installation

**Clone the repository and navigate to the project folder :**
```text
git clone https://github.com/Mathildec25/legal_GUI.git

cd legal_GUI
```
**Create the Conda environment using the environment.yml file :**
```text
conda env create -f environment.yml
conda activate rdkit-env
```
The environment uses conda-forge to install RDKit and other dependencies.

**Kekule.js : Molecular Editor**

The HTML file **assets/kekule_editor.html** includes a local version of the Kekule.js source code

(see: https://github.com/partridgejiang/Kekule.js).

No additional npm installation is required.

Optional manual installation via **npm**:
```text
npm install kekule
```
**Run the Dash application :**
```text
python app.py
```

## Key Features

- Search for substances by name, CAS number, or SMILES

- Draw molecular structures using an integrated editor (Kekule.js), with automatic SMILES generation

- Display legal information associated with each substance:

   - Classification under Belgian annexes (narcotics, psychotropics)

   - Precursor category (EU Regulations 273/2004 and 111/2005)

   - Dynamic legal messages with links to official sources

- Visualize the molecular structure from SMILES using RDKit

- Responsive interface with clear navigation between the search page and the database

- Interactive filtering and alphabetical sorting (A to Z or Z to A) in the database

- Export the database as a CSV file for external analysis


## User Guide
 ### Search Page
 When the application opens, the search page appears by default. You can perform several actions:

 **1. Search by name, CAS number, or SMILES**

- Use the central input field (dropdown) to type:

  - the molecule’s name (e.g., cannabidiol),

  - its CAS number (e.g., 13956-29-1),

  - or its SMILES representation.

- Click the Search button to display the legal and chemical information related to the selected substance.

**2. Draw a molecule manually (Draw)**

- Click the Draw button to open a modal window.

- On the left side, the integrated Kekule.js editor allows you to draw your molecule manually.

- Once finished, click Get SMILES: the generated SMILES will automatically be inserted into the bottom input field.

- Then click Search to display the legal and database-related information at the bottom of the window (just scroll down).

 You can also paste a canonical SMILES directly into the input bar, click Draw to view the structure on the right, and then click Search.

 **Warning: Always use a canonical SMILES to ensure accurate detection.**

- To close the drawing window, click the Close button.

### Page Navigation
  A sidebar on the left allows navigation between pages:

 - Search : default page (search input, drawing tool, results)

 - Database : view the complete database

### Database Page
  The Full Database page allows you to:

  - view all substances listed in donnees.csv,

  - sort columns alphabetically (A → Z),

  - dynamically filter results,

  - export the database as a CSV file via the Export button.
  ## Data Used

  The application relies on a structured CSV database '(data/donnees.csv)' with the following columns:
   - **`name`** : substance name
   - **`cas`** : CAS number
   - **`smiles` / `smiles_canonique`** : chemical representation
   - **`belgian_annex`** : relevant Belgian annex (I, II, III, IV)
   - **`EU Annex/Category 273/2004` / `EU Annex/Category 111/2005`** : precursor category under EU regulations
  
The database is not yet exhaustive: some narcotics, psychotropics, and precursors may be missing.
The database is automatically loaded at startup by the utils/data_loader.py script.

## Project Structure

```text
legal_GUI/
├── app.py                      # Main entry point – initializes the Dash app and registers callbacks
├── environment.yml             # Conda environment file
├── assets/                     # Static assets: CSS, icons, JS, Kekule
│   ├── kekule.js_master/       # Kekule.js source code (molecular editor)
│   ├── kekule_editor.html      # Custom integration of the Kekule editor
│   ├── clientside.js           # Client-side JS functions for Dash (e.g. retrieve SMILES from Kekule)
│   ├── style.css               # Custom CSS styles
│   ├── fontawesome.css         # Font Awesome icons
│   └── favicon.ico             # App icon
├── callbacks/                  # Dash callbacks (user interaction logic)
│   └── substance_callbacks.py  # Callbacks for search, drawing, and result display
├── data/
│   └── donnees.csv             # Substances database
├── layout/                     # Interface layout components
│   ├── main_layout.py          # Overall structure (navigation, sidebar, content)
│   ├── search_page.py          # Search page layout
│   └── database_page.py        # Interactive database display
├── utils/                      # Utility functions
│   ├── data_loader.py          # Loads and cleans the database (normalizes column names)
│   ├── draw_rdkit.py           # Molecule rendering from SMILES using RDKit
│   ├── Narc_psy_procedure.py   # Legal info for narcotics / psychotropics
│   └── Prec_procedure.py       # Legal info for precursors

```


## Technologies Used
 **VS Code**

**Python 3.9**

**Dash and Dash Bootstrap Components**

**RDKit – molecular visualization**

**Kekule.js – HTML/JS molecular editor**

**Pandas,HTML/CSS, Font Awesome**

## Notes
In the **Draw** section, it is crucial to use a **canonical SMILES** to ensure proper detection of substances and accurate display of legal information.

Detection and display in the interface are highly dependent on:

- **strictly matching column names**, including **spaces and capitalization**, in the CSV file
(e.g., `Name`, `belgian_annex`, `EU Annex/Category 273/2004`, `EU Annex/Category 111/2005`, etc.),

- **how the values are entered in these columns**(e.g., `Annexe I`, `1`, `I` may yield different results)

Currently, Category 4 precursors and Annex IV of the Belgian Royal Decree are not yet handled.

At load time, column names are automatically converted to lowercase, spaces are replaced with underscores (`_`), and extra spaces are stripped.

Therefore, it is essential to **strictly follow the structure** of the donnees.csv file to ensure the application functions correctly.

If the database is modified (e.g., columns added, names changed, new formats introduced),
**the Python code must be updated accordingly**, especially the detection conditions in the legal analysis functions.



