from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def draw_substructure_svg(smiles, smarts, width=500, height=500):
    """
    Génère un dessin SVG d'une molécule à partir d’un SMILES, 
    en mettant en surbrillance une sous-structure définie par un motif SMARTS.

    Paramètres :
    - smiles (str) : chaîne SMILES représentant la molécule complète.
    - smarts (str) : motif SMARTS à rechercher dans la molécule.
    - width (int) : largeur de l’image SVG.
    - height (int) : hauteur de l’image SVG.

    Retour :
    - svg (str) : code SVG avec la sous-structure mise en évidence.
    """
    # Crée la molécule et le motif depuis les chaînes SMILES/SMARTS
    mol = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts(smarts)

    # Vérifie si la molécule ou le motif est invalide
    if mol is None or patt is None:
        return "<svg><!-- Invalid molecule or pattern --></svg>"

    # Recherche des atomes correspondant au motif SMARTS
    match = mol.GetSubstructMatch(patt)
    if not match:
        return "<svg><!-- No match found --></svg>"

    hit_ats = list(match)  # Liste des indices des atomes correspondants
    hit_bonds = []         # Liste des indices des liaisons correspondantes

    # Identifie les liaisons du motif à mettre en surbrillance
    for bond in patt.GetBonds():
        b_idx = bond.GetBeginAtomIdx()  # index relatif au motif
        e_idx = bond.GetEndAtomIdx()

        # Vérifie que les indices sont valides dans la correspondance
        if b_idx >= len(hit_ats) or e_idx >= len(hit_ats):
            continue  # évite un IndexError

        # Récupère les indices des atomes dans la molécule d'origine
        aid1 = hit_ats[b_idx]
        aid2 = hit_ats[e_idx]

        # Récupère l'objet liaison entre les deux atomes
        bond_obj = mol.GetBondBetweenAtoms(aid1, aid2)
        if bond_obj:  # Ajoute si la liaison existe bien
            hit_bonds.append(bond_obj.GetIdx())

    # Création du dessinateur SVG avec dimensions personnalisées
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

    # Définition des couleurs de surbrillance (gris clair ici)
    highlight_colors = {idx: (0.9, 0.9, 0.9) for idx in hit_ats}
    bond_colors = {idx: (0.9, 0.9, 0.9) for idx in hit_bonds}

    # Dessin de la molécule avec les atomes/liaisons en surbrillance
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, mol,
        highlightAtoms=hit_ats,
        highlightBonds=hit_bonds,
        highlightAtomColors=highlight_colors,
        highlightBondColors=bond_colors
    )

    # Finalisation du dessin SVG
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Nettoyage du SVG : suppression de l'en-tête XML inutile
    return svg.replace('<?xml version="1.0" encoding="UTF-8"?>', '')
