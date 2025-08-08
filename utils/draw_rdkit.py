# utils/draw_rdkit.py
from typing import Optional, Tuple, List

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDraw2D


def _find_substructure(mol: Chem.Mol, smarts: Optional[str]) -> Tuple[List[int], List[int]]:
    """Retourne (atomIndices, bondIndices) de la première correspondance SMARTS."""
    if not smarts:
        return [], []
    patt = Chem.MolFromSmarts(smarts)
    if not patt:
        return [], []
    match = mol.GetSubstructMatch(patt)
    if not match:
        return [], []
    hit_atoms = list(match)
    hit_bonds = []
    hit_set = set(hit_atoms)
    for b in mol.GetBonds():
        if b.GetBeginAtomIdx() in hit_set and b.GetEndAtomIdx() in hit_set:
            hit_bonds.append(b.GetIdx())
    return hit_atoms, hit_bonds


def draw_substructure_svg(
    smiles: str,
    smarts: Optional[str] = None,
    size: Tuple[int, int] = (420, 260),
) -> str:
    """
    Dessine un SVG RDKit à partir d’un SMILES.
    - smarts: sous-structure à surligner (optionnelle).
    - size: (width, height).
    Retourne le SVG (string).
    """
    if not smiles:
        return "<svg><!-- empty --></svg>"

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "<svg><!-- invalid smiles --></svg>"

    # Calcul coords 2D (idempotent, rapide si déjà présentes)
    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    # Options “lightweight” pour aller vite
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = False
    opts.explicitMethyl = False
    opts.addAtomIndices = False
    opts.addBondIndices = False
    opts.fixedBondLength = 25.0  # un peu plus compact

    # Surlignage sous-structure si demandé
    hit_atoms, hit_bonds = _find_substructure(mol, smarts)
    if hit_atoms or hit_bonds:
        highlight_atom_map = {idx: (0.8, 0.0, 0.0) for idx in hit_atoms}  # couleur par défaut (neutre RDKit)
        highlight_bond_map = {idx: (0.8, 0.0, 0.0) for idx in hit_bonds}
        drawer.DrawMolecule(
            mol,
            highlightAtoms=hit_atoms,
            highlightBonds=hit_bonds,
            highlightAtomColors=highlight_atom_map,
            highlightBondColors=highlight_bond_map,
        )
    else:
        drawer.DrawMolecule(mol)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # Petit nettoyage
    svg = svg.replace("svg:","")
    return svg
