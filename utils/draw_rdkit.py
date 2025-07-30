from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def draw_substructure_svg(smiles, smarts, width=500, height=500):
    """
    Generates an SVG drawing of a molecule from a SMILES string,
    highlighting a substructure defined by a SMARTS pattern.

    Parameters:
    - smiles (str): SMILES string representing the full molecule.
    - smarts (str): SMARTS pattern to search for in the molecule.
    - width (int): width of the SVG image.
    - height (int): height of the SVG image.

    Returns:
    - svg (str): SVG code with the substructure highlighted.
    """
    # Create the molecule and pattern from SMILES/SMARTS strings
    mol = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts(smarts)

    # Check if the molecule or pattern is invalid
    if mol is None or patt is None:
        return "<svg><!-- Invalid molecule or pattern --></svg>"

    # Find atoms matching the SMARTS pattern
    match = mol.GetSubstructMatch(patt)
    if not match:
        return "<svg><!-- No match found --></svg>"

    hit_ats = list(match)  # List of indices of matching atoms
    hit_bonds = []         # List of indices of matching bonds

    # Identify bonds in the pattern to highlight
    for bond in patt.GetBonds():
        b_idx = bond.GetBeginAtomIdx()  # index relative to pattern
        e_idx = bond.GetEndAtomIdx()

        # Make sure indices are valid in the match
        if b_idx >= len(hit_ats) or e_idx >= len(hit_ats):
            continue  # avoid IndexError

        # Get the atom indices in the original molecule
        aid1 = hit_ats[b_idx]
        aid2 = hit_ats[e_idx]

        # Get the bond object between the two atoms
        bond_obj = mol.GetBondBetweenAtoms(aid1, aid2)
        if bond_obj:  # Add if the bond exists
            hit_bonds.append(bond_obj.GetIdx())

    # Create the SVG drawer with custom dimensions
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

    # Define highlight colors (light grey here)
    highlight_colors = {idx: (0.9, 0.9, 0.9) for idx in hit_ats}
    bond_colors = {idx: (0.9, 0.9, 0.9) for idx in hit_bonds}

    # Draw the molecule with highlighted atoms/bonds
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, mol,
        highlightAtoms=hit_ats,
        highlightBonds=hit_bonds,
        highlightAtomColors=highlight_colors,
        highlightBondColors=bond_colors
    )

    # Finalize the SVG drawing
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Clean the SVG: remove unnecessary XML header
    return svg.replace('<?xml version="1.0" encoding="UTF-8"?>', '')
