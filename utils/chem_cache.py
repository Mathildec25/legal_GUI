# utils/chem_cache.py
from __future__ import annotations

from functools import lru_cache
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import rdDepictor

# Rendu SVG (léger) via ton utilitaire existant
from utils.draw_rdkit import draw_substructure_svg


# ---------- Helpers internes (cachés) ----------

@lru_cache(maxsize=50000)
def _canonicalize_smiles_cached(raw_smiles: str) -> str:
    """Canonise un SMILES (via RDKit). Retourne '' si invalide."""
    if not raw_smiles:
        return ""
    mol = Chem.MolFromSmiles(raw_smiles)
    if not mol:
        return ""
    return Chem.MolToSmiles(mol, canonical=True)


@lru_cache(maxsize=20000)
def _mol_from_canonical_smiles(canonical_smiles: str):
    """Construit un Mol RDKit à partir d’un SMILES *canonique* (cache)."""
    if not canonical_smiles:
        return None
    mol = Chem.MolFromSmiles(canonical_smiles)
    if not mol:
        return None
    # Coordonnées 2D calculées une fois et mémorisées avec l'objet
    rdDepictor.Compute2DCoords(mol)
    return mol


def _minify_svg(svg: str) -> str:
    """Nettoie un peu le SVG retourné (facultatif)."""
    if not svg:
        return svg
    # RDKit insère parfois 'svg:' dans les balises
    return svg.replace("svg:", "").strip()


@lru_cache(maxsize=20000)
def _svg_from_canonical_smiles(
    canonical_smiles: str,
    width: int,
    height: int,
) -> str:
    """Rendu SVG (sans SMARTS) à partir d’un SMILES *canonique* (cache)."""
    if not canonical_smiles:
        return "<svg><!-- empty --></svg>"
    svg = draw_substructure_svg(canonical_smiles, smarts=None, size=(width, height))
    return _minify_svg(svg)


@lru_cache(maxsize=20000)
def _svg_from_canonical_with_smarts(
    canonical_smiles: str,
    smarts: Optional[str],
    width: int,
    height: int,
) -> str:
    """Rendu SVG (avec SMARTS) à partir d’un SMILES *canonique* (cache)."""
    if not canonical_smiles:
        return "<svg><!-- empty --></svg>"
    svg = draw_substructure_svg(canonical_smiles, smarts=smarts, size=(width, height))
    return _minify_svg(svg)


# ---------- API publique (à utiliser dans tes callbacks) ----------

def canonicalize_smiles(raw_smiles: str) -> str:
    """
    Canonise un SMILES en passant par un cache LRU.
    Retourne '' si invalide.
    """
    return _canonicalize_smiles_cached(raw_smiles or "")


def mol_from_smiles(raw_smiles: str):
    """
    Retourne un Mol RDKit à partir d’un SMILES quelconque.
    - Canonise d’abord
    - Retourne None si invalide
    """
    can = canonicalize_smiles(raw_smiles)
    if not can:
        return None
    return _mol_from_canonical_smiles(can)


def svg_from_smiles(raw_smiles: str, size: Tuple[int, int] = (420, 260)) -> str:
    """
    Retourne un SVG (string) à partir d’un SMILES quelconque (pas de SMARTS).
    - Canonise d’abord
    - Met en cache par SMILES canonique + dimensions
    """
    can = canonicalize_smiles(raw_smiles)
    if not can:
        return "<svg><!-- empty --></svg>"
    w, h = size
    return _svg_from_canonical_smiles(can, int(w), int(h))


def svg_from_smiles_with_smarts(
    raw_smiles: str,
    smarts: Optional[str],
    size: Tuple[int, int] = (420, 260),
) -> str:
    """
    Retourne un SVG (string) avec surlignage d’une sous-structure SMARTS.
    - Canonise d’abord
    - Met en cache par (SMILES canonique, SMARTS, dimensions)
    """
    can = canonicalize_smiles(raw_smiles)
    if not can:
        return "<svg><!-- empty --></svg>"
    w, h = size
    return _svg_from_canonical_with_smarts(can, smarts, int(w), int(h))


def clear_all_caches() -> None:
    """Vide tous les caches (utile en hot‑reload)."""
    _canonicalize_smiles_cached.cache_clear()
    _mol_from_canonical_smiles.cache_clear()
    _svg_from_canonical_smiles.cache_clear()
    _svg_from_canonical_with_smarts.cache_clear()
