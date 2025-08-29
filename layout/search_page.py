from dash import html, dcc
import dash_bootstrap_components as dbc
from utils.data_loader import df
import pandas as pd

def _pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

# Définis des alias possibles pour chaque champ
FIELD_CANDIDATES = {
    "NAME": ["NAME", "Name", "Substance", "Substance Name", "nom"],
    "IUPAC": ["IUPAC", "IUPAC Name", "Nom IUPAC"],
    "CAS": ["CAS", "CAS Number", "CAS_No", "CAS number"],
    "SMILES": ["SMILES", "Smiles"],
    "canonical_smiles": ["canonical_smiles", "Canonical_Smiles", "Canonical SMILES"],
}

options = []
for target, candidates in FIELD_CANDIDATES.items():
    col = _pick_col(df, candidates)
    if col is None:
        continue
    for val in df[col]:
        if pd.notna(val):
            options.append({"label": val, "value": val})

def get_search_page(hidden=False):
    """
    Page de recherche repensée avec design moderne
    """
    display = "none" if hidden else "block"
    
    return html.Div([
        
        # Section de recherche principale avec card moderne
        html.Section([
            dbc.Container([
                dbc.Card([
                    dbc.CardBody([
                        # Titre de section
                        html.H2([
                            html.I(className="fas fa-search me-2 text-primary"),
                            "Molecular Search"
                        ], className="search-section-title"),
                        
                        html.P(
                            "Choose your preferred search method and find detailed regulatory information",
                            className="search-section-subtitle"
                        ),
                        
                        # Méthodes de recherche
                        dbc.Row([
                            dbc.Col([
                                dbc.Card([
                                    dbc.CardBody([
                                        html.I(className="fas fa-search search-method-icon"),
                                        html.H5("Text Search", className="mt-2"),
                                        html.P("Search by name, CAS number, or SMILES notation", 
                                              className="text-muted small")
                                    ], className="text-center search-method-card active")
                                ], className="search-method", id="method-text")
                            ], md=4),
                            
                            dbc.Col([
                                dbc.Card([
                                    dbc.CardBody([
                                        html.I(className="fas fa-pencil-alt search-method-icon"),
                                        html.H5("Draw Molecule", className="mt-2"),
                                        html.P("Draw your molecule structure interactively", 
                                              className="text-muted small")
                                    ], className="text-center search-method-card")
                                ], className="search-method", id="method-draw")
                            ], md=4),
                            
                            dbc.Col([
                                dbc.Card([
                                    dbc.CardBody([
                                        html.I(className="fas fa-upload search-method-icon"),
                                        html.H5("Upload File", className="mt-2"),
                                        html.P("Import from MOL, SDF or other formats", 
                                              className="text-muted small")
                                    ], className="text-center search-method-card")
                                ], className="search-method", id="method-upload")
                            ], md=4)
                        ], className="mb-4 search-methods-row"),
                        
                        # Barre de recherche moderne
                        html.Div([
                            dbc.InputGroup([
                                dcc.Dropdown(
                                    id='input-substance',
                                    options=options,
                                    placeholder="Enter molecule name, CAS number, or SMILES...",
                                    searchable=True,
                                    className="modern-search-dropdown"
                                ),
                                dbc.Button([
                                    html.I(className="fas fa-search me-2"),
                                    html.Span("Search", className="d-none d-md-inline")
                                ],
                                id="search-button",
                                color="primary",
                                className="search-btn-modern")
                            ], className="search-input-group")
                        ], className="search-input-container")
                        
                    ])
                ], className="search-main-card shadow-lg")
            ], className="search-container")
        ], className="search-section-wrapper"),
        
        # Résultats de recherche
        html.Section([
            dbc.Container([
                html.Div(id="substance-details", className="results-container")
            ])
        ], className="results-section"),
        
        # Modal modernisé pour le dessin
        dbc.Modal([
            dbc.ModalHeader([
                html.H4([
                    html.I(className="fas fa-pencil-alt me-2"),
                    "Molecular Editor"
                ], className="modal-title-modern")
            ]),
            dbc.ModalBody([
                # Instructions utilisateur
                dbc.Alert([
                    html.I(className="fas fa-info-circle me-2"),
                    "Draw your molecule using the interactive editor below, then click 'Get SMILES' to extract the structure."
                ], color="info", className="mb-3"),
                
                # Interface de dessin
                dbc.Row([
                    dbc.Col([
                        html.H6([
                            html.I(className="fas fa-pen me-2"),
                            "Draw Molecule"
                        ], className="mb-3"),
                        html.Iframe(
                            src="/assets/kekule_editor.html",
                            className="kekule-editor-modern",
                            id="kekule-editor"
                        )
                    ], md=6),
                    
                    dbc.Col([
                        html.H6([
                            html.I(className="fas fa-eye me-2"),
                            "Preview & Visualization"
                        ], className="mb-3"),
                        html.Div(
                            id="display-rdkit-svg",
                            className="rdkit-display-modern"
                        )
                    ], md=6)
                ], className="drawing-interface-row"),
                
                # Input SMILES
                html.Div([
                    html.Label("Generated SMILES:", className="form-label fw-bold"),
                    dbc.InputGroup([
                        dbc.Input(
                            id="manual-smiles",
                            placeholder="SMILES will appear here after drawing...",
                            className="smiles-input-modern"
                        ),
                        dbc.Button([
                            html.I(className="fas fa-copy me-1"),
                            "Copy"
                        ], color="outline-secondary", id="copy-smiles")
                    ])
                ], className="mt-3 smiles-input-section"),
                
                # Boutons d'action
                html.Div([
                    dbc.ButtonGroup([
                        dbc.Button([
                            html.I(className="fas fa-flask me-2"),
                            "Get SMILES"
                        ],
                        id="get-smiles-btn",
                        color="info",
                        className="action-btn"),
                        
                        dbc.Button([
                            html.I(className="fas fa-play me-2"),
                            "Preview"
                        ],
                        id="load-smiles",
                        color="secondary",
                        className="action-btn"),
                        
                        dbc.Button([
                            html.I(className="fas fa-search me-2"),
                            "Search Molecule"
                        ],
                        id="modal-search-button",
                        color="success",
                        className="action-btn")
                    ], className="d-flex gap-2 justify-content-center")
                ], className="mt-4 action-buttons-section"),
                
                # Store pour les données
                dcc.Store(id="raw-smiles"),
                
                # Résultats du modal
                html.Div(id="modal-search-results", className="mt-4")
                
            ]),
            dbc.ModalFooter([
                dbc.Button([
                    html.I(className="fas fa-times me-2"),
                    "Close"
                ],
                id="close-draw",
                color="light",
                className="modal-close-btn")
            ])
        ],
        id="modal-draw",
        is_open=False,
        size="xl",
        className="modern-modal",
        backdrop="static")
        
    ], 
    style={"display": display},
    className="search-page-container")