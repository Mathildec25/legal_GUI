# search_page.py

from dash import html, dcc
import dash_bootstrap_components as dbc
from utils.data_loader import df
import pandas as pd

# Options du Dropdown : tous les identifiants disponibles dans la base
options = []
for col in ['name', 'cas', 'smiles', 'smiles_canonique']:
    for val in df[col]:
        if pd.notna(val):
            options.append({"label": val, "value": val})

#  Fonction qui retourne le layout de la page de recherche
def get_search_page(hidden=False):
    display = "none" if hidden else "block"

    return html.Div([

        # Message d’introduction centré
        html.Div(
            html.H5(
                "Enter a chemical name or Cas_NO or Smiles to begin",
                style={"color": "#003366", "marginBottom": "20px"}
            ),
            style={"marginLeft": "190px"}
        ),

        # Barre de recherche et bouton "Draw"
        dbc.Row([
            dbc.Col(
                dbc.Button([
                    html.I(className="fas fa-pencil-alt", style={"marginRight": "6px"}),  # ✏️ Font Awesome
                    "Draw"
                ],
                id="draw-button",
                n_clicks=0,
                style={"background-color": "#066ab1", "color": "white"}),
                width="auto"
            ),

            dbc.Col(
                dcc.Dropdown(
                    id='input-substance',
                    options=options,
                    placeholder="Search by name, CAS number, or SMILES...",
                    searchable=True
                ),
                width=8
            ),

            dbc.Col(
                dbc.Button([
                    html.I(className="fas fa-search", style={"marginRight": "6px"}),  # Font Awesome
                    "Search"
                ],
                id="search-button",
                n_clicks=0,
                style={"background-color": "#066ab1", "color": "white"}),
                width="auto"
            )
        ],
        className="align-items-center",
        style={"marginTop": "40px", "marginLeft": "20px", "gap": "5px"}),

        # Détails de la substance sélectionnée
        html.Div(id="substance-details", className="mt-4"),

        # Modal de dessin Kekule + zone SMILES
        dbc.Modal([
            dbc.ModalBody([

                # Instruction utilisateur
                html.P([
                    html.I(className="fas fa-vial me-2"),# Font Awesome
                    "You can manually draw your molecule here and Press Get Smiles."
                ], style={"fontWeight": "bold", "color": "#003366"}),

                # Kekule + Affichage RDKit
                dbc.Row([
                    dbc.Col(
                        html.Iframe(
                            src="/assets/kekule_editor.html",
                            style={"width": "100%", "height": "500px", "border": "none"},
                            id="kekule-editor"
                        ),
                        width=6
                    ),
                    dbc.Col(
                        html.Div(id="display-rdkit-svg", style={"height": "500px"}),
                        width=6
                    )
                ]),

                html.Br(),

                # Input SMILES manuel
                dcc.Input(
                    id="manual-smiles",
                    placeholder="Paste canonical SMILES here Press Draw and Search to get legal information .",
                    type="text",
                    style={"width": "100%"}
                ),

                html.Br(),

                # Boutons : Get SMILES | Enter | Search | Close
                dbc.Row([
                    dbc.Col([
                        dbc.Button([
                            html.I(className="fas fa-flask me-2"),# Font Awesome
                            "Get SMILES"
                        ],
                        id="get-smiles-btn",
                        n_clicks=0,
                        className="mt-2 me-2",
                        style={"background-color": "#066ab1"},
                        size="sm"),

                        dbc.Button([
                            html.I(className="fas fa-play me-2"),# Font Awesome
                            "draw"
                        ],
                        id="load-smiles",
                        n_clicks=0,
                        className="mt-2 me-2",
                        style={"background-color": "#066ab1"},
                        size="sm"),

                        dbc.Button([
                            html.I(className="fas fa-search me-2"),  # Font Awesome
                            "Search"
                        ],
                        id="modal-search-button",
                        n_clicks=0,
                        className="mt-2",
                        style={"background-color": "#066ab1"},
                        color="warning",
                        size="sm")
                    ], width="auto"),

                    dbc.Col(
                        dbc.Button([
                            html.I(className="fas fa-times me-2"),  # Font Awesome
                            "Close"
                        ],
                        id="close-draw",
                        n_clicks=0,
                        className="mt-2",
                        style={"background-color": "#066ab1", "color": "white"},
                        size="sm"),
                        className="text-end"
                    )
                ]),

                # Stocke les SMILES (bruts) récupérés de Kekule
                dcc.Store(id="raw-smiles"),

                # Résultats du bouton "Enter" (affichage dynamique)
                html.Div(id="modal-search-results", className="mt-4")
            ])
        ],
        id="modal-draw",
        is_open=False,
        size="xl",
        fullscreen=True)

    ],
    style={
        "display": display,
        "margin-left": "-70px",   # Décalage global
        "padding-left": "70px",   # Compensation du décalage
        "width": "calc(100% + 20px)"
    })
