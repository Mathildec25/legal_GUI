# search_page.py
from dash import html, dcc
import dash_bootstrap_components as dbc
from utils.data_loader import df
import pandas as pd

# Options pour le dropdown à partir de plusieurs colonnes
options = []
for col in ['substance_name', 'cas_number', 'smiles']:
    for val in df[col]:
        if pd.notna(val):
            options.append({"label": val, "value": val})
 
# Fonction qui retourne le layout de la page de recherche
def get_search_page(hidden=False):
    display = "none" if hidden else "block"
    return html.Div([
        # Barre de recherche et boutons
        dbc.Row([
            dbc.Col(
                dbc.Button([
           html.I(className="fas fa-pencil-alt", style={"marginRight": "6px"}),
               "Draw"
             ], id="draw-button", n_clicks=0,
              style={"background-color": "#066ab1", "color": "white"})
,
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
                    html.I(className="fas fa-search", style={"marginRight": "6px"}),
                        "Search"
                    ], id="search-button", n_clicks=0,
                     style={"background-color": "#066ab1", "color": "white"})
                       ,
                     width="auto"
                    )
                    ], className="align-items-center", style={"marginTop": "40px", "marginLeft": "20px", "gap": "5px"}),
            # Détails de la substance recherchée
            html.Div(id="substance-details", className="mt-4"),
            
            # Modal pour dessiner une structure moléculaire
            dbc.Modal([
                dbc.ModalHeader(dbc.ModalTitle("Draw Structure")),
                dbc.ModalBody(html.Div("Drawing area", style={
                "border": "2px dashed #aaa", "height": "200px",
                "width": "100%", "backgroundColor": "#F0F2F5"
            })),
            dbc.ModalFooter(
                dbc.Button("Close", id="close-draw", className="ms-auto", n_clicks=0,style={"background-color": "#066ab1"})
            )
        ], id="modal-draw", is_open=False, size="lg")
    ], style={"display": display})
