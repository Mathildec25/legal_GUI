from dash import html, dcc
import dash_bootstrap_components as dbc
import pandas as pd

#Charger la base de données
df = pd.read_csv("data/donnees.csv", sep=";", encoding="utf-8")
# Nettoyer les noms de colonnes : tout en minuscules, sans espaces
df.columns = [col.lower().replace(" ", "_").strip() for col in df.columns]

# Fonction principale pour construire l'interface utilisateur
def get_main_layout():
# Barre de recherche avec bouton Draw, champ de saisie et bouton Search
    search_bar = dbc.Row(
        [
            # Bouton pour ouvrir la zone de dessin
            dbc.Col(
                dbc.Button("Draw", color="success", className="ms-2", id="draw-button"),
                width="auto",
                style={"marginTop": "35px","marginLeft": "30px"}
            ),
            # Zone de saisie du nom de la substance
            dbc.Col(
                dcc.Input(
                    id="input-substance",
                    type="text",
                    placeholder="Enter substance name...",
                    debounce=True,
                    className="form-control"
                ),
                style={"marginTop": "35px", "marginLeft": "10px", "width": "90px"}
            ),
            # Bouton de recherche
            dbc.Col(
                dbc.Button("Search", color="success", className="ms-2", id="search-button", n_clicks=0),
                width="auto",
                style={"marginTop": "35px","marginLeft": "30px"}
            )
        ],
        className="g-0 flex-nowrap mt-3 mt-md-0",
        align="center"
    )
# Barre latérale avec menu vertical
    sidebar = dbc.Nav([
        html.H2("Menu", style={"color": "white", "fontSize": "36px"}), # Titre du menu
        html.Hr(),# Ligne de séparation
        dbc.NavLink("DataBase", href="#", active=True, style={"color": "white","backgroundColor":"#1E5631"}),
        dbc.NavLink("Help", href="#", active=True, style={"color": "white","backgroundColor":"#1E5631"})
    ],
        vertical=True,
        pills=True,
        className="sidebar" #  Classe personnalisée définie dans sidebar.css
    )
 # Fenêtre modale pour dessiner une molécule
    draw_modal = dbc.Modal(
        [
            dbc.ModalHeader("Draw Structure"), # Titre de la modale
            dbc.ModalBody( # Zone vide 
                html.Div("Drawing area", style={
                    "border": "2px dashed #aaa",
                    "height": "200px",
                    "width": "100%",
                    "backgroundColor": "#F0F2F5"
                })
            ),
            dbc.ModalFooter( 
                dbc.Button("Close", id="close-draw", className="ms-auto", n_clicks=0, color="success")
            )
        ],
        id="modal-draw",
        is_open=False,# Fermé par défaut
        size="lg" # Grande taille
    )
# Retourne l'organisation de la page en deux colonnes : sidebar et contenu principal
    return dbc.Row([
        dbc.Col(sidebar, width=2),# Colonne pour la barre latérale
        dbc.Col([# Colonne principale
            search_bar,
            html.Div(id="substance-details", className="mt-4"),
            draw_modal
        ], width=10)
    ])
