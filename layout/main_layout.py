from dash import dcc, html
import dash_bootstrap_components as dbc
from layout.search_page import get_search_page

# Fonction principale du layout
# ===========================
# Cette fonction retourne la structure complète de l'application Dash
# avec une sidebar à gauche et un contenu dynamique à droite.

def get_main_layout():

    #
    #  Sidebar avec navigation
    # ===========================

    sidebar = html.Div([
        # Logo de l'application (icône)
        html.Div(
            html.Img(src="assets/favicon.ico", className="logo-img"),
            className="logo-container"
        ),

        # Titre de l'application
        html.Div("SafeLab", className="logo-text"),

        html.Hr(),  # Séparation visuelle

        #  Navigation verticale (Font Awesome + texte)
        dbc.Nav([
            # Lien vers la page de recherche
            dbc.NavLink([
                html.I(className="fas fa-search nav-icon"),
                html.Span("Search", className="nav-text")
            ],
            href="/search",
            active="exact",
            className="navlink"),

            # Lien vers la base de données
            dbc.NavLink([
                html.I(className="fas fa-database", style={"marginRight": "10px"}),
                html.Span("Database", className="nav-text")
            ],
            href="/database",
            active="exact",
            className="navlink")
        ],
        vertical=True,
        pills=True)
    ],
    className="sidebar sidebar-open")

    # Structure de page principale (2 colonnes : sidebar + contenu)
    # ===========================

    return dbc.Row([
        # Colonne gauche : sidebar (2/12)
        dbc.Col(sidebar, width=2),

        # Colonne droite : contenu dynamique (10/12)
        dbc.Col([
            dcc.Location(id="url", refresh=False),              # Gère l'URL pour la navigation
            get_search_page(hidden=True),                       # Charge la page Search mais masquée par défaut
            html.Div(id="page-content")                         # Affiche dynamiquement le contenu selon l'URL
        ], width=10)
    ],
    className="g-0")  # "g-0" supprime l'espacement horizontal entre les colonnes Bootstrap
