from dash import dcc, html
import dash_bootstrap_components as dbc
from layout.search_page import get_search_page

# Fonction qui retourne le layout global de l'application avec la sidebar
def get_main_layout():
    sidebar = html.Div([
    html.Div(
        html.Img(src="/assets/logo.png", className="logo-img"),
        className="logo-container"
    ),
    html.Div("SafeLab", className="logo-text"),
    html.Hr(),
    dbc.Nav([
        # Lien vers la page de recherche
    dbc.NavLink([
        html.I(className="fas fa-search nav-icon"),
        html.Span("Search",className="nav-text")
    ], href="/search", active="exact", className="navlink"),
        # Lien vers la base de données
    dbc.NavLink([
        html.I(className="fas fa-database", style={"marginRight": "10px"}),
        html.Span("Database",className="nav-text")
    ], href="/database", active="exact", className="navlink")
], vertical=True, pills=True)
], className="sidebar sidebar-open")  
 


    return dbc.Row([
        dbc.Col(sidebar, width=2),
        dbc.Col([
            dcc.Location(id="url", refresh=False),# Pour la navigation
            get_search_page(hidden=True), # Page de recherche (initialement cachée)
            html.Div(id="page-content")# Contenu dynamique de la page
        ], width=10)
    ])