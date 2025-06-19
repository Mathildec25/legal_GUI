from dash import html, dcc, Input, Output, State, callback, dash_table
import pandas as pd
import dash_bootstrap_components as dbc

# Charger la base de données 
df = pd.read_csv("data/donnees.csv", sep=";", encoding="utf-8")

# Renommer les colonnes pour éviter les espaces
df.columns = [col.lower().replace(" ", "").strip() for col in df.columns]

def get_main_layout():
    return html.Div([
        html.H2("Main Interface", className="text-success text-center mt-4"),

        html.Div([
            html.Label("Please enter a substance:", className="fw-bold mb-2"),
        # Zone de saisie pour la recherche d'informations
            dcc.Input(
                id="input-substance",
                type="text",
                placeholder="e.g., acetone",
                debounce=True,
                className="form-control"
            ),
        ], className="mb-4"),
       #Zone d'affichage des résultats
        html.Div(id="substance-details", className="mt-3"),

        html.H4("Draw your structure here", className="text-primary mt-4"),
       #Zone pour dessiner la molécule
        html.Div(className="border border-2 border-secondary rounded w-100 h-100", id="drawing-box"),
       #Bouton pour retourner à la page d'accueil
         dbc.Button("Back to Homepage", id="back-button", color="secondary", className="mt-4"),

        html.Hr(),

        dbc.DropdownMenu(
            label="Menu",
            children=[
                dbc.DropdownMenuItem("View the database", id="menu-show-table")
            ],
            color="info",
            className="mb-3"
        ),

        dbc.Collapse(
            dash_table.DataTable(
                id='table-database',
                columns=[{"name": col, "id": col} for col in df.columns],
                data=df.to_dict('records'),
                page_size=10,
                style_table={
                    'maxHeight': '300px',
                    'overflowY': 'auto',
                    'overflowX': 'auto'
                },
                style_cell={
                    'textAlign': 'left',
                    'padding': '5px'
                },
                style_header={
                    'backgroundColor': '#f8f9fa',
                    'fontWeight': 'bold'
                }
            ),
            id="collapse-table",
            is_open=False
        )
    ], className="container")

@callback(
    Output("substance-details", "children"),
    Input("input-substance", "value")
)
# Fonction pour afficher les résultats de la molécule saisie dans la zone de recherche
def display_info(substance):
    if not substance:
        return ""

    resultats = df[df['substance_name'].str.lower() == substance.lower()]

    if resultats.empty:
        return dbc.Alert("Substance not found.", color="danger", className="mt-2")
# sélectionne la première ligne de resultats
    ligne = resultats.iloc[0]
    nom = ligne['substance_name']

    return html.Div([
        html.H5(f"Here are the details of {nom} :", className="text-info"),
        html.Ul([
            html.Li(f"CAS number : {ligne['cas_number']}"),
            html.Li(f"SMILES : {ligne['smiles']}"),
            html.Li(f"Legal status : {ligne['legal_status']}"),
            html.Li(f"Convention : {ligne['convention']}"),
            html.Li(f"Schedule or Table : {ligne['schedule_or_table']}")
        ], className="ms-3")
    ])

@callback(
    Output("collapse-table", "is_open"),
    Input("menu-show-table", "n_clicks"),
    State("collapse-table", "is_open")
)
# Fonction pour afficher ou cacher la base de donnée
def  show_or_hide_table(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open
