from dash import Dash, Input, Output, html, State
import dash_bootstrap_components as dbc
from Main_window import get_main_layout,df # Import de l'interface principale et de la base de données

# Application Dash avec un thème Bootstrap
external_stylesheets = [dbc.themes.CERULEAN]
app = Dash(__name__, external_stylesheets=external_stylesheets)


app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            get_main_layout() # Appel de la fonction qui crée l’interface principale 
        ])
    ])
], fluid=True)# fluid=True pour que le layout s’adapte à la taille de l’écran

# Callback pour afficher les détails d'une substance saisie dans la zone de recherche
@app.callback(
    Output("substance-details", "children"),
    Input("input-substance", "value")
)
def display_info(substance):
    if not substance:
        return "" # Ne rien afficher si aucun texte n’est entré

    # Recherche exacte de la substance dans la base de données
    results = df[df['substance_name'].str.lower() == substance.lower()]

    if results.empty:
        return dbc.Alert("Substance not found.", color="danger", className="mt-2") # Alerte si non trouvée
    
# Récupère la première ligne correspondante
    ligne = results.iloc[0]
    nom = ligne['substance_name']

# Affiche les informations de la substance trouvée
    return dbc.Col([
        dbc.Row(html.H5(f"Here are the details of {nom} :", className="text-info")),
        html.Ul([ # Liste des propriétés
            html.Li(f"CAS number : {ligne['cas_number']}"),
            html.Li(f"SMILES : {ligne['smiles']}"),
            html.Li(f"Legal status : {ligne['legal_status']}"),
            html.Li(f"Convention : {ligne['convention']}"),
            html.Li(f"Schedule or Table : {ligne['schedule_or_table']}"),
        ], className="ms-3")
    ])
#  Callback pour afficher ou masquer la zone de dessin (modal)
@app.callback(
    Output("modal-draw", "is_open"),# ouverture ou fermeture du modal
    Input("draw-button", "n_clicks"),# Clic sur le bouton "Draw"
    Input("close-draw", "n_clicks"),# Clic sur le bouton "Close"
    State("modal-draw", "is_open")# État actuel du modal (ouvert/fermé)
)
def toggle_modal(draw_clicks, close_clicks, is_modal_open):
    if close_clicks and close_clicks > 0:
        return False  # Fermer le modal(Zone de dessin) si on clique sur "Close"
    if draw_clicks and draw_clicks > 0:
        return True  # Ouvrir le modal si on clique sur "Draw"
    return is_modal_open
if __name__ == "__main__":
    app.run(debug=True)
