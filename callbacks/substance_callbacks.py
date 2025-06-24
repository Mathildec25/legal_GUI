from dash import Input, Output, State
import dash_bootstrap_components as dbc
from utils.data_loader import load_clean_data

df = load_clean_data()

def register_main_callbacks(app):
    @app.callback(
        Output("substance-details", "children"),
        Input("input-substance", "value")
    )
    def display_info(substance):
        if not substance:
            return ""

        resultats = df[df['substance_name'].str.lower() == substance.lower()]
        if resultats.empty:
            return dbc.Alert("Substance not found.", color="danger", className="mt-2")

        ligne = resultats.iloc[0]
        nom = ligne['substance_name']

        return dbc.Container([
            dbc.Container(f"Here are the details of {nom} :", className="text-info h5 mb-3"),
            dbc.ListGroup([
                dbc.ListGroupItem(f"CAS number : {ligne['cas_number']}"),
                dbc.ListGroupItem(f"SMILES : {ligne['smiles']}"),
                dbc.ListGroupItem(f"Legal status : {ligne['legal_status']}"),
                dbc.ListGroupItem(f"Convention : {ligne['convention']}"),
                dbc.ListGroupItem(f"Schedule or Table : {ligne['schedule_or_table']}")
            ], flush=True)
        ])

    @app.callback(
        Output("collapse-table", "is_open"),
        Input("menu-show-table", "n_clicks"),
        State("collapse-table", "is_open")
    )
    def show_or_hide_table(n_clicks, is_open):
        if n_clicks:
            return not is_open
        return is_open
