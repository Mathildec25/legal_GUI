import pandas as pd
from dash import Input, Output, State, callback, no_update
import dash_bootstrap_components as dbc
from utils.data_loader import load_clean_data

df = load_clean_data()

def register_main_callbacks(app):
    @app.callback(
        Output("modal-draw", "is_open"),
        [Input("draw-button", "n_clicks"), Input("close-draw", "n_clicks"), Input("confirm-draw", "n_clicks")],
        [State("modal-draw", "is_open")]
    )
    def toggle_modal(n1, n2, n3, is_open):
        if n1 or n2 or n3:
            return not is_open
        return is_open

    @app.callback(
        Output("substance-details", "children"),
        Input("store-drawn-smiles", "data"),
        prevent_initial_call=True
    )
    def display_info_from_drawn_smiles(smiles):
        if not smiles:
            return dbc.Alert("No SMILES data received.", color="warning")

        resultats = df[df['smiles'].str.strip() == smiles.strip()]
        if resultats.empty:
            return dbc.Alert("No match found in database.", color="danger")

        ligne = resultats.iloc[0]
        return dbc.Container([
            dbc.Container(f"Here are the details of {ligne['substance_name']} :", className="text-info h5 mb-3"),
            dbc.ListGroup([
                dbc.ListGroupItem(f"CAS number : {ligne['cas_number']}"),
                dbc.ListGroupItem(f"SMILES : {ligne['smiles']}"),
                dbc.ListGroupItem(f"Legal status : {ligne['legal_status']}"),
                dbc.ListGroupItem(f"Convention : {ligne['convention']}"),
                dbc.ListGroupItem(f"Schedule or Table : {ligne['schedule_or_table']}")
            ])
        ])
