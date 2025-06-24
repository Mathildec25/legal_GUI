from dash import Input, Output, State, dash_table
import dash_bootstrap_components as dbc
from utils.data_loader import df

def register_callbacks(app, df):
    @app.callback(
        Output("substance-details", "children"),
        Input("input-substance", "value")
    )
    def display_info(substance):
        if not substance:
            return ""

        results = df[df['substance_name'].str.lower() == substance.lower()]
        if results.empty:
            return dbc.Alert("Substance not found.", color="danger", className="mt-2")

        row = results.iloc[0]
        return dbc.Card([
            dbc.CardHeader(f"Here are the details of {row['substance_name']} :"),
            dbc.CardBody(
                dbc.ListGroup([
                    dbc.ListGroupItem(f"CAS number : {row['cas_number']}"),
                    dbc.ListGroupItem(f"SMILES : {row['smiles']}"),
                    dbc.ListGroupItem(f"Legal status : {row['legal_status']}"),
                    dbc.ListGroupItem(f"Convention : {row['convention']}"),
                    dbc.ListGroupItem(f"Schedule or Table : {row['schedule_or_table']}"),
                ])
            )
        ], className="mt-3")

    @app.callback(
        Output("modal-draw", "is_open"),
        Input("draw-button", "n_clicks"),
        Input("close-draw", "n_clicks"),
        State("modal-draw", "is_open")
    )
    def toggle_modal(n_draw, n_close, is_open):
        if n_close and n_close > 0:
            return False
        if n_draw and n_draw > 0:
            return True
        return is_open

    @app.callback(
    Output("modal-database", "is_open"),
    Input("open-database", "n_clicks"),
    Input("close-database", "n_clicks"),
    State("modal-database", "is_open"),
)
    def toggle_database_modal(n_open, n_close, is_open):
         if n_close:
            return False
         if n_open:
            return True
         return is_open
