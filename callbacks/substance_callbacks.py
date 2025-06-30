from dash import Input, Output, State, html
import dash_bootstrap_components as dbc
from utils.data_loader import df
from layout.database_page import get_database_page
from layout.search_page import get_search_page


def register_callbacks(app,df):
    # Affichage des détails d'une substance sélectionnée
    @app.callback(
        Output("substance-details", "children"),
        Input("input-substance", "value")
    )
    def display_info(substance):
        if not substance:
            return ""

        substance_lower = substance.lower()
        results = df[
            (df['substance_name'].str.lower() == substance_lower) |
            (df['cas_number'].str.lower() == substance_lower) |
            (df['smiles'].str.lower() == substance_lower)
        ]

        if results.empty:
            return dbc.Alert("Substance not found.", color="danger", className="mt-2")

        row = results.iloc[0]
        return dbc.Card([
            dbc.CardHeader(f"Here are the details of {row['substance_name']} :", style={"backgroundColor": "#ffcb47","color":"#000000"}),
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
    # Ouverture/fermeture de la fenêtre modale "Draw"
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
    # Gestion du routage entre les pages (search / database)
    @app.callback(
    Output("page-content", "children"),
    Input("url", "pathname")
    )
    def display_page(pathname):
       if pathname in ["/", "/search"]:
        return get_search_page()
       elif pathname == "/database":
        return get_database_page()
       return html.Div("Page not found.", className="p-4")
