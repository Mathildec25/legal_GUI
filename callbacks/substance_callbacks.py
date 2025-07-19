from dash import Input, Output, State, html
import dash_bootstrap_components as dbc
from utils.data_loader import df
from layout.database_page import get_database_page
from layout.search_page import get_search_page
from utils.legal_procedure import get_legal_procedure
from utils.Narc_psy_procedure import get_stup_psy_procedure


def register_callbacks(app, df):

    # 📌 Callback pour afficher les détails d'une substance
    @app.callback(
        Output("substance-details", "children"),
        Input("input-substance", "value")
    )
    def display_info(substance):
        if not substance:
            return ""

        substance_lower = substance.lower()
        results = df[
            (df['name'].str.lower() == substance_lower) |
            (df['cas'].str.lower() == substance_lower) |
            (df['smiles'].str.lower() == substance_lower)
        ]

        if results.empty:
            return dbc.Alert("Substance not found.", color="danger", className="mt-2")

        row = results.iloc[0]

        # 🔍 Recherche réglementaire
        reg_row = df[df['cas'] == row['cas']]
        if not reg_row.empty:
            reg = reg_row.iloc[0]

            # On regarde si on traite un stupéfiant/psychotrope
            stup_psy_info = get_stup_psy_procedure(reg)

            if stup_psy_info:
                legal_info = stup_psy_info
            else:
                legal_info = get_legal_procedure(reg)
        else:
            legal_info = "No regulatory information available."

        return dbc.Card([
            dbc.CardHeader(
                f"Here are the details of {row['name']} :",
                style={"backgroundColor": "#ffcb47", "color": "#000000"}
            ),
            dbc.CardBody(
                dbc.Row([
                    # Bloc chimique
                    dbc.Col([
                        html.H5("Chemical Information", className="mb-3"),
                        dbc.ListGroup([
                            dbc.ListGroupItem(f"CAS number : {row['cas']}"),
                            dbc.ListGroupItem(f"SMILES : {row['smiles']}"),
                            dbc.ListGroupItem(f"UNO : {row['uno']}"),
                            dbc.ListGroupItem(f"Convention : {row['convention']}"),
                            dbc.ListGroupItem(f"Belgian Annex : {row['belgian_annex']}"),
                        ])
                    ], md=6),

                    # Bloc juridique
                    dbc.Col([
                        html.H5("Legal Information", className="mb-3"),
                        dbc.Card(
                            dbc.CardBody(
                                html.Pre(legal_info, style={"whiteSpace": "pre-wrap"})
                            ),
                            color="light",
                            className="shadow-sm"
                        )
                    ], md=6)
                ])
            )
        ], className="mt-3")

    # 📌 Callback pour ouvrir/fermer la modale "Draw"
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

    # 📌 Callback pour gérer la navigation entre les pages
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
