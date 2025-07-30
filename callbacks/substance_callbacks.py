from dash import Input, Output, State, html, ClientsideFunction
import dash_bootstrap_components as dbc
from rdkit import Chem

# Internal imports
from utils.Prec_procedure import get_legal_procedure
from utils.Narc_psy_procedure import get_stup_psy_procedure
from utils.draw_rdkit import draw_substructure_svg
from layout.database_page import get_database_page
from layout.search_page import get_search_page

# Callback registration function
# ===========================================
def register_callbacks(app, df):

    # Result card generator
    # ========================================
    def generate_result_card(search, is_modal=False):
        match = df[
            (df['name'].str.lower() == search.lower()) |
            (df['cas'].str.lower() == search.lower()) |
            (df['smiles'].str.lower() == search.lower()) |
            (df['iupac'].str.lower() == search.lower()) |
            (df['smiles_canonique'].str.lower() == search.lower())
        ]

        if match.empty:
            return html.Div("No molecule found.", style={"color": "red"})

        row = match.iloc[0]

        # Secure access to regulatory row (via CAS)
        match_reg = df[df['cas'] == row['cas']]
        reg = match_reg.iloc[0] if not match_reg.empty else row

        # Display legal procedure (narcotics OR precursors)
        legal_info = get_stup_psy_procedure(reg) or get_legal_procedure(reg)

        return dbc.Card([
            dbc.CardHeader(
                f"Here are the details of {row['name']} :",
                style={"backgroundColor": "#ffcb47", "color": "#000000"}
            ),
            dbc.CardBody(
                dbc.Row([
                    # Chemical info column
                    dbc.Col([
                        html.H5("Chemical Information", className="mb-3"),
                        dbc.ListGroup([
                            dbc.ListGroupItem(f"Name : {row['name']}"),
                            dbc.ListGroupItem(f"Chemical Name (IUPAC) : {row['iupac']}"),
                            dbc.ListGroupItem(f"CAS number : {row['cas']}"),
                            dbc.ListGroupItem(f"SMILES : {row['smiles']}"),
                            dbc.ListGroupItem(f"UNO_Convention : {row['uno']}"),
                            dbc.ListGroupItem(f"Convention : {row['convention']}"),
                            dbc.ListGroupItem(f"EU : {row['eu']}"),
                            dbc.ListGroupItem(f"EU Annex/Category 273/2004 : {row['eu_annex/category_273/2004']}"),
                            dbc.ListGroupItem(f"EU Annex/Category 111/2005 : {row['eu_annex/category_111/2005']}"),
                            dbc.ListGroupItem(f"Belgian Annex : {row['belgian_annex']}"),
                            dbc.ListGroupItem(f"Belgian Code : {row['belgian_code']}"),
                        ])
                    ], md=6),

                    # Legal info column
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

    # Callbacks
    # ============================

    # Callback: search via main dropdown
    @app.callback(
        Output("substance-details", "children"),
        Input("input-substance", "value")
    )
    def update_from_main_input(input_value):
        if input_value:
            return generate_result_card(input_value, is_modal=False)
        return ""

    # Callback: search from modal after "Search" button
    @app.callback(
        Output("modal-search-results", "children"),
        Input("modal-search-button", "n_clicks"),
        State("manual-smiles", "value")
    )
    def update_from_modal(n_clicks, smiles_value):
        if n_clicks and smiles_value:
            return generate_result_card(smiles_value, is_modal=True)
        return ""

    # Callback: open / close the modal
    @app.callback(
        Output("modal-draw", "is_open"),
        [Input("draw-button", "n_clicks"),
         Input("close-draw", "n_clicks")],
        State("modal-draw", "is_open")
    )
    def toggle_modal(n_draw, n_info, is_open):
        if n_draw or n_info:
            return not is_open
        return is_open

    # Callback: page change (routing)
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

    # Callback: render SVG from SMILES in Kekule input
    @app.callback(
        Output("display-rdkit-svg", "children"),
        Input("load-smiles", "n_clicks"),
        State("manual-smiles", "value")
    )
    def draw_from_editor(n_clicks, smiles):
        if not smiles:
            return ""
        svg = draw_substructure_svg(smiles, smarts=smiles)
        return html.Iframe(srcDoc=svg, style={"width": "100%", "height": "500px", "border": "none"})

    # Client-side JS callback: retrieve SMILES from Kekule
    app.clientside_callback(
        ClientsideFunction(namespace="clientside", function_name="getKekuleSmiles"),
        Output("raw-smiles", "data"),
        Input("get-smiles-btn", "n_clicks")
    )

    # Canonicalize SMILES with RDKit (for SMILES generated by Kekule)
    def canonicalize_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ''
        return Chem.MolToSmiles(mol, canonical=True)

    # Callback: inject canonical SMILES after JS retrieval
    @app.callback(
        Output("manual-smiles", "value"),
        Input("raw-smiles", "data")
    )
    def update_manual_smiles(smiles):
        if not smiles:
            return ""
        return canonicalize_smiles(smiles)
