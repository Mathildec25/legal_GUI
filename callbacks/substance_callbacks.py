from dash import Input, Output, State, html, ClientsideFunction
import dash_bootstrap_components as dbc

# Internal imports
from utils.Prec_procedure import get_legal_procedure
from utils.Narc_psy_procedure import get_stup_psy_procedure
from utils.chem_cache import mol_from_smiles, svg_from_smiles  # cache LRU
from layout.database_page import get_database_page
from layout.search_page import get_search_page


def register_callbacks(app, df):
    """Register all callbacks for the substance search and display."""

    # Choisir le bon bloc légal
    def pick_legal_block(row):
        annex = str(row.get('belgian_annex', '')).strip().lower()
        cat_273 = str(row.get('eu_annex/category_273/2004', '')).strip().lower()
        cat_111 = str(row.get('eu_annex/category_111/2005', '')).strip().lower()

        if annex:  # stup/psy
            return get_stup_psy_procedure(row)
        if cat_273 or cat_111:  # précurseurs
            return get_legal_procedure(row)
        return html.Div("Aucune procédure identifiée pour cet enregistrement.", className="text-muted")

    # Recherche tolérante
    def find_matches(search):
        if not search:
            return df.iloc[0:0]  # empty DataFrame

        s = str(search).strip().lower()

        # Exact match sur name, cas, smiles, iupac, canonical_smiles
        exact = df[
            df['name'].str.lower() == s
            | df['cas'].str.lower() == s
            | df['smiles'].str.lower() == s
            | df['iupac'].str.lower() == s
            | df['canonical_smiles'].str.lower() == s
        ]
        if not exact.empty:
            return exact

        # Fallback contains sur name et iupac
        return df[df['name'].str.lower().str.contains(s, na=False)
                  | df['iupac'].str.lower().str.contains(s, na=False)]

    # Carte résultat
    def generate_result_card(search, is_modal=False):
        matches = find_matches(search)
        if matches.empty:
            return html.Div("No molecule found.", style={"color": "red"})

        row = matches.iloc[0]

        # Sécurisation via CAS
        cas_val = str(row.get('cas', '')).strip()
        match_reg = df[df['cas'].astype(str).str.strip() == cas_val] if cas_val else matches
        reg = match_reg.iloc[0] if not match_reg.empty else row

        # Bloc légal
        legal_info = pick_legal_block(reg)

        return dbc.Card([
            dbc.CardHeader(
                f"Here are the details of {row.get('name','(unknown)')} :",
                style={"backgroundColor": "#ffcb47", "color": "#000000"}
            ),
            dbc.CardBody(
                dbc.Row([
                    dbc.Col([
                        html.H5("Chemical Information", className="mb-3"),
                        dbc.ListGroup([
                            dbc.ListGroupItem(f"Name : {row.get('name','')}"),
                            dbc.ListGroupItem(f"Chemical Name (IUPAC) : {row.get('iupac','')}"),
                            dbc.ListGroupItem(f"CAS number : {row.get('cas','')}"),
                            dbc.ListGroupItem(f"SMILES : {row.get('smiles','')}"),
                            dbc.ListGroupItem(f"UNO_Convention : {row.get('uno','')}"),
                            dbc.ListGroupItem(f"Convention : {row.get('convention','')}"),
                            dbc.ListGroupItem(f"EU : {row.get('eu','')}"),
                            dbc.ListGroupItem(f"EU Annex/Category 273/2004 : {row.get('eu_annex/category_273/2004','')}"),
                            dbc.ListGroupItem(f"EU Annex/Category 111/2005 : {row.get('eu_annex/category_111/2005','')}"),
                            dbc.ListGroupItem(f"Belgian Annex : {row.get('belgian_annex','')}"),
                            dbc.ListGroupItem(f"Belgian Code : {row.get('belgian_code','')}"),
                        ])
                    ], md=6),

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

    # ============================
    # Callbacks
    # ============================

    @app.callback(
        Output("substance-details", "children"),
        Input("input-substance", "value")
    )
    def update_from_main_input(input_value):
        if input_value:
            return generate_result_card(input_value, is_modal=False)
        return ""

    @app.callback(
        Output("modal-search-results", "children"),
        Input("modal-search-button", "n_clicks"),
        State("manual-smiles", "value")
    )
    def update_from_modal(n_clicks, smiles_value):
        if n_clicks and smiles_value:
            return generate_result_card(smiles_value, is_modal=True)
        return ""

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

    @app.callback(
        Output("display-rdkit-svg", "children"),
        Input("load-smiles", "n_clicks"),
        State("manual-smiles", "value")
    )
    def draw_from_editor(n_clicks, smiles):
        if not smiles:
            return ""
        mol = mol_from_smiles(smiles)
        if mol is None:
            return html.Div("Invalid SMILES.", style={"color": "red"})
        svg = svg_from_smiles(smiles)
        return html.Iframe(srcDoc=svg, style={"width": "100%", "height": "500px", "border": "none"})

    app.clientside_callback(
        ClientsideFunction(namespace="clientside", function_name="getKekuleSmiles"),
        Output("raw-smiles", "data"),
        Input("get-smiles-btn", "n_clicks")
    )

    @app.callback(
        Output("manual-smiles", "value"),
        Input("raw-smiles", "data")
    )
    def update_manual_smiles(smiles):   
        if not smiles:
            return ""
        return canonicalize_smiles(smiles) or ""
