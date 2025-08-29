from dash import Input, Output, State, html, ClientsideFunction, callback_context
import dash_bootstrap_components as dbc
from rdkit import Chem
import pandas as pd
import json

# Internal imports
from utils.Prec_procedure import get_legal_procedure
from utils.Narc_psy_procedure import get_stup_psy_procedure
from utils.draw_rdkit import draw_substructure_svg
from layout.database_page import get_database_page
from layout.search_page import get_search_page

def register_callbacks(app, df):
    """
    Callbacks refactorisés pour la nouvelle interface moderne
    """
    
    # ===================================
    # Générateur de cartes de résultats moderne
    # ===================================
    def generate_modern_result_card(search, is_modal=False):
        """Génère une carte de résultat avec le nouveau design"""
        match = df[
            (df['name'].str.lower() == search.lower()) |
            (df['cas'].str.lower() == search.lower()) |
            (df['smiles'].str.lower() == search.lower()) |
            (df['iupac'].str.lower() == search.lower()) |
            (df['canonical_smiles'].str.lower() == search.lower())
        ]

        if match.empty:
            return create_no_results_card(search)

        row = match.iloc[0]
        
        # Accès sécurisé aux données réglementaires
        match_reg = df[df['cas'] == row['cas']]
        reg = match_reg.iloc[0] if not match_reg.empty else row

        # Informations légales
        legal_info = get_stup_psy_procedure(reg) or get_legal_procedure(reg)

        return create_molecule_card(row, reg, legal_info, is_modal)
    
    def create_no_results_card(search):
        """Carte moderne pour aucun résultat"""
        return dbc.Card([
            dbc.CardBody([
                html.Div([
                    html.I(className="fas fa-search-minus", 
                          style={"fontSize": "3rem", "color": "#94a3b8", "marginBottom": "1rem"}),
                    html.H4("No Results Found", className="text-muted"),
                    html.P(f"No molecule found matching '{search}'", className="text-muted"),
                    dbc.Alert([
                        html.I(className="fas fa-lightbulb me-2"),
                        html.Strong("Search Tips:"),
                        html.Ul([
                            html.Li("Try searching by exact chemical name"),
                            html.Li("Use the complete CAS number format (e.g., 50-00-0)"),
                            html.Li("Verify your SMILES notation is correct"),
                            html.Li("Try the molecular drawing tool for complex structures")
                        ], className="mb-0 mt-2")
                    ], color="info", className="mt-3"),
                    dbc.Button([
                        html.I(className="fas fa-pencil-alt me-2"),
                        "Try Drawing Instead"
                    ], 
                    color="primary", 
                    outline=True,
                    id="try-drawing-btn",
                    className="mt-2")
                ], className="text-center p-4")
            ])
        ], className="result-card no-results-card fade-in")

    def create_molecule_card(row, reg, legal_info, is_modal=False):
        """Crée une carte molécule moderne"""
        
        # Détermination du statut réglementaire
        regulatory_status = determine_regulatory_status(reg)
        
        return dbc.Card([
            # Header avec gradient moderne
            dbc.CardHeader([
                html.Div([
                    html.H3([
                        html.I(className="fas fa-atom me-2"),
                        row['name']
                    ], className="molecule-title text-white mb-1"),
                    html.P(
                        f"CAS: {row['cas']} | Formula: {extract_molecular_formula(row)}",
                        className="molecule-subtitle text-white-50 mb-0"
                    )
                ])
            ], className="molecule-header-modern"),
            
            dbc.CardBody([
                # Badges de statut en haut
                html.Div([
                    create_status_badges(reg)
                ], className="status-badges-container mb-4"),
                
                # Contenu principal en colonnes
                dbc.Row([
                    # Colonne informations chimiques
                    dbc.Col([
                        html.H5([
                            html.I(className="fas fa-flask me-2 text-primary"),
                            "Chemical Information"
                        ], className="section-title"),
                        
                        create_modern_info_list([
                            ("Chemical Name (IUPAC)", row['iupac']),
                            ("CAS Number", row['cas']),
                            ("SMILES", row['smiles'], "monospace"),
                            ("Canonical SMILES", row['canonical_smiles'], "monospace"),
                            ("UNO Convention", row['uno']),
                            ("Convention Status", row['convention'])
                        ])
                    ], lg=6),
                    
                    # Colonne informations réglementaires
                    dbc.Col([
                        html.H5([
                            html.I(className="fas fa-gavel me-2 text-primary"),
                            "Regulatory Information"
                        ], className="section-title"),
                        
                        create_regulatory_info_panel(reg, legal_info)
                    ], lg=6)
                ]),
                
                # Section actions
                html.Hr(className="my-4"),
                html.Div([
                    dbc.ButtonGroup([
                        dbc.Button([
                            html.I(className="fas fa-eye me-2"),
                            "View 3D Structure"
                        ], color="primary", outline=True, size="sm"),
                        
                        dbc.Button([
                            html.I(className="fas fa-download me-2"),
                            "Export Data"
                        ], color="secondary", outline=True, size="sm"),
                        
                        dbc.Button([
                            html.I(className="fas fa-share-alt me-2"),
                            "Share"
                        ], color="info", outline=True, size="sm")
                    ], className="molecule-actions")
                ], className="text-center")
            ])
        ], className="molecule-result-card fade-in shadow-lg")

    def create_status_badges(reg):
        """Crée les badges de statut colorés"""
        badges = []
        
        # Badge EU 273/2004
        eu_273 = str(reg.get('eu_annex/category_273/2004', '')).strip().lower()
        if eu_273 and eu_273 not in ['', 'nan', 'none']:
            badge_color, badge_text = get_273_badge_info(eu_273)
            badges.append(
                dbc.Badge([
                    html.I(className="fas fa-shield-alt me-1"),
                    f"EU 273/2004: {badge_text}"
                ], color=badge_color, className="me-2 mb-2 status-badge")
            )
        
        # Badge EU 111/2005
        eu_111 = str(reg.get('eu_annex/category_111/2005', '')).strip().lower()
        if eu_111 and eu_111 not in ['', 'nan', 'none']:
            badge_color, badge_text = get_111_badge_info(eu_111)
            badges.append(
                dbc.Badge([
                    html.I(className="fas fa-globe-europe me-1"),
                    f"EU 111/2005: {badge_text}"
                ], color=badge_color, className="me-2 mb-2 status-badge")
            )
        
        # Badge Belgian Annex
        belgian = str(reg.get('belgian_annex', '')).strip().lower()
        if belgian and belgian not in ['', 'nan', 'none']:
            badges.append(
                dbc.Badge([
                    html.I(className="fas fa-flag me-1"),
                    f"Belgian: {belgian.title()}"
                ], color="warning", className="me-2 mb-2 status-badge")
            )
        
        if not badges:
            badges.append(
                dbc.Badge([
                    html.I(className="fas fa-check-circle me-1"),
                    "Not Controlled"
                ], color="success", className="me-2 mb-2 status-badge")
            )
        
        return badges

    def create_modern_info_list(items):
        """Crée une liste d'informations moderne"""
        list_items = []
        for item in items:
            if len(item) == 2:
                label, value = item
                font_family = "inherit"
            else:
                label, value, font_family = item
            
            if pd.notna(value) and str(value).strip():
                list_items.append(
                    html.Div([
                        html.Small(label, className="text-muted fw-bold"),
                        html.Div(
                            str(value),
                            style={"fontFamily": font_family} if font_family != "inherit" else {},
                            className="info-value"
                        )
                    ], className="info-item mb-3")
                )
        
        return html.Div(list_items, className="modern-info-list")

    def create_regulatory_info_panel(reg, legal_info):
        """Crée le panneau d'informations réglementaires"""
        if not legal_info:
            return dbc.Alert([
                html.I(className="fas fa-info-circle me-2"),
                "No specific regulatory procedures detected. Please verify with competent authorities."
            ], color="light", className="regulatory-alert")
        
        return dbc.Card([
            dbc.CardBody([
                legal_info
            ])
        ], color="light", outline=True, className="regulatory-panel")

    # Fonctions utilitaires
    def determine_regulatory_status(reg):
        """Détermine le statut réglementaire global"""
        # Logique de détermination du statut
        pass

    def get_273_badge_info(category):
        """Retourne les infos de badge pour EU 273/2004"""
        if "category 1" in category:
            return "danger", "Cat 1 - License Required"
        elif "subcategory 2a" in category:
            return "warning", "Cat 2A - Registration"
        elif "subcategory 2b" in category:
            return "warning", "Cat 2B - Threshold"
        elif "category 3" in category:
            return "info", "Cat 3 - Vigilance"
        return "secondary", "Regulated"

    def get_111_badge_info(category):
        """Retourne les infos de badge pour EU 111/2005"""
        if "category 1" in category:
            return "danger", "Cat 1 - Import/Export Auth"
        elif "category 2" in category:
            return "warning", "Cat 2 - Export Auth"
        elif "category 3" in category:
            return "info", "Cat 3 - Limited Export"
        return "secondary", "Regulated"

    def extract_molecular_formula(row):
        """Extrait ou génère la formule moléculaire"""
        # Si on a la formule directement
        if 'molecular_formula' in row and pd.notna(row['molecular_formula']):
            return row['molecular_formula']
        
        # Sinon essayer de la générer depuis SMILES
        smiles = row.get('smiles') or row.get('canonical_smiles')
        if smiles and pd.notna(smiles):
            try:
                mol = Chem.MolFromSmiles(str(smiles))
                if mol:
                    return Chem.rdMolDescriptors.CalcMolFormula(mol)
            except:
                pass
        
        return "Formula not available"

    # ===================================
    # CALLBACKS PRINCIPAUX
    # ===================================
    
    # Callback: Navigation et routing
    @app.callback(
        [Output("page-content", "children"),
         Output("hero-section", "style")],
        Input("url", "pathname")
    )
    def display_page(pathname):
        """Gestion du routing avec masquage conditionnel du hero"""
        if pathname in ["/", "/search"]:
            return get_search_page(), {"display": "block"}
        elif pathname == "/database":
            return get_database_page(), {"display": "none"}
        elif pathname == "/analytics":
            return create_analytics_page(), {"display": "none"}
        else:
            return create_404_page(), {"display": "none"}

    # Callback: Recherche principale
    @app.callback(
        Output("substance-details", "children"),
        Input("input-substance", "value"),
        prevent_initial_call=True
    )
    def update_search_results(input_value):
        """Mise à jour des résultats de recherche avec loading state"""
        if not input_value:
            return create_welcome_card()
        
        return generate_modern_result_card(input_value, is_modal=False)

    # Callback: Gestion des méthodes de recherche
    @app.callback(
        [Output("method-text", "className"),
         Output("method-draw", "className"),
         Output("method-upload", "className"),
         Output("modal-draw", "is_open")],
        [Input("method-text", "n_clicks"),
         Input("method-draw", "n_clicks"),
         Input("method-upload", "n_clicks")],
        prevent_initial_call=True
    )
    def handle_search_methods(text_clicks, draw_clicks, upload_clicks):
        """Gestion des différentes méthodes de recherche"""
        ctx = callback_context
        if not ctx.triggered:
            return "search-method", "search-method", "search-method", False
        
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        base_class = "search-method"
        active_class = "search-method active"
        
        classes = [base_class, base_class, base_class]
        modal_open = False
        
        if triggered_id == "method-text":
            classes[0] = active_class
        elif triggered_id == "method-draw":
            classes[1] = active_class
            modal_open = True
        elif triggered_id == "method-upload":
            classes[2] = active_class
            # Ici on pourrait déclencher un file upload
        
        return classes[0], classes[1], classes[2], modal_open

    # Callback: Modal de dessin
    @app.callback(
        Output("modal-draw", "is_open", allow_duplicate=True),
        Input("close-draw", "n_clicks"),
        State("modal-draw", "is_open"),
        prevent_initial_call=True
    )
    def close_draw_modal(n_clicks, is_open):
        """Fermeture du modal de dessin"""
        if n_clicks:
            return False
        return is_open

    # Callback: Recherche depuis le modal
    @app.callback(
        Output("modal-search-results", "children"),
        Input("modal-search-button", "n_clicks"),
        State("manual-smiles", "value"),
        prevent_initial_call=True
    )
    def search_from_modal(n_clicks, smiles_value):
        """Recherche depuis le modal de dessin"""
        if n_clicks and smiles_value:
            return generate_modern_result_card(smiles_value, is_modal=True)
        return ""

    # Callback: Rendu SVG depuis SMILES
    @app.callback(
        Output("display-rdkit-svg", "children"),
        Input("load-smiles", "n_clicks"),
        State("manual-smiles", "value"),
        prevent_initial_call=True
    )
    def render_molecule_svg(n_clicks, smiles):
        """Rendu de la molécule en SVG"""
        if not smiles or not n_clicks:
            return html.Div([
                html.I(className="fas fa-molecule", 
                      style={"fontSize": "3rem", "color": "#cbd5e1"}),
                html.P("Draw molecule and click 'Preview'", 
                      className="text-muted mt-2")
            ], className="text-center d-flex flex-column align-items-center justify-content-center h-100")
        
        try:
            svg = draw_substructure_svg(smiles, smarts=smiles)
            return html.Iframe(
                srcDoc=svg,
                style={"width": "100%", "height": "100%", "border": "none"}
            )
        except Exception as e:
            return dbc.Alert([
                html.I(className="fas fa-exclamation-triangle me-2"),
                f"Error rendering molecule: {str(e)}"
            ], color="danger")

    # Client-side callback pour Kekule
    app.clientside_callback(
        ClientsideFunction(namespace="clientside", function_name="getKekuleSmiles"),
        Output("raw-smiles", "data"),
        Input("get-smiles-btn", "n_clicks")
    )

    # Callback: Canonicalisation SMILES
    @app.callback(
        Output("manual-smiles", "value"),
        Input("raw-smiles", "data"),
        prevent_initial_call=True
    )
    def update_manual_smiles(smiles_data):
        """Canonicalisation des SMILES récupérés"""
        if not smiles_data:
            return ""
        
        try:
            mol = Chem.MolFromSmiles(smiles_data)
            if mol is None:
                return smiles_data  # Retourner tel quel si invalide
            return Chem.MolToSmiles(mol, canonical=True)
        except:
            return smiles_data

    # Callback: Copie des SMILES
    app.clientside_callback(
        """
        function(n_clicks, smiles_value) {
            if (n_clicks && smiles_value) {
                navigator.clipboard.writeText(smiles_value);
                return "Copied!";
            }
            return "Copy";
        }
        """,
        Output("copy-smiles", "children"),
        Input("copy-smiles", "n_clicks"),
        State("manual-smiles", "value")
    )

    # ===================================
    # Fonctions utilitaires pour les pages
    # ===================================
    
    def create_welcome_card():
        """Carte d'accueil quand aucune recherche"""
        return dbc.Card([
            dbc.CardBody([
                html.Div([
                    html.I(className="fas fa-search", 
                          style={"fontSize": "4rem", "color": "#2563eb", "marginBottom": "1rem"}),
                    html.H4("Start Your Search", className="text-primary"),
                    html.P("Enter a molecule name, CAS number, or SMILES to get started", 
                          className="text-muted"),
                    html.Hr(),
                    html.H6("Popular Searches:", className="text-secondary"),
                    dbc.ButtonGroup([
                        dbc.Button("Caffeine", size="sm", color="outline-primary"),
                        dbc.Button("Aspirin", size="sm", color="outline-primary"),
                        dbc.Button("Glucose", size="sm", color="outline-primary"),
                    ], className="example-searches")
                ], className="text-center p-4")
            ])
        ], className="welcome-card fade-in")

    def create_analytics_page():
        """Page Analytics (placeholder)"""
        return html.Div([
            dbc.Container([
                html.H2("Analytics Dashboard", className="mb-4"),
                dbc.Alert(
                    "Analytics dashboard coming soon...",
                    color="info"
                )
            ])
        ])

    def create_404_page():
        """Page 404"""
        return html.Div([
            dbc.Container([
                html.Div([
                    html.H1("404", className="display-1 text-muted"),
                    html.H3("Page Not Found"),
                    html.P("The page you're looking for doesn't exist."),
                    dbc.Button("Go Home", href="/", color="primary")
                ], className="text-center py-5")
            ])
        ])