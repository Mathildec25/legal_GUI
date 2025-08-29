from dash import html, dash_table, dcc
import dash_bootstrap_components as dbc
from utils.data_loader import df
import pandas as pd

def get_database_page():
    """
    Page Database repensée avec design moderne et fonctionnalités avancées
    """
    
    # Statistiques de la base de données
    # Statistiques de la base de données
    total_substances = len(df)
    regulated_substances = len(df[df['Belgian annex'].notna()])
    controlled_273 = len(df[df['EU_Annex_Category_273_2004'].notna()])
    controlled_111 = len(df[df['EU_Annex_Category_111_2005'].notna()])

    return html.Div([
        
        # Header avec statistiques
        html.Section([
            dbc.Container([
                html.H2([
                    html.I(className="fas fa-database me-3 text-primary"),
                    "Regulatory Database"
                ], className="database-title mb-4"),
                
                # Cards de statistiques
                dbc.Row([
                    dbc.Col([
                        create_stat_card(
                            "Total Substances",
                            total_substances,
                            "fas fa-flask",
                            "primary"
                        )
                    ], md=3),
                    dbc.Col([
                        create_stat_card(
                            "Belgian Regulated",
                            regulated_substances,
                            "fas fa-flag",
                            "warning"
                        )
                    ], md=3),
                    dbc.Col([
                        create_stat_card(
                            "EU 273/2004",
                            controlled_273,
                            "fas fa-shield-alt",
                            "danger"
                        )
                    ], md=3),
                    dbc.Col([
                        create_stat_card(
                            "EU 111/2005",
                            controlled_111,
                            "fas fa-globe-europe",
                            "info"
                        )
                    ], md=3)
                ], className="stats-row mb-4")
            ])
        ], className="database-header-section"),
        
        # Contrôles et filtres
        html.Section([
            dbc.Container([
                dbc.Card([
                    dbc.CardBody([
                        dbc.Row([
                            # Barre de recherche rapide
                            dbc.Col([
                                html.Label("Quick Search", className="form-label fw-bold"),
                                dbc.InputGroup([
                                    dbc.Input(
                                        id="db-quick-search",
                                        placeholder="Search in all columns...",
                                        className="database-search-input"
                                    ),
                                    dbc.Button([
                                        html.I(className="fas fa-search")
                                    ], id="db-search-btn", color="primary")
                                ])
                            ], md=4),
                            
                            # Filtre par catégorie
                            dbc.Col([
                                html.Label("Regulatory Category", className="form-label fw-bold"),
                                dcc.Dropdown(
                                    id="db-category-filter",
                                    options=[
                                        {"label": "All Categories", "value": "all"},
                                        {"label": "Belgian Annexes", "value": "belgian"},
                                        {"label": "EU 273/2004", "value": "eu_273"},
                                        {"label": "EU 111/2005", "value": "eu_111"},
                                        {"label": "UN Conventions", "value": "un"},
                                        {"label": "Non-Regulated", "value": "none"}
                                    ],
                                    value="all",
                                    className="database-dropdown"
                                )
                            ], md=4),
                            
                            # Options d'affichage
                            dbc.Col([
                                html.Label("View Options", className="form-label fw-bold"),
                                html.Div([
                                    dbc.Checklist(
                                        id="db-view-options",
                                        options=[
                                            {"label": "Show SMILES", "value": "smiles"},
                                            {"label": "Compact View", "value": "compact"}
                                        ],
                                        value=["smiles"],
                                        className="database-checklist"
                                    )
                                ])
                            ], md=4)
                        ])
                    ])
                ], className="controls-card mb-4")
            ])
        ], className="database-controls-section"),
        
        # Table de données moderne
        html.Section([
            dbc.Container([
                dbc.Card([
                    dbc.CardHeader([
                        dbc.Row([
                            dbc.Col([
                                html.H5([
                                    html.I(className="fas fa-table me-2"),
                                    "Database Contents"
                                ], className="mb-0")
                            ], md=6),
                            dbc.Col([
                                html.Div([
                                    dbc.ButtonGroup([
                                        dbc.Button([
                                            html.I(className="fas fa-download me-2"),
                                            "Export CSV"
                                        ], id="export-csv-btn", size="sm", color="success", outline=True),
                                        
                                        dbc.Button([
                                            html.I(className="fas fa-file-excel me-2"),
                                            "Export Excel"
                                        ], id="export-excel-btn", size="sm", color="info", outline=True),
                                        
                                        dbc.Button([
                                            html.I(className="fas fa-sync me-2"),
                                            "Refresh"
                                        ], id="refresh-btn", size="sm", color="secondary", outline=True)
                                    ])
                                ], className="export-controls")
                            ], md=6, className="text-end")
                        ])
                    ], className="database-card-header"),
                    
                    dbc.CardBody([
                        # Loading indicator
                        dbc.Spinner(
                            html.Div(id="database-table-container"),
                            color="primary",
                            type="border",
                            size="sm"
                        )
                    ], className="p-0")
                ], className="database-main-card")
            ])
        ], className="database-table-section")
        
    ], className="database-page-container")

def create_stat_card(title, value, icon, color):
    """Crée une carte de statistique"""
    return dbc.Card([
        dbc.CardBody([
            html.Div([
                html.Div([
                    html.I(className=f"{icon} stat-icon"),
                    html.H3(f"{value:,}", className="stat-value"),
                    html.P(title, className="stat-label")
                ], className="stat-content")
            ])
        ])
    ], className=f"stat-card stat-card-{color}")

def create_modern_datatable(df_filtered, view_options):
    """
    Crée une DataTable moderne avec styling avancé
    """
    
    # Configuration des colonnes selon les options
    columns_config = get_columns_config(view_options)
    
    return dash_table.DataTable(
        id="database-table-modern",
        columns=columns_config,
        data=df_filtered.to_dict("records"),
        
        # Configuration générale
        editable=False,
        filter_action="native",
        sort_action="native",
        page_action="native",
        page_current=0,
        page_size=25,

        
        # Sélection de lignes
        row_selectable="multi",
        selected_rows=[],
        
        # Style de la table
        style_table={
            "overflowX": "auto",
            "minHeight": "400px"
        },
        
        # Style des cellules
        style_cell={
            "textAlign": "left",
            "padding": "12px 16px",
            "fontFamily": "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
            "fontSize": "14px",
            "border": "1px solid #e2e8f0",
            "whiteSpace": "normal",
            "height": "auto",
            "minWidth": "120px",
            "maxWidth": "300px"
        },
        
        # Style des headers
        style_header={
            "backgroundColor": "#f8fafc",
            "color": "#374151",
            "fontWeight": "600",
            "border": "1px solid #d1d5db",
            "textAlign": "center"
        },
        
        # Style des données
        style_data={
            "backgroundColor": "white",
            "color": "#374151",
        },
        
        # Style conditionnel
        style_data_conditional=get_conditional_styles(),
        
        # Configuration avancée
        virtualization=True,
        fixed_rows={"headers": True},
        
        # Export
        export_format="xlsx",
        export_headers="display",
        
        # CSS classes personnalisées
        css=[{
            "selector": ".dash-table-tooltip",
            "rule": "background-color: #1f2937; color: white; border-radius: 6px;"
        }],
        
        # Tooltip
        tooltip_data=create_tooltips(df_filtered),
        tooltip_delay=0,
        tooltip_duration=None
    )

def get_columns_config(view_options):
    """Configuration des colonnes selon les options d'affichage"""
    
    base_columns = [
        {
            "name": "Name", 
            "id": "NAME", 
            "type": "text",
            "presentation": "markdown"
        },
        {
            "name": "CAS Number", 
            "id": "CAS", 
            "type": "text"
        },
        {
            "name": "IUPAC Name", 
            "id": "IUPAC", 
            "type": "text"
        },
        {
            "name": "Belgian Annex", 
            "id": "Belgian Annex", 
            "type": "text"
        },
        {
            "name": "EU 273/2004", 
            "id": "EU_Annex_Category_273_2004", 
            "type": "text"
        },
        {
            "name": "EU 111/2005", 
            "id": "EU_Annex_Category_111_2005", 
            "type": "text"
        },
        {
            "name": "UN Convention", 
            "id": "UNO", 
            "type": "text"
        }
    ]
    
    # Ajouter SMILES si demandé
    if "smiles" in view_options:
        base_columns.extend([
            {
                "name": "SMILES", 
                "id": "smiles", 
                "type": "text",
                "presentation": "markdown"
            },
            {
                "name": "Canonical SMILES", 
                "id": "canonical_smiles", 
                "type": "text",
                "presentation": "markdown"
            }
        ])
    
    return base_columns

def get_conditional_styles():
    """Styles conditionnels pour la table"""
    return [
        # Lignes alternées
        {
            "if": {"row_index": "odd"},
            "backgroundColor": "#f9fafb"
        },
        
        # Ligne sélectionnée
        {
            "if": {"state": "selected"},
            "backgroundColor": "#dbeafe",
            "border": "1px solid #3b82f6"
        },
        
        # Ligne survolée
        {
            "if": {"state": "active"},
            "backgroundColor": "#f0f9ff",
            "border": "1px solid #0ea5e9"
        },
        
        # Cellules Belgian Annex
        {
            "if": {
                "filter_query": "{belgian_annex} contains annex",
                "column_id": "belgian_annex"
            },
            "backgroundColor": "#fef3c7",
            "color": "#92400e"
        },
        
        # Cellules EU 273/2004 - Category 1
        {
            "if": {
                "filter_query": "{eu_annex/category_273/2004} contains category 1",
                "column_id": "eu_annex/category_273/2004"
            },
            "backgroundColor": "#fee2e2",
            "color": "#dc2626"
        },
        
        # Cellules EU 273/2004 - Category 2
        {
            "if": {
                "filter_query": "{eu_annex/category_273/2004} contains category 2",
                "column_id": "eu_annex/category_273/2004"
            },
            "backgroundColor": "#fed7aa",
            "color": "#ea580c"
        },
        
        # Cellules EU 111/2005
        {
            "if": {
                "filter_query": "{eu_annex/category_111/2005} contains category",
                "column_id": "eu_annex/category_111/2005"
            },
            "backgroundColor": "#dbeafe",
            "color": "#2563eb"
        },
        
        # SMILES columns - police monospace
        {
            "if": {"column_id": ["smiles", "canonical_smiles"]},
            "fontFamily": "Monaco, Consolas, 'Courier New', monospace",
            "fontSize": "12px"
        }
    ]

def create_tooltips(df_filtered):
    """Crée les tooltips pour la table"""
    tooltip_data = []
    
    for index, row in df_filtered.iterrows():
        tooltip_row = {}
        
        # Tooltip pour le nom
        if pd.notna(row.get('iupac')):
            tooltip_row['name'] = {
                'value': f"IUPAC: {row['iupac']}",
                'type': 'markdown'
            }
        
        # Tooltip pour CAS
        tooltip_row['CAS'] = {
            'value': f"CAS Registry Number: {row['CAS']}",
            'type': 'text'
        }
        
        # Tooltip pour les annexes belges
        if pd.notna(row.get('belgian_annex')):
            tooltip_row['belgian_annex'] = {
                'value': f"Belgian regulation: {row['belgian_annex']}",
                'type': 'text'
            }
        
        tooltip_data.append(tooltip_row)
    
    return tooltip_data