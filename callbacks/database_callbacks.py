from dash import Input, Output, State, callback_context
import dash_bootstrap_components as dbc
from dash import html
import pandas as pd
from layout.database_page import create_modern_datatable

def register_database_callbacks(app, df):
    """
    Callbacks spécifiques à la page Database
    """
    
    # Callback principal: gestion de la table avec filtres
    @app.callback(
        Output("database-table-container", "children"),
        [Input("db-quick-search", "value"),
         Input("db-category-filter", "value"),
         Input("db-view-options", "value"),
         Input("refresh-btn", "n_clicks")],
        prevent_initial_call=False
    )
    def update_database_table(search_value, category_filter, view_options, refresh_clicks):
        """
        Met à jour la table selon les filtres appliqués
        """
        df_filtered = df.copy()
        
        # Filtre par recherche rapide
        if search_value:
            mask = df_filtered.apply(
                lambda row: row.astype(str).str.contains(
                    search_value, case=False, na=False
                ).any(), axis=1
            )
            df_filtered = df_filtered[mask]
        
        # Filtre par catégorie réglementaire
        if category_filter and category_filter != "all":
            df_filtered = apply_category_filter(df_filtered, category_filter)
        
        # Options d'affichage
        view_options = view_options or []
        
        # Créer la table moderne
        return create_modern_datatable(df_filtered, view_options)
    
    def apply_category_filter(df_filtered, category):
        """Applique les filtres par catégorie"""
        if category == "belgian":
            return df_filtered[df_filtered['belgian_annex'].notna()]
        elif category == "eu_273":
            return df_filtered[df_filtered['eu_annex/category_273/2004'].notna()]
        elif category == "eu_111":
            return df_filtered[df_filtered['eu_annex/category_111/2005'].notna()]
        elif category == "un":
            return df_filtered[df_filtered['uno'].notna()]
        elif category == "none":
            return df_filtered[
                (df_filtered['belgian_annex'].isna()) &
                (df_filtered['eu_annex/category_273/2004'].isna()) &
                (df_filtered['eu_annex/category_111/2005'].isna())
            ]
        return df_filtered
    
    # Callback: Export CSV
    @app.callback(
        Output("export-csv-btn", "children"),
        Input("export-csv-btn", "n_clicks"),
        State("database-table-modern", "derived_virtual_data"),
        prevent_initial_call=True
    )
    def export_to_csv(n_clicks, table_data):
        """Export des données filtrées en CSV"""
        if n_clicks and table_data:
            df_export = pd.DataFrame(table_data)
            # Ici vous pouvez implémenter l'export réel
            return [
                html.I(className="fas fa-check me-2"),
                "Exported!"
            ]
        return [
            html.I(className="fas fa-download me-2"),
            "Export CSV"
        ]
    
    # Callback: Sélection de lignes
    @app.callback(
        Output("export-controls", "children"),
        Input("database-table-modern", "selected_rows"),
        prevent_initial_call=True
    )
    def update_selection_info(selected_rows):
        """Met à jour les informations de sélection"""
        if selected_rows:
            return dbc.Alert(
                f"{len(selected_rows)} row(s) selected",
                color="info",
                size="sm",
                className="mb-0"
            )
        return ""