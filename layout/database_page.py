from dash import html, dash_table
from utils.data_loader import df

# Fonction qui génère la page "Database"
def get_database_page():
    return html.Div([

        # Titre de la page
        html.H2("Full Database", style={
            "marginTop": "20px",
            "marginBottom": "20px",
            "marginRight": "100px",
            "textAlign": "center",
            "fontSize": "36px",
            "fontWeight": "bold",
            "color": "#003366",
            "fontFamily": "Arial, sans-serif"
        }),

        # Table de données interactive
        html.Div([
            dash_table.DataTable(
                id="database-table",
                columns=[{"name": col, "id": col} for col in df.columns],
                data=df.to_dict("records"),

                # Table non modifiable par l'utilisateur
                editable=False,

                # Filtres et tri natifs intégrés
                filter_action="native",
                sort_action="native",

                #  Pagination désactivée pour scroll fluide
                page_action="none",
                fixed_rows={"headers": True},  # Header fixe au scroll

                # Style global du conteneur de la table
                style_table={
                    "width": "calc(100% - 80px)",       # Ajusté pour ne pas déborder de la sidebar
                    "height": "calc(100vh - 160px)",    # Pleine hauteur visible
                    "overflowX": "auto",
                    "overflowY": "auto",
                    "marginLeft": "0px",
                    "marginRight": "20px",
                    "border": "1px solid #ddd",
                    "boxShadow": "0 1px 3px rgba(0, 0, 0, 0.1)"
                },

                # Style des cellules (contenu)
                style_cell={
                    "textAlign": "left",
                    "padding": "10px 14px",
                    "minWidth": "120px",
                    "maxWidth": "300px",
                    "whiteSpace": "normal",
                    "border": "1px solid #f0f0f0",
                    "overflow": "hidden",
                    "textOverflow": "ellipsis",
                    "fontFamily": "Segoe UI, sans-serif",
                    "fontSize": "14px"
                },

                # Style de l’en-tête (sticky + jaune)
                style_header={
                    "backgroundColor": "#ffcb47",
                    "color": "black",
                    "fontWeight": "bold",
                    "position": "sticky",
                    "top": 0,
                    "zIndex": 1,
                    "border": "1px solid #ddd"
                },

                # Style des données
                style_data={
                    "backgroundColor": "white",
                    "color": "black"
                },

                # Mise en valeur : lignes impaires, ligne active
                style_data_conditional=[
                    {
                        "if": {"row_index": "odd"},
                        "backgroundColor": "#fafafa"
                    },
                    {
                        "if": {"state": "active"},
                        "backgroundColor": "#d7f0e2",
                        "border": "1px solid #999"
                    }
                ],

                # Performances + export CSV
                virtualization=True,
                export_format="csv"
            )
        ])
    ],
    style={
        "width": "100%",
        "overflow": "hidden",
        "height": "100vh"
    })
