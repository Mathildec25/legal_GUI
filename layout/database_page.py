from dash import html, dash_table
from utils.data_loader import df

# Fonction qui retourne la page de la base de données
def get_database_page():
    return html.Div([
        html.H2("Full Database", style={"marginTop": "20px", "marginBottom": "20px"}),
    
    # Tableau interactif de la base de données
        dash_table.DataTable(
            id="database-table",
            columns=[
                {"name": col, "id": col}
                for col in df.columns
            ],
            data=df.to_dict("records"),
            editable=False,              
            filter_action="native",      # filtres 
            sort_action="native",        # tri
            page_action="native",        # pagination
            page_size=10,
            style_table={
                "overflowX": "auto",
                "overflowY": "scroll",
                "maxHeight": "600px"
            },
            style_cell={
                "textAlign": "left",
                "padding": "6px",
                "minWidth": "120px",
                "maxWidth": "300px",
                "whiteSpace": "normal"
            },
            style_header={
                "backgroundColor": "#ffcb47",
                 "color": "black",
                "fontWeight": "bold"
            },
            style_data_conditional=[
                {
                    "if": {"state": "active"},
                    "backgroundColor": "#d7f0e2",
                    "border": "1px solid #999"
                }
            ]
        )
    ], style={"padding": "20px"})
