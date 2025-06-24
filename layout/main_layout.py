from dash import dcc, html,dash_table
import dash_bootstrap_components as dbc
from utils.data_loader import df

def get_main_layout():
    
    # Barre de recherche centrée, rapprochée et décalée vers la gauche
    search_bar = dbc.Row(
        [
    dbc.Col(
                dbc.Button("Draw", id="draw-button", n_clicks=0,style={"backgroundColor": "#107c7e", "borderColor": "#0e7c7e", "color": "white"}),
                width="auto"
            ),
            dbc.Col(
                dcc.Dropdown(
                    id='input-substance',
                    options=[{"label": nom, "value": nom} for nom in sorted(df['substance_name'].unique())],
                    placeholder="Search or select a substance...",
                    searchable=True,
                    style={"width": "100%"}
                ),
                width=8
            ),
            dbc.Col(
                dbc.Button("Search", id="search-button", n_clicks=0,style={"backgroundColor": "#107c7e","color": "white"})
               ,width="auto"
            )
        ],
        className="align-items-center",
        style={
            "marginTop": "120px",
            "marginLeft": "60px",   # Décalage vers la gauche
            "gap": "2.5px"           # Proximité entre les composants
        }
    )

    # Barre latérale (menu vertical)
    sidebar = html.Div(
    children=[
        html.Div(
            html.Img(src="/assets/ChemCheck.png", style={"height": "70px", "width": "auto"}),
            style={
                "display": "flex",
                "justifyContent": "center",
                "alignItems": "center",
                "paddingTop": "20px",
                "paddingBottom": "20px"
            }
        ),
        dbc.Nav(
            children=[
                dbc.NavItem(dbc.NavLink("DataBase", href="#", active=True,id="open-database",style={"color": "white"})),
                dbc.NavItem(dbc.NavLink("Help", href="#", active=True, style={"color": "white"}))
            ],
            vertical=True,
            pills=True,
        )
    ],
    className="sidebar",  
    style={"height": "100vh"}
)

# Modale contenant la DataTable
    database_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Full Database")),
        dbc.ModalBody(
            dash_table.DataTable(
                id="database-table",
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                page_size=10,
                style_table={"overflowX": "auto", "padding": "20px"},
                style_cell={"textAlign": "left", "padding": "5px"},
                style_header={
                    "backgroundColor": "#0e7c7e",
                    "color": "white",
                    "fontWeight": "bold"
                }
            )
        ),
        dbc.ModalFooter(
            dbc.Button("Close", id="close-database", className="ms-auto", n_clicks=0, color="secondary")
        ),
    ],
    id="modal-database",
    is_open=False,
    size="xl",
)


    # Zone modale pour dessiner
    draw_modal = dbc.Modal(
        children=[
            dbc.ModalHeader(dbc.ModalTitle("Draw Structure")),
            dbc.ModalBody(
                html.Div("Drawing area", style={
                    "border": "2px dashed #aaa",
                    "height": "200px",
                    "width": "100%",
                    "backgroundColor": "#F0F2F5"
                })
            ),
            dbc.ModalFooter(
                dbc.Button("Close", id="close-draw", className="ms-auto", n_clicks=0, color="success")
            )
        ],
        id="modal-draw",
        is_open=False,
        size="lg"
    )

    # Composition finale
    return dbc.Row([
        dbc.Col(sidebar, width=2),
        dbc.Col([
            search_bar,
            html.Div(id="substance-details", className="mt-4"),
            draw_modal,
           database_modal 
        ], width=10)
    ])
