# layout/main_layout.py

from dash import dcc, html
import dash_bootstrap_components as dbc

def get_main_layout():
    # Barre de recherche
    search_bar = dbc.Row(
        [
            dbc.Col(
                dbc.Button("Draw", color="success", className="ms-2", id="draw-button"),
                width="auto",
                style={"marginTop": "35px", "marginLeft": "30px"}
            ),
            dbc.Col(
                dcc.Input(
                    id="input-substance",
                    type="text",
                    placeholder="Enter substance name...",
                    debounce=True,
                    className="form-control"
                ),
                style={"marginTop": "35px", "marginLeft": "10px", "width": "90px"}
            ),
            dbc.Col(
                dbc.Button("Search", color="success", className="ms-2", id="search-button", n_clicks=0),
                width="auto",
                style={"marginTop": "35px", "marginLeft": "30px"}
            )
        ],
        className="g-0 flex-nowrap mt-3 mt-md-0",
        align="center"
    )

    # Sidebar (composants Dash uniquement)
    sidebar = dbc.Nav(
        children=[
            dbc.NavItem(
                dbc.NavLink("Menu", disabled=True, style={"color": "white", "fontSize": "36px"})
            ),
            dbc.NavItem(dbc.NavLink("DataBase", href="#", active=True,
                                     style={"color": "white", "backgroundColor": "#1E5631"})),
            dbc.NavItem(dbc.NavLink("Help", href="#", active=True,
                                     style={"color": "white", "backgroundColor": "#1E5631"}))
        ],
        vertical=True,
        pills=True,
        className="sidebar"
    )

    # Modal de dessin
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

    # Assemblage final
    return dbc.Row([
        dbc.Col(sidebar, width=2),
        dbc.Col([
            search_bar,
            html.Div(id="substance-details", className="mt-4"),
            draw_modal
        ], width=10)
    ])
