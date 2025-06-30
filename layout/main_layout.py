
from dash import dcc, html
import dash_bootstrap_components as dbc

def get_main_layout():
    search_bar = dbc.Row([
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
    ], className="g-0 flex-nowrap mt-3 mt-md-0", align="center")

    sidebar = dbc.Nav([
        dbc.NavItem(dbc.NavLink("Menu", disabled=True, style={"color": "white", "fontSize": "36px"})),
        dbc.NavItem(dbc.NavLink("DataBase", href="#", active=True, style={"color": "white", "backgroundColor": "#1E5631"})),
        dbc.NavItem(dbc.NavLink("Help", href="#", active=True, style={"color": "white", "backgroundColor": "#1E5631"}))
    ], vertical=True, pills=True, className="sidebar")

    draw_modal = dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle("Draw Structure")),
        dbc.ModalBody(
            WebView(src="/assets/ketcher_embed.html", id="ketcher-frame", style={"height": "500px", "width": "100%"})
        ),
        dbc.ModalFooter([
            dbc.Button("Use structure", id="confirm-draw", color="primary", className="me-2"),
            dbc.Button("Close", id="close-draw", className="ms-auto", n_clicks=0, color="secondary")
        ])
    ], id="modal-draw", is_open=False, size="xl")

    return dbc.Row([
        dbc.Col(sidebar, width=2),
        dbc.Col([
            search_bar,
            html.Div(id="substance-details", className="mt-4"),
            draw_modal,
            dcc.Store(id="store-drawn-smiles")
        ], width=10)
    ])
