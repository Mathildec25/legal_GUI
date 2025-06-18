# Homepage
from dash import html, dcc
import dash_bootstrap_components as dbc

def get_home_layout():
    return html.Div(id="homepage", children=[
        html.H1("REGULATORY ANALYSIS TOOL", className="text-primary mb-4"),

        html.P(
            "This graphical interface allows you to determine the legal status of organic compounds "
            "and evaluate the safest synthetic pathway based on regulatory and safety data.",
            className="lead"
        ),

        dbc.Button("Start", id="start-button", n_clicks=0, color="primary", size="lg", className="mt-4")
    ], style={
        "marginTop": "100px",
        "padding": "40px"
    })
