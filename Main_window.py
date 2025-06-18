# Mainwindow
from dash import html, dcc
import dash_bootstrap_components as dbc

def get_main_layout():
    return [
        html.H2("Main Interface", className="text-success my-4"),

        dbc.Label("Please enter a substance:", className="mb-2"),

        dcc.Input(
            id="substance-input",
            type="text",
            placeholder="e.g., acetone",
            className="form-control mb-4"
        ),

        html.Div("Draw your structure here", className="fw-bold mb-2 mt-4"),

        html.Div(style={
            "border": "2px dashed #aaa",
            "height": "250px",
            "width": "250px",
            "backgroundColor": "#f9f9f9"
        })
    ]