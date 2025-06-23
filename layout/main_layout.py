import dash_bootstrap_components as dbc
from dash import dcc, dash_table
import pandas as pd

# Chargement des données
df = pd.read_csv("data/donnees.csv", sep=";", encoding="utf-8")
df.columns = [col.lower().replace(" ", "").strip() for col in df.columns]

def get_main_layout():
    return dbc.Container([
        dbc.Row([
            dbc.Col(
                dbc.Card([
                    dbc.CardHeader("Main Interface", className="text-success text-center"),
                    dbc.CardBody([
                        dbc.Label("Please enter a substance:", className="fw-bold mb-2"),
                        dcc.Input(
                            id="input-substance",
                            type="text",
                            placeholder="e.g., acetone",
                            debounce=True,
                            className="form-control mb-4"
                        ),

                        dbc.Spinner(dbc.Container(id="substance-details", className="mt-3")),

                        dbc.Container("Draw your structure here", className="h5 text-primary mt-4"),
                        dbc.Container(id="drawing-box", className="border border-2 border-secondary rounded w-100 h-100 mt-2", style={"height": "200px"}),

                        dbc.Button("Back to Homepage", id="back-button", color="secondary", className="mt-4"),
                        html_hr(),

                        dbc.DropdownMenu(
                            label="Menu",
                            children=[
                                dbc.DropdownMenuItem("View the database", id="menu-show-table")
                            ],
                            color="info",
                            className="mb-3"
                        ),

                        dbc.Collapse(
                            dash_table.DataTable(
                                id='table-database',
                                columns=[{"name": col, "id": col} for col in df.columns],
                                data=df.to_dict('records'),
                                page_size=10,
                                style_table={
                                    'maxHeight': '300px',
                                    'overflowY': 'auto',
                                    'overflowX': 'auto'
                                },
                                style_cell={
                                    'textAlign': 'left',
                                    'padding': '5px'
                                },
                                style_header={
                                    'backgroundColor': '#f8f9fa',
                                    'fontWeight': 'bold'
                                }
                            ),
                            id="collapse-table",
                            is_open=False
                        )
                    ])
                ]),
            width=12)
        ])
    ], fluid=True)

def html_hr():
    """Remplace html.Hr() par un séparateur visuel Bootstrap-compatible"""
    return dbc.Container(style={"borderTop": "1px solid #dee2e6", "marginTop": "2rem", "marginBottom": "2rem"})
