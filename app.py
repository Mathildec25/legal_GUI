from dash import Dash, Input, Output, html
import dash_bootstrap_components as dbc

from Homepage import get_home_layout
from Main_window import get_main_layout

external_stylesheets = [dbc.themes.CERULEAN]
app = Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            get_home_layout(),
            html.Div(id="main-interface", style={"display": "none"})
        ], width=8, className="offset-md-2")
    ])
], fluid=True)

@app.callback(
    Output("homepage", "style"),
    Output("main-interface", "style"),
    Output("main-interface", "children"),
    Input("start-button", "n_clicks")
)
def show_main_interface(n_clicks):
    if n_clicks and n_clicks > 0:
        return {"display": "none"}, {"display": "block"}, get_main_layout()
    return {"display": "block"}, {"display": "none"}, []

if __name__ == "__main__":
    app.run(debug=True)

