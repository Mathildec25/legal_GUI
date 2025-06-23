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
            html.Div(id="main-interface", children=get_main_layout(), style={"display": "none"})
        ], width=8, className="offset-md-2")
    ])
], fluid=True)

@app.callback(
    Output("homepage", "style"),
    Output("main-interface", "style"),
    Input("start-button", "n_clicks"),
    Input("back-button", "n_clicks")
)
def switch_home_and_main(start_clicks, back_clicks):
    if start_clicks and (not back_clicks or start_clicks > back_clicks):
        return {"display": "none"}, {"display": "block"}
    elif back_clicks and (not start_clicks or back_clicks >= start_clicks):
        return {"display": "block"}, {"display": "none"}
#Show by default the homePage
    return {"display": "block"}, {"display": "none"}

if __name__ == "__main__":
    app.run(debug=True)

