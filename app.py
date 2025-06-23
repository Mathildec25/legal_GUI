from dash import Dash, html
import dash_bootstrap_components as dbc
from layout.homepage_layout import get_home_layout
from layout.main_layout import get_main_layout
from callbacks.navigation import register_navigation_callbacks
from callbacks.main_callbacks import register_main_callbacks

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

# Enregistrement des callbacks
register_navigation_callbacks(app)
register_main_callbacks(app)

if __name__ == "__main__":
    app.run(debug=True)
