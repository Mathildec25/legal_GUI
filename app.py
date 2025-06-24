# app.py

from dash import Dash
import dash_bootstrap_components as dbc
from layout.main_layout import get_main_layout
from callbacks.substance_callbacks import register_callbacks
from utils.data_loader import df

# Instanciation de l'app
app = Dash(__name__, external_stylesheets=[dbc.themes.CERULEAN])

# Enregistrement des callbacks
register_callbacks(app, df)

# Layout principal directement affiché
app.layout = get_main_layout()

if __name__ == "__main__":
    app.run(debug=True)
