from dash import Dash
import dash_bootstrap_components as dbc
# Import du layout principal de l'application
from layout.main_layout import get_main_layout
# Import des fonctions de callbacks
from callbacks.substance_callbacks import register_callbacks
# Chargement des données à partir du fichier utils
from utils.data_loader import df

app = Dash(
    __name__,
    title="SafeLab",  # Le titre qui apparaîtra dans l’onglet
    external_stylesheets=[dbc.themes.FLATLY],# Thème Bootstrap utilisé
    suppress_callback_exceptions=True # Autorise la gestion de pages dynamiques (multi-pages)
)



# Enregistrement des callbacks
register_callbacks(app,df)

# Définition du layout principal (structure visuelle de l'app)
app.layout = get_main_layout()
# Lancement de l'application 
if __name__ == "__main__":
    app.run(debug=True)
