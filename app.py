from dash import Dash
import dash_bootstrap_components as dbc
# Import the main layout of the application
from layout.main_layout import get_main_layout
# Import the callback registration functions
from callbacks.substance_callbacks import register_callbacks
# Load data from the utils module
from utils.data_loader import df

app = Dash(
    __name__,
    title="SafeLab",  # Title that will appear in the browser tab
    external_stylesheets=[dbc.themes.FLATLY],# Bootstrap theme used
    suppress_callback_exceptions=True # Allows handling of dynamic (multi-page) content
)



# Register the callbacks
register_callbacks(app,df)

# Set the main layout (visual structure of the app)
app.layout = get_main_layout()
# Launch the application
if __name__ == "__main__":
    app.run(debug=True)
