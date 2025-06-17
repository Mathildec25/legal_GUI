from dash import Dash, html

app = Dash(__name__)

app.layout = html.Div([
    html.H1("interface GUI"),
    html.Div("Bienvenue dans votre interface Dash.")
])

if __name__ == "__main__":
    app.run(debug=True)
# premier test