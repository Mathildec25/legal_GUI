from dash import Dash, html, Input, Output, dcc

app = Dash(__name__)

app.layout = html.Div([
    html.Div(id="homepage", children=[
        html.H1("REGULATORY ANALYSIS TOOL", style={"fontSize": "32px"}),

        html.P(
            "This graphical interface allows you to determine the legal status of organic compounds "
            "and evaluate the safest synthetic pathway based on regulatory and safety data.",
            style={
                "fontSize": "16px",
                "lineHeight": "1.6",
                "maxWidth": "500px",
                "marginBottom": "30px",
                "marginTop": "20px"
            }
        ),

        html.Button("Start", id="start-button", n_clicks=0,
                    style={"padding": "10px 20px", "fontSize": "16px"})
    ], style={
        "maxWidth": "600px",
        "margin": "100px auto 0 auto",
        "textAlign": "left",
        "padding": "40px"
    }),

    html.Div(id="main-interface", style={"display": "none"})
])

@app.callback(
    Output("homepage", "style"),
    Output("main-interface", "style"),
    Output("main-interface", "children"),
    Input("start-button", "n_clicks")
)
def show_main_interface(n_clicks):
    if n_clicks and n_clicks > 0:
        return {"display": "none"}, {"display": "block"}, [

            html.H1("Main Interface", style={"fontSize": "28px"}),

            html.Label("Please enter a substance", style={
                "fontSize": "16px",
                "marginTop": "20px"
            }),

            dcc.Input(
                id="substance-input",
                type="text",
                placeholder="e.g., acetone",
                style={
                    "width": "100%",
                    "padding": "10px",
                    "marginTop": "10px",
                    "fontSize": "16px"
                }
            ),

            html.Div("Draw your structure here", style={
                "fontWeight": "bold",
                "marginTop": "40px",
                "marginBottom": "10px"
            }),

            html.Div(style={
                "border": "2px dashed #aaa",
                "height": "250px",
                 "width":"250px",
                "backgroundColor": "#f9f9f9"
            })
        ]
    return {"display": "block"}, {"display": "none"}, []

if __name__ == "__main__":
    app.run(debug=True)
