from dash import Input, Output

def register_navigation_callbacks(app):
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
        return {"display": "block"}, {"display": "none"}
