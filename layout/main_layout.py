
from dash import dcc, html
import dash_bootstrap_components as dbc
from layout.search_page import get_search_page

# Main layout function
# ===========================
# This function returns the full structure of the Dash application
# with a sidebar on the left and dynamic content on the right.
def get_main_layout():


    #
    #  Sidebar with navigation
    # ===========================

    sidebar = html.Div([
        # Application logo (icon)
        html.Div(
            html.Img(src="assets/favicon.ico", className="logo-img"),
            className="logo-container"
        ),

        # Application title
        html.Div("SafeLab", className="logo-text"),

        html.Hr(),  # Visual separator

        #  Vertical navigation (Font Awesome icons + text)
        dbc.Nav([
            # Link to the search page
            dbc.NavLink([
                html.I(className="fas fa-search nav-icon"),
                html.Span("Search", className="nav-text")
            ],
            href="/search",
            active="exact",
            className="navlink"),

            # Link to the database page
            dbc.NavLink([
                html.I(className="fas fa-database", style={"marginRight": "10px"}),
                html.Span("Database", className="nav-text")
            ],
            href="/database",
            active="exact",
            className="navlink")
        ],
        vertical=True,
        pills=True)
    ],
    className="sidebar sidebar-open")

    # Main page structure (2 columns: sidebar + content)
    # ===========================


    return dbc.Row([
        # Left column: sidebar (2/12)
        dbc.Col(sidebar, width=2),

        # Right column: dynamic content (10/12)
        dbc.Col([
            dcc.Location(id="url", refresh=False),              # Handles the URL for navigation
            get_search_page(hidden=True),                       # Loads the Search page, hidden by default
            html.Div(id="page-content")                         # Displays the content dynamically based on URL
        ], width=10)
    ],
    className="g-0") # "g-0" removes horizontal spacing between Bootstrap columns
