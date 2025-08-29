from dash import dcc, html
import dash_bootstrap_components as dbc
from layout.search_page import get_search_page

def get_main_layout():
    """
    Layout principal repensé avec design moderne
    Header horizontal + contenu pleine largeur
    """
    
    # Header moderne horizontal
    header = html.Header([
        dbc.Container([
            dbc.Row([
                # Logo section
                dbc.Col([
                    html.A([
                        html.I(className="fas fa-atom me-2"),
                        html.Span("SafeLab", className="logo-text")
                    ], href="/", className="logo-link")
                ], width="auto"),
                
                # Navigation pills
                dbc.Col([
                    html.Nav([
                        dbc.ButtonGroup([
                            dbc.Button([
                                html.I(className="fas fa-search me-2"),
                                "Search"
                            ], 
                            id="nav-search",
                            href="/search",
                            color="primary",
                            outline=True,
                            className="nav-pill"),
                            
                            dbc.Button([
                                html.I(className="fas fa-database me-2"),
                                "Database"
                            ],
                            id="nav-database", 
                            href="/database",
                            color="primary",
                            outline=True,
                            className="nav-pill"),
                            
                            dbc.Button([
                                html.I(className="fas fa-chart-bar me-2"),
                                "Analytics"
                            ],
                            id="nav-analytics",
                            href="/analytics", 
                            color="primary",
                            outline=True,
                            className="nav-pill disabled"),
                        ], className="nav-pills")
                    ])
                ], className="d-flex justify-content-center"),
                
                # User menu
                dbc.Col([
                    html.Div([
                        dbc.Button([
                            html.I(className="fas fa-user-circle")
                        ], 
                        color="light",
                        className="user-avatar-btn",
                        id="user-menu-btn"),
                        
                        dbc.Tooltip(
                            "User Menu",
                            target="user-menu-btn"
                        )
                    ])
                ], width="auto", className="d-flex justify-content-end")
            ], className="align-items-center")
        ], fluid=True, className="header-container")
    ], className="modern-header")
    
    # Hero section (visible uniquement sur la page d'accueil)
    hero_section = html.Section([
        html.Div([
            html.H1([
                "Regulatory Compliance Platform",
                html.Div(className="hero-decoration")
            ], className="hero-title"),
            html.P(
                "Search and analyze controlled substances with advanced molecular visualization and real-time regulatory information",
                className="hero-subtitle"
            )
        ], className="hero-content")
    ], className="hero-section", id="hero-section")
    
    # Main content wrapper
    main_content = html.Main([
        dcc.Location(id="url", refresh=False),
        html.Div(id="page-content", className="page-content-wrapper")
    ], className="main-content")
    
    # Footer moderne
    footer = html.Footer([
        dbc.Container([
            dbc.Row([
                dbc.Col([
                    html.P([
                        "© 2024 SafeLab - Regulatory Compliance Platform",
                        html.Br(),
                        html.Small("Built for regulatory compliance and molecular analysis", 
                                  className="text-muted")
                    ])
                ], md=8),
                dbc.Col([
                    html.Div([
                        html.A([html.I(className="fab fa-github")], 
                              href="#", className="footer-link me-3"),
                        html.A([html.I(className="fas fa-book")], 
                              href="#", className="footer-link me-3"),
                        html.A([html.I(className="fas fa-envelope")], 
                              href="#", className="footer-link")
                    ], className="footer-links")
                ], md=4, className="text-end")
            ])
        ], fluid=True)
    ], className="modern-footer")
    
    return html.Div([
        header,
        hero_section,
        main_content,
        footer
    ], className="app-container")