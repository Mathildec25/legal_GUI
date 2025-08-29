"""
SafeLab 2.0 - Regulatory Compliance Platform
Application principale avec interface moderne
"""

import dash
from dash import html
import dash_bootstrap_components as dbc
import pandas as pd

# Imports locaux
from layout.main_layout import get_main_layout
from callbacks.substance_callbacks import register_callbacks
from callbacks.database_callbacks import register_database_callbacks
from utils.data_loader import df

# Configuration de l'application
app = dash.Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.BOOTSTRAP,
        dbc.icons.FONT_AWESOME,
        "/assets/modern_styles.css"  # Nos styles personnalisés
    ],
    suppress_callback_exceptions=True,
    title="SafeLab - Regulatory Compliance Platform",
    update_title="SafeLab - Loading...",
    meta_tags=[
        {
            "name": "viewport",
            "content": "width=device-width, initial-scale=1.0"
        },
        {
            "name": "description",
            "content": "Advanced regulatory compliance platform for controlled substances"
        }
    ]
)

# Configuration du serveur
server = app.server

# Layout principal de l'application
app.layout = html.Div([
    
    # Contenu principal
    get_main_layout()
], className="app-wrapper")

# Enregistrement des callbacks
register_callbacks(app, df)
register_database_callbacks(app, df)

# JavaScript côté client pour les fonctionnalités avancées
app.clientside_callback(
    """
    function() {
        // Animation d'entrée pour les éléments
        window.addEventListener('DOMContentLoaded', function() {
            // Ajouter les classes d'animation
            const elements = document.querySelectorAll('.fade-in, .slide-up');
            elements.forEach((el, index) => {
                el.style.animationDelay = `${index * 0.1}s`;
            });
            
            // Smooth scrolling
            document.querySelectorAll('a[href^="#"]').forEach(anchor => {
                anchor.addEventListener('click', function (e) {
                    e.preventDefault();
                    const target = document.querySelector(this.getAttribute('href'));
                    if (target) {
                        target.scrollIntoView({
                            behavior: 'smooth',
                            block: 'start'
                        });
                    }
                });
            });
        });
        
        return window.dash_clientside.no_update;
    }
    """,
    dash.dependencies.Output('url', 'search'),
    dash.dependencies.Input('url', 'pathname')
)

if __name__ == "__main__":
    # Configuration pour le développement
    app.run(
        debug=True,
        host="127.0.0.1",   # uniquement accessible en local
        port=8050,
        dev_tools_hot_reload=True,
        dev_tools_ui=True
    )
