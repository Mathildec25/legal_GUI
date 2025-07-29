from dash import html

def get_legal_procedure(reg):
    """
    Génère les blocs d'information juridique (Dash HTML) en fonction
    des régulations européennes (273/2004, 111/2005) et belges (annexes).

    Paramètre :
    - reg (dict) : dictionnaire contenant les informations de régulation pour une substance.

    Retour :
    - html.Div : ensemble de blocs HTML décrivant les obligations légales.
    """

    # Récupération et normalisation des catégories d'annexes
    cat_273       = str(reg.get('eu_annex/category_273/2004', '')).strip().lower()
    cat_111       = str(reg.get('eu_annex/category_111/2005', '')).strip().lower()
    belgian_annex = str(reg.get('belgian_annex', '')).strip().lower() # retourne nan

    blocks = []

    #  Fonction utilitaire pour créer des éléments <li> avec icône
    # Les icônes sont fournies par Font Awesome via la classe CSS "fas fa-..."
    # Exemple : fa-user-check, fa-database, etc.
    def li(icon, text):
        return html.Li([
            html.I(className=f"fas {icon} me-2 text-primary"),
            text
        ])

    # EU 273/2004 – Possession
    # ===============================

    if "annex i, category 1" in cat_273:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-shield-alt me-2 text-danger"),
                "EU 273/2004 – Category 1"
            ]),
            html.P([
                html.I(className="fas fa-file-alt me-2 text-secondary"),
                "Possession requires prior authorization. Client declaration is mandatory."
            ]),
            html.Ul([
                li("fa-certificate", "License required for possession."),
                li("fa-user-check", "Client declaration required for all transactions."),
                li("fa-exclamation-circle", "Declaration if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))

    elif "annex i,ii; subcategory 2a" in cat_273:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-folder-open me-2 text-warning"),
                "EU 273/2004 – SubCategory 2A"
            ]),
            html.P("Possession requires registration. Client declaration is mandatory."),
            html.Ul([
                li("fa-user-edit", "End users must register before possessing the substance."),
                li("fa-user-check", "Client declaration required."),
                li("fa-exclamation-circle", "Declaration if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))

    elif "annex i,ii; subcategory 2b" in cat_273:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-folder-minus me-2 text-warning"),
                "EU 273/2004 – SubCategory 2B"
            ]),
            html.P("Registration not required, but client declaration is mandatory if quantity exceeds threshold."),
            html.Ul([
                li("fa-ban", "No registration required for end users."),
                li("fa-user-check", "Client declaration required if threshold exceeded."),
                li("fa-exclamation-circle", "Declaration if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))

    elif "annex i, category 2" in cat_273:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-folder me-2 text-warning"),
                "EU 273/2004 – Category 2"
            ]),
            html.P("Possession requires registration. Client declaration is mandatory."),
            html.Ul([
                li("fa-id-badge", "Registration required."),
                li("fa-user-check", "Client declaration required."),
                li("fa-exclamation-circle", "Declaration if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))

    elif "annex i, category 3" in cat_273:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-eye me-2 text-muted"),
                "EU 273/2004 – Category 3"
            ]),
            html.P("No license or registration required. Vigilance is required."),
            html.Ul([
                li("fa-check-circle", "No authorization or registration needed."),
                li("fa-exclamation-triangle", "Client declaration not required."),
                li("fa-exclamation-circle", "Declaration if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))

    #  EU 111/2005 – Import/Export
    # ===============================

    if "annex, category 1" in cat_111:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-globe-europe me-2 text-danger"),
                "EU 111/2005 – Category 1"
            ]),
            html.P("Import/export authorization is required."),
            html.Ul([
                li("fa-arrow-down", "Import authorization is required for entering the EU."),
                li("fa-arrow-up", "Export authorization is required to leave the EU."),
                li("fa-user-check", "Client declaration is mandatory."),
                li("fa-exclamation-circle", "Declaration if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))

    elif "annex, category 2" in cat_111:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-exchange-alt me-2 text-warning"),
                "EU 111/2005 – Category 2"
            ]),
            html.P("Export authorization is required. Import is unrestricted."),
            html.Ul([
                li("fa-arrow-up", "Export authorization required."),
                li("fa-check", "No import authorization required."),
                li("fa-user-check", "Client declaration is mandatory."),
                li("fa-exclamation-circle", "Declaration required if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))

    elif "annex, category 3" in cat_111:
        blocks.append(html.Div([
            html.H5([
                html.I(className="fas fa-route me-2 text-muted"),
                "EU 111/2005 – Category 3"
            ]),
            html.P("Export authorization only required for certain countries."),
            html.Ul([
                li("fa-check", "No import authorization required."),
                li("fa-exclamation-triangle", "Export authorization required only for specific destinations."),
                li("fa-exclamation-triangle", "Client declaration not required."),
                li("fa-exclamation-circle", "Declaration if suspicious transaction.")
            ])
        ], style={"marginBottom": "20px"}))


    # Aucune régulation détectée
    # ===============================

    if not blocks:
        return html.Div([
            html.P([
                html.I(className="fas fa-info-circle me-2 text-muted"),
                "No regulatory procedure detected under EU 273/2004 or EU 111/2005."
            ]),
            html.P([
                html.I(className="fas fa-triangle-exclamation me-2 text-warning"),
                "Please verify manually with the competent authority if you're planning to handle or transport this substance."
            ])
        ])

    # Liens utiles (affichés dans tous les cas)
    # ===============================

    shared_links = html.Div([
        html.H5([
            html.I(className="fas fa-link me-2 text-primary"),
            "Useful Links"
        ]),
        html.Ul([
            html.Li(html.A(
                "Apply for license or registration",
                href="https://www.afmps.be/fr/humain/produits_particuliers/subst_specialement_reglementees/precurseurs/agrement_enregistrement",
                target="_blank"
            )),
            html.Li(html.A(
                "Client declaration",
                href="https://www.afmps.be/fr/usage_humain/produits_particuliers/substances_specialement_reglementees/precurseurs/declaration_du",
                target="_blank"
            )),
            html.Li(html.A(
                "Import/Export authorization",
                href="https://www.afmps.be/fr/humain/produits_particuliers/subst_specialement_reglementees/precurseurs/autorisation_import_export",
                target="_blank"
            ))
        ])
    ], style={"marginTop": "30px"})

    return html.Div(blocks + [shared_links])
