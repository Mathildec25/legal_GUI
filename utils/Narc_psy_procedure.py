from dash import html

def get_stup_psy_procedure(reg):
    """
    Returns the legal procedure applicable to narcotic/psychotropic substances,
    including icons (via Font Awesome), official AFMPS links, and storage rules,
    based on the Belgian annex (Annex I, II, III...).

    Parameter:
    - reg (dict): dictionary containing the 'belgian_annex' attribute.

    Returns:
    - html.Div: Dash HTML block describing the applicable legal obligations.
    """
    # Retrieve and normalize the Belgian annex
    annex = str(reg.get("belgian_annex", "")).strip().lower()

    # If no annex is specified, return nothing
    if annex in ["", "nan", "none", "null"]:
        return ""

    # Normalization: keep only the first two words of the annex
    normalized = annex.split()[0:2]
    joined = " ".join(normalized)

    #  Belgian Annex I
    # ============================
    if joined in ["annex i", "annex ia", "annex ib", "annex ic"]:
        return html.Div([
            html.H5([
                html.I(className="fas fa-gavel me-2 text-warning"),  # Icône via Font Awesome
                "Belgian Annex I"
            ], className="fw-bold"),

            html.P([
                html.I(className="fas fa-exclamation-triangle me-2 text-danger"),
                "Strict control: classified as a narcotic (UN Convention of 1961 – Schedules I to IV)"
            ]),

            html.P([
                html.I(className="fas fa-warehouse me-2 text-secondary"),
                "Storage: Must be kept in a closed space secured against theft (Article 40, §1 of the Royal Decree 2017)."
            ]),

            html.P([
                html.I(className="fas fa-id-badge me-2 text-info"),
                "The activity license is valid for 3 years from the date of issuance and can be renewed up to 3 months before expiration."
            ]),

            # Authorization links (activity / end-user)
            html.Ul([
                html.Li([
                    html.I(className="fas fa-file-signature me-2 text-primary"),
                    html.A("Activity Authorization", href="https://www.afmps.be/fr/humain/produits_particuliers/subst_specialement_reglementees/stupefiants_et_psychotropes/autorisation_d_activites", target="_blank")
                ]),
                html.Li([
                    html.I(className="fas fa-user-check me-2 text-primary"),
                    html.A("End-User Authorization", href="https://www.afmps.be/fr/humain/produits_particuliers/subst_specialement_reglementees/stupefiants_et_psychotropes/autorisation_utilisateur_final", target="_blank")
                ])
            ]),

            # Mandatory registration in Narcoreg
            html.P([
                html.I(className="fas fa-database me-2 text-success"),
                "All transactions involving narcotics or psychotropic substances listed in Annexes I, II and IV (except Ic and IVc preparations) of the Royal Decree of 6 September 2017 must be recorded in the AFMPS online system: "
            ]),
            html.A("AFMPS Narcotic Registration Platform", href="https://www.afmps.be/fr/narcoreg", target="_blank", className="btn btn-outline-info btn-sm mt-2")
        ])

    #  Belgian Annex II
    # ============================
    elif joined in ["annex ii", "annex iia", "annex iib","annex iic"]:
        return html.Div([
            html.H5([
                html.I(className="fas fa-gavel me-2 text-warning"),
                "Belgian Annex II"
            ], className="fw-bold"),

            html.P([
                html.I(className="fas fa-info-circle me-2 text-primary"),
                "Psychotropic substances under control (UN Convention of 1971 – Schedules I and II + Council Decision 2004/757/JHA)"
            ]),

            html.P([
                html.I(className="fas fa-warehouse me-2 text-secondary"),
                "Storage: Must be kept in a closed space secured against theft (Article 40, §1 of the Royal Decree 2017)."
            ]),

            html.P([
                html.I(className="fas fa-id-badge me-2 text-info"),
                "The activity license is valid for 3 years from the date of issuance and can be renewed up to 3 months before expiration."
            ]),

            html.Ul([
                html.Li([
                    html.I(className="fas fa-file-signature me-2 text-primary"),
                    html.A("Activity Authorization", href="https://www.afmps.be/fr/humain/produits_particuliers/subst_specialement_reglementees/stupefiants_et_psychotropes/autorisation_d_activites", target="_blank")
                ]),
                html.Li([
                    html.I(className="fas fa-user-check me-2 text-primary"),
                    html.A("End-User Authorization", href="https://www.afmps.be/fr/humain/produits_particuliers/subst_specialement_reglementees/stupefiants_et_psychotropes/autorisation_utilisateur_final", target="_blank")
                ])
            ]),

            html.P([
                html.I(className="fas fa-database me-2 text-success"),
                "All transactions involving narcotics or psychotropic substances listed in Annexes I, II and IV (except Ic and IVc preparations) of the Royal Decree of 6 September 2017 must be recorded in the AFMPS online system: "
            ]),
            html.A("AFMPS Narcotic Registration Platform", href="https://www.afmps.be/fr/narcoreg", target="_blank", className="btn btn-outline-info btn-sm mt-2")
        ])

    # Belgian Annex III
    # ============================
    elif joined == "annex iii":
        return html.Div([
            html.H5([
                html.I(className="fas fa-gavel me-2 text-warning"),
                "Belgian Annex III"
            ], className="fw-bold"),

            html.P([
                html.I(className="fas fa-exclamation-circle me-2 text-warning"),
                "Less-restricted psychotropics (UN Convention of 1971 – Schedules III and IV)"
            ]),

            html.P([
                html.I(className="fas fa-box-open me-2 text-secondary"),
                "Storage: Must be kept in a suitably arranged space offering protection against breakage or theft (Article 40, §2 of the Royal Decree 2017)."
            ])
        ])
    #  Unsupported case
    # ============================
    else:
        return html.Div([
            html.H5([
                html.I(className="fas fa-question-circle me-2 text-muted"),
                f"Belgian Annex: {annex}"
            ]),
            html.P("– Not yet handled")
        ])
