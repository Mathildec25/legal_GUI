def get_legal_procedure(reg):
    """
    Retourne la procédure juridique applicable selon les colonnes du fichier réglementaire.

    :param reg: Ligne du DataFrame réglementaire (une substance)
    :return: Texte indiquant la procédure à suivre
    """
    cat_273 = str(reg.get('eu_annex/category_273/2004', '')).strip().lower()
    cat_111 = str(reg.get('eu_annex/category_111/2005', '')).strip().lower()

    msg_273 = "✔ No specific regulatory procedure identified (EU 273/2004)"
    msg_111 = "✔ No specific regulatory procedure identified (EU 111/2005)"

    # -------- EU 273/2004 --------
    procedures_273 = []

    if "annex i, category 1" in cat_273:
        procedures_273.append(" License required (Annex I, Category 1)")
    if "annex i, subcategory 2a" in cat_273 or "annex i, subcategory 2b" in cat_273 or "annex i, category 2" in cat_273:
        procedures_273.append("ℹRegistration required (Annex I, Category 2)")
    if "annex i, category 3" in cat_273:
        procedures_273.append("Declaration required (Annex I, Category 3)")

    if procedures_273:
        msg_273 = "**EU 273/2004**\n" + "\n".join(procedures_273)

    # -------- EU 111/2005 --------
    procedures_111 = []

    if "annex, category 1" in cat_111:
        procedures_111.append(" Import/Export License required (Annex, Category 1)")
    if "annex, category 2" in cat_111:
        procedures_111.append(" Registration required (Annex, Category 2)")
    if "annex, category 3" in cat_111:
        procedures_111.append(" Declaration required (Annex, Category 3)")

    if procedures_111:
        msg_111 = "**EU 111/2005**\n" + "\n".join(procedures_111)

    return f"{msg_273}\n\n{msg_111}"


