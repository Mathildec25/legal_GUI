def get_stup_psy_procedure(reg):
    """
    Retourne la procédure juridique applicable pour les stupéfiants/psychotropes.
    Ignore les valeurs vides ou nulles de 'belgian_annex'.
    """
    annex = str(reg.get("belgian_annex", "")).strip().lower()

    if annex in ["", "nan", "none", "null"]:
        return ""  # Rien à signaler : ce n’est pas un stupéfiant

    normalized = annex.split()[0:2]
    joined = " ".join(normalized)

    if joined in ["annex i", "annex ia", "annex ib", "annex ic"]:
        return "**Belgian Annex I**\n❗ Strict control: classified as a narcotic (UN Convention of 1961 – Schedules I to IV)"

    elif joined in ["annex ii", "annex iia", "annex iib"]:
        return "**Belgian Annex II**\nℹ Psychotropic substances under control (UN Convention of 1971 – Schedules I and II + Council Decision 2004/757/JHA)"

    elif joined == "annex iii":
        return "**Belgian Annex III**\n⚠ Less-restricted psychotropics (UN Convention of 1971 – Schedules III and IV)"

    else:
        return f"**Belgian Annex**: {annex} – Not yet handled"
