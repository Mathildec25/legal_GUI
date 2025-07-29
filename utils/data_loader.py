# utils/data_loader.py

import pandas as pd

def load_data(path: str = "C:/Users/user/legal_GUI/data/donnees.csv") -> pd.DataFrame:

    """
    Charge le fichier CSV des substances, nettoie les noms de colonnes.

    :param path: Chemin vers le fichier CSV
    :return: DataFrame nettoyé
    """
    df = pd.read_csv(path, sep=",", encoding="utf-8")

    # Nettoyage des noms de colonnes
    df.columns = [col.lower().replace(" ", "_").strip() for col in df.columns]

    return df

# Chargement immédiat à l'import
df = load_data()
