# utils/data_loader.py

import pandas as pd

def load_data(path: str = "C:/Users/user/legal_GUI/data/donnees.csv") -> pd.DataFrame:

    """
    Loads the CSV file of substances and cleans the column names.

    :param path: Path to the CSV file
    :return: Cleaned DataFrame
    """
    df = pd.read_csv(path, sep=",", encoding="utf-8")

    # Clean column names: lowercase, replace spaces with underscores, strip extra spaces
    df.columns = [col.lower().replace(" ", "_").strip() for col in df.columns]

    return df

# Automatically load data on import
df = load_data()
