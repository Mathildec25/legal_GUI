import pandas as pd

def load_clean_data(path="data/donnees.csv"):
    df = pd.read_csv(path, sep=";", encoding="utf-8")
    df.columns = [col.lower().replace(" ", "").strip() for col in df.columns]
    return df
