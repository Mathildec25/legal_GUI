import pandas as pd
from pathlib import Path
import csv

def load_data(path="data/donnees.xlsx"):
    # lit la première feuille du fichier Excel
    df = pd.read_excel(path, dtype=str,engine="openpyxl")  
    return df

df = load_data()
