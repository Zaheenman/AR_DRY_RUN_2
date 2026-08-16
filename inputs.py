from pathlib import Path
import pandas as pd
BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "inputs.csv"

def load_inputs():
    data = pd.read_csv(INPUT_FILE)
    return data
