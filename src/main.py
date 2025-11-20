import pandas as pd
import pymzml
from matplotlib import pyplot as plt

from src.ionome_core import Ionome


def main():
    pd.set_option('display.max_columns', None)
    first = Ionome(run_id="SL005", samples="samples_SL011.yaml")
    first.load_data()
    first.extract_xic()
    first.plot_chromatogram("tic")

if __name__ == "__main__":
    main() 