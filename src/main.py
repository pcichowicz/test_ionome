import pandas as pd
import pymzml
from matplotlib import pyplot as plt

from preprocess import MzmlParser
from ionome_core import Ionome


def main():
    pd.set_option('display.max_columns', None)
    first = Ionome(run_id="SL005", samples="samples_SL011.yaml")
    first.load_data()
    first.extract_xic()
    # first.correct_baseline()
    first.peak_detection_xic()


    # df_24hr_coumerate = first.df_xic["SL005_24hr"]["p-coumerate"]
    # df_48hr_coumerate = first.df_xic["SL005_48hr"]["p-coumerate"]
    # df_72hr_coumerate = first.df_xic["SL005_72hr"]["p-coumerate"]
    #
    # df_24hr_urocanate = first.df_xic["SL005_24hr"]["urocanate"]
    # df_48hr_urocanate = first.df_xic["SL005_48hr"]["urocanate"]
    # df_72hr_urocanate = first.df_xic["SL005_72hr"]["urocanate"]
    #
    # df_24hr_caffeate = first.df_xic["SL005_24hr"]["caffeate"]
    # df_48hr_caffeate = first.df_xic["SL005_48hr"]["caffeate"]
    # df_72hr_caffeate = first.df_xic["SL005_72hr"]["caffeate"]

    # fig,ax = plt.subplots(figsize=(12,8))
    # ax.plot(df_24hr_coumerate['retention_time'], df_24hr_coumerate['intensity'], color='green', label="Intensity")
    # ax.plot(df_24hr_coumerate['retention_time'], df_24hr_coumerate['corrected'], color='black', label="Corrected")
    # ax.plot(df_24hr_coumerate['retention_time'], df_24hr_coumerate['baseline'], color='red', label="Baseline")
    #
    # ax.set_xlim((10,11))
    # fig, ax = plt.subplots(1,3)
    # ax[0].plot(df_24hr_coumerate['retention_time'], df_24hr_coumerate['corrected'])
    # ax[0].set_title("24hr coumerate")
    # ax[1].plot(df_48hr_coumerate['retention_time'], df_48hr_coumerate['corrected'])
    # ax[1].set_title("48hr coumerate")
    # ax[2].plot(df_72hr_coumerate['retention_time'], df_72hr_coumerate['corrected'])
    # ax[2].set_title("72hr coumerate")
    #
    # ax[1, 0].plot(df_24hr_urocanate['retention_time'], df_24hr_urocanate['corrected'])
    # ax[1,0].set_title("24hr urocanate")
    # ax[1, 1].plot(df_48hr_urocanate['retention_time'], df_48hr_urocanate['corrected'])
    # ax[1,1].set_title("48hr urocanate")
    # ax[1, 2].plot(df_72hr_urocanate['retention_time'], df_72hr_urocanate['corrected'])
    # ax[1,2].set_title("72hr urocanate")
    #
    # ax[2, 0].plot(df_24hr_caffeate['retention_time'], df_24hr_caffeate['corrected'])
    # ax[2, 0].set_title("24hr caffeate")
    # ax[2, 1].plot(df_48hr_caffeate['retention_time'], df_48hr_caffeate['corrected'])
    # ax[2, 1].set_title("48hr caffeate")
    # ax[2, 2].plot(df_72hr_caffeate['retention_time'], df_72hr_caffeate['corrected'])
    # ax[2, 2].set_title("72hr caffeate")

    # plt.tight_layout()
    # plt.show()
if __name__ == "__main__":
    main() 