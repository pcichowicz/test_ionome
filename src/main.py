import pprint

import pandas as pd
import pymzml
from matplotlib import pyplot as plt
from xml.etree import ElementTree as ET

import numpy as np
from src.ionome_core import Ionome

from scipy.sparse.linalg import spsolve
from scipy import sparse


def main():
    # ----- Init -----
    pd.set_option('display.max_columns', None)
    first = Ionome(run_id="SL2031", samples="samples_SL2031.yaml")

    #----------------------------------------------------------------

    # ----- Load data -----
    first.load_data()


    ## ----------------------- ##
    ## Metadata from mzML
    ## ------------------------##
    file_xml = "../data/raw_SL2031/018_20230825_SL2031__EL-cat_MS1_neg.mzML"
    # ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}
    #
    # tree = ET.parse(file_xml)
    # root = tree.getroot()
    #
    # metadata = {}
    # mmzl_root = root.find('mzml:mzML', ns)
    # file_description = root.find('mzml:fileDescription', ns)
    #
    # if mmzl_root is not None:
    #     file_description = mmzl_root.find('mzml:fileDescription', ns)
    #
    #     if file_description is not None:
    #         metadata['file_description'] = {}
    #         for cv_param in file_description.findall('.//mzml:cvParam', ns):
    #             name = cv_param.attrib.get('name')
    #             value = cv_param.attrib.get('value', 'N/A')
    #             metadata['file_description'][name] = value
    #
    #
    # pprint.pp(metadata)


    # configs = []
    # for inst in root.findall('.//mzml:instrumentConfiguration', ns):
    #     config_id = inst.attrib.get('id')
    #     components = {}
    #
    #     for comp in inst.findall('.//mzml:component', ns):
    #         comp_type = comp.attrib.get('order')
    #         cv = comp.find('mzml:cvParam', ns)
    #         if cv is not None:
    #             components[comp_type] = cv.attrib.get('name')
    #
    #     configs.append({
    #         "id": config_id,
    #         "components": components
    #     })
    # print(configs)
    #
    # source = []
    # for src in root.findall('.//mzml:source', ns):
    #     params = {}
    #     for cv in src.findall('mzml:cvParam', ns):
    #         params[cv.attrib['name']] = cv.attrib.get('value')
    #     source.append(params)
    # print(source)
    #
    # analyzers = []
    # for an in root.findall('.//mzml:analyzer', ns):
    #     params = {}
    #     for cv in an.findall('mzml:cvParam', ns):
    #         params[cv.attrib['name']] = cv.attrib.get('value')
    #     analyzers.append(params)
    # print(analyzers)
    #
    # detectors = []
    # for det in root.findall('.//mzml:detector', ns):
    #     params = {}
    #     for cv in det.findall('mzml:cvParam', ns):
    #         params[cv.attrib['name']] = cv.attrib.get('value')
    #     detectors.append(params)
    # print(detectors)
    # ----------------------------------------------------------------

    # ----- Quality control -----
    first.extract_quality_control()

    # __ TIC __
    plt.figure(figsize=(10, 6))
    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            plt.plot(sampleData.quality_control["retention_time"], sampleData.quality_control["tic"], alpha=0.5, label=sampleData.unique_id)

    plt.title(f"Total Ion Chromatograms - Overlay")
    plt.xlabel("Retention time (min)")
    plt.ylabel("Total Ion Intensity")
    plt.legend()
    plt.tight_layout()
    plt.show()
    #
    # # __ Peaks per scan - richness/noise __
    plt.figure(figsize=(10, 3))
    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            plt.plot(sampleData.quality_control['retention_time'], sampleData.quality_control['peaks_per_scan'], label = sampleData.unique_id)

    plt.title("Peaks per scan")
    plt.xlabel("Retention time (min)")
    plt.ylabel("Number of peaks")
    plt.grid(True, linestyle="-", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # __ BPC & bpc_mz overlay __
    fig, ax = plt.subplots(2,1, figsize=(10,6), sharex=True)
    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            ax[0].plot(sampleData.quality_control["retention_time"], sampleData.quality_control['bpc'], label = sampleData.unique_id)
            ax[0].set_ylabel("Base peak intensity")
            ax[1].scatter(sampleData.quality_control["retention_time"], sampleData.quality_control["bpc_mz"], c=sampleData.quality_control["bpc"], cmap='viridis', s=6)
            ax[1].set_ylabel("Base peak m/z")

            ax[1].set_xlabel("Retention time (min)")
    # fig.colorbar(im, ax=ax[1], label='BPC intensity')
    plt.legend()
    plt.tight_layout()
    plt.show()
    # ----------------------------------------------------------------

    # ----- XIC -----
    first.extract_ion_chromatograms()

    # __ XIC Plot __
    plt.figure(figsize=(10, 6))
    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            plt.plot(sampleData.xic['catechin']["retention_time"], sampleData.xic['catechin']['intensity'], alpha=0.5,
                     label=sampleData.unique_id)

    plt.title(f"XIC - Overlay")
    plt.xlabel("Retention time (min)")
    plt.ylabel("Intensity")
    plt.legend()
    plt.tight_layout()
    plt.show()
if __name__ == "__main__":
    main()