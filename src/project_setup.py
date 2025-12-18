import os
import shutil
import sys
import argparse
from pprint import pprint

import pandas as pd
import yaml
from pathlib import Path
from src.helpers import log_method_entry
from src.utils import samples_df_to_yaml, load_samples_table

## ------------------- ##
## Helpers
## ------------------- ##

def load_yaml(path):
    with open(path, 'r') as file:
        return yaml.safe_load(file)

def write_yaml_config(path, data):
    with open(path, 'w') as file:
        yaml.dump(data,file,default_flow_style=False, sort_keys=False)

def write_yaml_samples(path, data):
    with open(path, 'w') as file:
        return yaml.dump(data,default_flow_style=False, sort_keys=False)


def write_samples_yaml(df: pd.DataFrame, out_path: Path):
    data = samples_df_to_yaml(df)
    sample_yaml = write_yaml_samples(out_path, data)
    sample_yaml = sample_yaml.replace('\n-', '\n\n-')

    with open(out_path, 'w') as file:
        file.write(sample_yaml)

## ------------------------------- ##
## Project Setup Class
## ------------------------------- ##

class IonomeProjectSetup:
    def __init__(self,
                 project_dir: str | Path = 'projects',
                 configs_dir:str | Path = 'config'
                 ):

        log_method_entry()

        self.base_dir = Path(__file__).resolve().parent.parent
        self.project_dir = self.base_dir / project_dir
        self.templates_dir = self.base_dir / configs_dir

        self.project_dir.mkdir(parents=True, exist_ok=True)

        if not self.templates_dir.exists():
            raise FileNotFoundError(f'Template directory for pipeline configurations does not exist: {self.templates_dir}')

        # print(f"\t[>] Project setup initialized")

    def create_project(self,
                       project_name: str,
                       overwrite: bool = False,
                       data_path:str | Path | None = None,
                       sample_table: str | Path | None = None) -> Path:
        """
        Creates a new project directory with configuration templates.

        """
        log_method_entry()

        project_root = self.project_dir / project_name

        if project_root.exists():
            if not overwrite:
                raise FileExistsError(f"Project directory already exists: {project_root}\n"
                                      f"If you want to overwrite this project, delete the existing project directory or set the overwrite flag to 'True'.")
            shutil.rmtree(project_root)


        project_root.mkdir(parents=True)
        (project_root / "raw_data").mkdir()
        (project_root / "processed").mkdir()
        (project_root / "results").mkdir()
        (project_root / "logs").mkdir()

        config_template = load_yaml(self.templates_dir / f"config.yaml")
        cfg_path = project_root / f"config_{project_name}.yaml"
        write_yaml_config(cfg_path, config_template)

        # THIS SHOULD BE sample_table, however for testing using a dir within project

        if not sample_table:
            sample_template = load_yaml(self.templates_dir / f"samples_template.yaml")
            samp_path = project_root / f"samples_{project_name}.yaml"
            write_yaml_samples(samp_path, sample_template)

        sample_table_path = self.base_dir / "input_out" / sample_table
        samp_out = project_root / f"samples_{project_name}.yaml"
        sample_template = load_samples_table(sample_table_path)
        write_samples_yaml(sample_template, samp_out)

        if data_path:
            print(f"Checking for data files to copy.")

        print(f"\t[>] Project '{project_name}' created successfully.")
        print(f"\t Follow these steps:")
        print(f"\t [1.] Before running 'Ionome', manually copy data files into the raw_data folder/directory.")
              # f"\t\t  Or provide full path and file extention where data files are located to be copied automatically.")
        print(f"\t [2.] Check and edit the config_{project_name}.yaml file if necessary.")
        print(f"\t [3.] Confirm the samples_{project_name}.yaml file has all the samples data.")
        print(f"\t [4.] Run Ionome('{project_name}', '{samp_out.stem}')")
        return project_root

    def download_data(self):
        print(f"Downloads data from database.")
##==============================##
## Command Line Interface (CLI) ##
##==============================##

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create a new Ionome LC-MS project directory.'
    )

    parser.add_argument('-p', '--project', type=str, help="Project directory name")

    args = parser.parse_args()

    setup = IonomeProjectSetup()
    setup.create_project('SL2031')
