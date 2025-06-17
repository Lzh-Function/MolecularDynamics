import openmm as mm
from openmm import *
from openmm.app import *
import openmm.app as app
from openmm import unit
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
from rdkit import Chem
from rdkit.Chem import AllChem
import pdbfixer
import mdtraj as md
import MDAnalysis as mda
import nglview as nv

import numpy as np
import time
import os
from pathlib import Path
from tqdm.notebook import tqdm
import pickle
from argparse import ArgumentParser, FileType
import yaml

from preprocess import *
from simulation import *
from postprocess import *

def main(args):
    
    ts1 = time.perf_counter()
    print("preprocess start...")
    
    complex_models, FFs = mkmodels(args)
    
    with open(args.complex_saving_dir_path / "complex_models.pkl", "wb") as f:
        pickle.dump(complex_models, f)
    with open(args.complex_saving_dir_path / "FFs.pkl", "wb") as f:
        pickle.dump(FFs, f)

    ts1 = time.perf_counter()
    
    print("simulation start...")
    print("-" * 60)
    for complex_model, ff, lig_name in zip(complex_models, FFs, args.ligs_name):
        ts2 = time.perf_counter()
        simulate(complex_model, ff, lig_name, args)
        tf2 = time.perf_counter()
        elapse2 = tf2 - ts2
        h = elapse2 // 3600
        m = (elapse2 % 3600) // 60
        s = elapse2 % 60
        print(f"elapsed time: {h} h {m} min {s} sec")
        print("-" * 60)
    
    print("postprocess start...")
    os.makedirs(args.results_dir_path / "centered", exist_ok=True)
    pdbs = []
    dcds = []
    for lig_name in args.ligs_name:
        pdbs.append(args.results_dir_path / f"{lig_name}" / f"{lig_name}_equilibrated.pdb")
        dcds.append(args.results_dir_path / f"{lig_name}" / f"{lig_name}_trajectory.dcd")
    for lig_name, pdb, dcd in zip(args.ligs_name, pdbs, dcds):
        centering(lig_name, pdb, dcd, args)
    print("-" * 60)
    
    tf1 = time.perf_counter()
    elapse1 = tf1 - ts1
    h = elapse1 // 3600
    m = (elapse1 % 3600) // 60
    s = elapse1 % 60
    print(f"all finished. elapsed time: {h} h {m} min {s} sec")
    
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--config", type=FileType(mode="r"), default=None)
    args = parser.parse_args()
    config_dict = yaml.load(args.config, Loader=yaml.FullLoader)
    arg_dict = args.__dict__
    for key, value in config_dict.items():
        arg_dict[key] = value
    args.config = Path(args.config.name)
    
    experiment_dir = args.config.parent
    args.complex_saving_dir_path = experiment_dir / "complexes"
    os.makedirs(args.complex_saving_dir_path, exist_ok=True)
    args.config_dir_path = experiment_dir / "configs"
    os.makedirs(args.config_dir_path, exist_ok=True)
    args.results_dir_path = experiment_dir / "results"
    os.makedirs(args.results_dir_path, exist_ok=True)
    args.pbsa_input_dir_path = experiment_dir / "pbsa" / "input"
    os.makedirs(args.pbsa_input_dir_path, exist_ok=True)
    
    print("===== config =====")
    for key, value in vars(args).items():
        print(f"{key}: {value}")
    print("===== config =====")
    
    main(args)