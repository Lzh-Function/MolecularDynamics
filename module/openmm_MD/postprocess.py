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
import MDAnalysis.transformations as trans
import nglview as nv

import numpy as np
import time
import os
from tqdm.notebook import tqdm
import pickle


def centering(lig_name, pdb, dcd, args):
    """
    input: ligs_name, pdb(path), dcd(path), args
    output(to dir): protein-centered pdb/dcd
    """
    u = mda.Universe(pdb, dcd)
    prot = u.select_atoms("protein")
    lig = u.select_atoms("resname UNK")
    if args.metal != None:
        metal = u.select_atoms(f"resname {args.metal}")
        transforms = [trans.center_in_box(prot), trans.wrap(lig), trans.wrap(metal)]
        u.trajectory.add_transformations(*transforms)

        prot = u.select_atoms("protein")
        lig = u.select_atoms("resname UNK")
        metal = u.select_atoms(f"resname {args.metal}")
        complex = prot + lig + metal
        complex.write(args.results_dir_path / "centered" / f"{lig_name}.pdb")
        complex.write(args.results_dir_path / "centered" / f"{lig_name}.dcd", frames="all")
    else:
        transforms = [trans.center_in_box(prot), trans.wrap(lig)]
        u.trajectory.add_transformations(*transforms)

        prot = u.select_atoms("protein")
        lig = u.select_atoms("resname UNK")
        complex = prot + lig
        complex.write(args.results_dir_path / "centered" / f"{lig_name}.pdb")
        complex.write(args.results_dir_path / "centered" / f"{lig_name}.dcd", frames="all")
    