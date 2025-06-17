import openmm as mm
from openmm import *
from openmm.app import *
import openmm.app as app
from openmm import unit
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from rdkit import Chem
from rdkit.Chem import AllChem
from pdbfixer import PDBFixer
import mdtraj as md
import MDAnalysis as mda
import nglview as nv

import numpy as np
import time
import os
from tqdm.notebook import tqdm
import pickle
from argparse import Namespace


def prot_prep(args:Namespace):
    """
    input: args
    output: pdbfixer(preprocessed protein), mod_metal(metal)
    """
    protein = PDBFixer(args.protein_path)
    protein.removeHeterogens()
    protein.findMissingResidues()

    chains = list(protein.topology.chains())
    keys = protein.missingResidues.keys()
    for key in list(keys):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del protein.missingResidues[key]
    protein.findNonstandardResidues()
    protein.replaceNonstandardResidues()
    protein.findMissingAtoms()
    protein.addMissingAtoms()
    protein.addMissingHydrogens(args.ph)
    protein = app.Modeller(protein.topology, protein.positions)

    if args.metal != None:
        fixer = PDBFixer(args.protein_path)
        mod_metal = app.Modeller(fixer.topology, fixer.positions)
        toDelete = []
        for v in fixer.topology.residues():
            if v.name != args.metal:
                toDelete.append(v)
        mod_metal.delete(toDelete)
        protein.add(mod_metal.topology, mod_metal.positions)
    
    return protein

def listup_ligs(args:Namespace):
    """
    input: args
    output: ligs(list of the ligand path), ligs_name(list of the ligand name)
    """
    ligs = os.listdir(args.ligand_dir_path)
    ligs = sorted([os.path.join(args.ligand_dir_path, lig) for lig in ligs])
    ligs_name = sorted([lig.replace(".sdf", "") for lig in os.listdir(args.ligand_dir_path)])

    return ligs, ligs_name

def lig_prep(lig_path, lig_name, protein:Modeller):
    """
    input: ligand path(sdf), ligand name
    output: omm_mol(openmm object), off_mol(openff object)
    """
    mol = Chem.SDMolSupplier(lig_path)[0]
    mol = Chem.RemoveHs(mol)
    mol = Chem.rdmolops.AddHs(mol, addCoords=True)
    mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))

    off_mol = Molecule.from_rdkit(mol)
    off_mol.name = lig_name
    element_counter_dict = {}
    for off_atom, rdkit_atom in zip(off_mol.atoms, mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        if element in element_counter_dict.keys():
            element_counter_dict[element] += 1
        else:
            element_counter_dict[element] = 1
        off_atom.name = element + str(element_counter_dict[element])
    
    mol_top = off_mol.to_topology().to_openmm()
    mol_pos = off_mol.conformers[0].to_openmm()
    protein.add(mol_top, mol_pos)

    return protein, off_mol

"""
def merging(fixer, omm_mol, mod_metal=None):

    md_protein_top = md.Topology.from_openmm(fixer.topology)
    md_ligand_top = md.Topology.from_openmm(omm_mol.topology)
    if mod_metal != None:
        md_metal_top = md.Topology.from_openmm(mod_metal.topology)
    md_complex_top = md_protein_top.join(md_ligand_top)
    if md_metal_top:
        md_complex_top = md_complex_top.join(md_metal_top)
    complex_top = md_complex_top.to_openmm()

    if md_metal_top:
        total_atoms = len(fixer.positions) + len(omm_mol.positions) + len(mod_metal.positions)
        complex_pos = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)
        complex_pos[:len(fixer.positions)] = fixer.positions
        complex_pos[len(fixer.positions):-len(mod_metal.positions)] = omm_mol.positions
        complex_pos[-len(mod_metal.positions):] = mod_metal.positions
    else:
        total_atoms = len(fixer.positions) + len(omm_mol.positions)
        complex_pos = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)
        complex_pos[:len(fixer.positions)] = fixer.positions
        complex_pos[len(fixer.positions):] = omm_mol.positions
    
    return complex_top, complex_pos
"""

def defineFF(off_mol, args):
    """
    input: off_mol(ligand openff object), args
    output: FF(ForceField)
    """
    FF = app.ForceField(args.protein_ff, args.solvent_ff)
    gaff = GAFFTemplateGenerator(molecules=off_mol)
    FF.registerTemplateGenerator(gaff.generator)

    return FF

def mkmodels(args):
    protein = prot_prep(args)
    ligs, ligs_name = listup_ligs(args)
    args.ligs_name = ligs_name
    
    complex_models = []
    FFs = []

    for lig, lig_name in zip(ligs, ligs_name):
        complex_model, off_mol = lig_prep(lig, lig_name, protein)
        ff = defineFF(off_mol, args)
        water_model = "tip3p" if "tip3p" in args.solvent_ff else "tip4pew"
        complex_model.addSolvent(ff, model=water_model, padding=1.2*unit.nanometers, ionicStrength=0.15*unit.molar)
        top = complex_model.getTopology()
        pos = complex_model.getPositions()
        app.PDBFile.writeFile(top, pos, open(args.complex_saving_dir_path / f"{lig_name}_complex.pdb", "w"))
        complex_models.append(complex_model)
        FFs.append(ff)
    
    return complex_models, FFs