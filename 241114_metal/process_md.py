import os
import glob
from argparse import ArgumentParser, FileType
import yaml
import time
from copy import deepcopy

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import numpy as np
import pandas as pd

import openmm as mm
import openmm.app as app
from openmm import unit
import pdbfixer
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
import mdtraj as md

from mymodeller import MyModeller
#from extract_ligand_pymol import extract_ligand


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--config",type=FileType(mode="r"),default=None)
    args = parser.parse_args()
    config_dict = yaml.load(args.config,Loader=yaml.FullLoader)
    arg_dict = args.__dict__
    for key, value in config_dict.items():
        arg_dict[key] = value
    args.config = args.config.name
    args.experiment_dir = "/".join(args.config.split("/")[:-1])

    args.nonbondedMethod = app.PME
    args.constraints = app.HBonds 
    args.temperature = args.temperature * unit.kelvin
    args.friction = args.friction / unit.picoseconds
    args.pressure = args.pressure * unit.atmospheres
    args.padding = (args.padding / 1e-9) * unit.nanometers
    args.dt = (args.dt / 1e-15) * unit.femtoseconds
    args.nonbondedCutoff = (args.nonbondedCutoff / 1e-9) * unit.nanometers
    return args

def prepare_protein(protein_file,ph,metal_file=None):
    fixer = pdbfixer.PDBFixer(protein_file)
    fixer.removeHeterogens()
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in list(keys):
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    modeller = MyModeller(fixer.topology,fixer.positions)

    if metal_file != None:
        metal_fixer = pdbfixer.PDBFixer(metal_file)
        metal_modeller = MyModeller(metal_fixer.topology,metal_fixer.positions)
        new_top, new_pos = add_modeller_withCharge(modeller,metal_modeller)
        modeller = MyModeller(new_top,new_pos)
    return modeller

def _add_modeller_withCharge(topology,positions,addModeller):
    atop = addModeller.topology
    apos = addModeller.positions

    newAtoms = {}
    for chain in atop.chains():
        newChain = topology.addChain(chain.id)
        for residue in chain.residues():
            newResidue = topology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
            for atom in residue.atoms():
                newAtom = topology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                newAtoms[atom] = newAtom
                positions.append(deepcopy(apos[atom.index]))
    for bond in atop.bonds():
        topology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)
    
    return topology, positions

def add_modeller_withCharge(*addModeller):
    new_top = app.Topology()
    new_pos = [] * unit.nanometers

    for i,v in enumerate(addModeller):
        if i == 0:
            new_top.setPeriodicBoxVectors(v.topology.getPeriodicBoxVectors())
        new_top, new_pos = _add_modeller_withCharge(new_top,new_pos,v)
    return new_top, new_pos

def prepare_ligand(pdb_file,ligand_file):
    rdkit_mol = Chem.MolFromPDBFile(pdb_file)
    prepared_ligand = Chem.SDMolSupplier(ligand_file)
    prepared_ligand = [v for v in prepared_ligand if v != 0][0]
    prepared_ligand = Chem.rdmolops.AddHs(prepared_ligand,addCoords=True)

    off_mol = Molecule.from_rdkit(prepared_ligand)
    element_counter_dict = {}
    for off_atom, rdkit_atom in zip(off_mol.atoms,rdkit_mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        if element in element_counter_dict.keys():
            element_counter_dict[element] += 1
        else:
            element_counter_dict[element] = 1
        off_atom.name = element + str(element_counter_dict[element])     
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    mol_positions = off_mol.conformers[0].to_openmm()
    mod_mol = MyModeller(mol_topology,mol_positions)
    return mod_mol, off_mol   

def prep_model(args,top,pos,off_mol,write_file="complex_model"):
    FF = app.ForceField(args.protein_ff,args.solvent_ff)
    smff = SMIRNOFFTemplateGenerator(molecules=off_mol)
    FF.registerTemplateGenerator(smff.generator)
    model = MyModeller(top,pos)
    model.addSolvent(FF,padding=args.padding,ionicStrength=args.ionicStrength * unit.molar)
    if len(write_file) > 0:
        app.PDBFile.writeFile(model.getTopology(),model.getPositions(),open("data/{}.pdb".format(write_file),"w"))
    return model, FF

def process(args,ligand_file):
    mod_pro = prepare_protein(args.protein_file,args.ph,args.metal_file)
    mod_mol, off_mol = prepare_ligand(args.protein_file,ligand_file)
    top, pos = add_modeller_withCharge(mod_pro,mod_mol)
    model, FF = prep_model(args,top,pos,off_mol)
    return model, FF

def run(args,model,FF,result_dir):
    platform = mm.Platform.getPlatformByName("CUDA" if args.cuda else "cpu")
    platformProperties = {"Precision":"single"}

    system = FF.createSystem(model.topology,
                             nonbondedMethod=args.nonbondedMethod,
                             nonbondedCutoff=args.nonbondedCutoff,
                             constraints=args.constraints,
                             rigidWater=args.rigidWater,
                             ewaldErrorTolerance=args.ewaldErrorTolerance)
    system.addForce(mm.MonteCarloBarostat(args.pressure,args.temperature,args.barostatInterval))

    integrator = mm.LangevinMiddleIntegrator(args.temperature,args.friction,args.dt)
    integrator.setConstraintTolerance(args.constraintTolerance)
    simulation = app.Simulation(model.topology,system,integrator,platform,platformProperties)
    simulation.context.setPositions(model.positions)

    with open(os.path.join(result_dir,"sysmtem.xml"),mode="w") as f:
        f.write(mm.XmlSerializer.serialize(system))
    with open(os.path.join(result_dir,"integrator.xml"),mode="w") as f:
        f.write(mm.XmlSerializer.serialize(integrator))

    print("{} simulation start".format(result_dir.split("/")[-1]))
    simulation.minimizeEnergy()
    with open(os.path.join(result_dir,"minimized.pdb"),"w") as f:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True,enforcePeriodicBox=False).getPositions(),
            file=f,
            keepIds=True
        )
    ts = time.perf_counter()
    simulation.context.setVelocitiesToTemperature(args.temperature)
    simulation.step(args.equilibrationSteps)

    dcdReporter = app.DCDReporter(os.path.join(result_dir,"trajectory.dcd"),args.reporterStep,enforcePeriodicBox=True)
    dataReporter = app.StateDataReporter(os.path.join(result_dir,"log.txt"),args.reporterStep,totalSteps=args.steps,step=True,speed=True,progress=True,potentialEnergy=True,temperature=True,separator="\t")
    checkpointReporter = app.CheckpointReporter(os.path.join(result_dir,"checkpoint.chk"),args.checkpointStep)
    simulation.reporters.append(dcdReporter)
    simulation.reporters.append(dataReporter)
    simulation.reporters.append(checkpointReporter)
    simulation.currentStep = 0

    simulation.step(args.steps)
    tg = time.perf_counter()
    simulation.saveState(os.path.join(result_dir,"final_state.xml"))
    state = simulation.context.getState(getPositions=True,enforcePeriodicBox=True)
    with open(os.path.join(result_dir,"final_state.cif"),mode="w") as file:
        app.PDBxFile.writeFile(simulation.topology,state.getPositions(),file)
    print("{} simulation finish".format(result_dir.split("/")[-1]))
    elapse = tg - ts
    h = elapse // 3600
    m = (elapse % 3600) // 60
    s = elapse % 60
    print(f"equilibration - simulation elapsed time: {h} h {m} min {s} sec")



def main():
    # make metal ion file in advance
    args = get_args()
    ligand_path = glob.glob("ligand/*.sdf")
    idx = [v.split(".")[0].split("_")[-1] for v in ligand_path]
    for v,w in zip(ligand_path,idx):
        result_dir = "result/{}".format(w)
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        print("{} process start".format(w))
        model, FF = process(args,v)
        run(args,model,FF,result_dir)


if __name__ == "__main__":
    main()