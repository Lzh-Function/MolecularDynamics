import os
import glob
from argparse import ArgumentParser, FileType
import yaml
import time

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

def process_protein(protein_file):
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
    ph = 7.0
    fixer.addMissingHydrogens(ph)
    return fixer

def extract_metal(protein_file,metal_name):
    fixer = pdbfixer.PDBFixer(protein_file)
    mod_metal = app.Modeller(fixer.topology,fixer.positions)
    toDelete = []
    for v in fixer.topology.residues():
        if v.name != metal_name:
            toDelete.append(v)
    mod_metal.delete(toDelete)
    return mod_metal

def process_ligand(pdb_file,ligand_file):
    """
    rdkit_mol = Chem.MolFromPDBFile(pdb_file)
    ligand = Chem.rdmolops.SplitMolByPDBResidues(rdkit_mol)["LG1"]
    ligand = Chem.RemoveHs(ligand)
    reference = Chem.MolFromSmiles(smiles)
    prepared_ligand = AllChem.AssignBondOrdersFromTemplate(reference,ligand)
    prepared_ligand.AddConformer(ligand.GetConformer(0))
    prepared_ligand = Chem.rdmolops.AddHs(prepared_ligand,addCoords=True)
    prepared_ligand = Chem.MolFromMolBlock(Chem.MolToMolBlock(prepared_ligand))
    """
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
    mol_positions = off_mol.conformers[0].to("nanometers")
    mod_mol = app.Modeller(mol_topology,mol_positions)
    return mod_mol, off_mol         

def complex(*object):
    topology = [md.Topology.from_openmm(v.topology) for v in object]
    complex_top = topology[0]
    for v in topology[1:]:
        complex_top = complex_top.join(v)
    atoms = [len(v.positions) for v in object]
    complex_pos = unit.Quantity(np.zeros([sum(atoms),3]), unit=unit.nanometers)
    k = 0
    for v,w in zip(object,atoms):
        complex_pos[k:k+w] = v.positions
        k += w
    return complex_top, complex_pos

def prep_model(args,complex_top,complex_pos,off_mol,write_file=""):
    FF = app.ForceField(*args.forcefield)
    smff = SMIRNOFFTemplateGenerator(molecules=off_mol)
    FF.registerTemplateGenerator(smff.generator)
    model = app.Modeller(complex_top.to_openmm(),complex_pos)
    model.addSolvent(FF,padding=args.padding,ionicStrength=args.ionicStrength * unit.molar)
    if len(write_file) > 0:
        app.PDBFile.writeFile(model.getTopology(),model.getPositions(),open("data/{}.pdb".format(write_file),"w"))
    return model, FF

def process(args,ligand_file):
    fixer = process_protein(args.protein_file)
    mod_mol, off_mol = process_ligand(args.protein_file,ligand_file)
    if args.metal:
        mod_metal = extract_metal(args.protein_file,args.metal_name)
        top, pos = complex(fixer,mod_mol,mod_metal)
    else:
        top, pos = complex(fixer,mod_mol)

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
    args = get_args()
    ligand_path = glob.glob("ligand/*.sdf")[5:]
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