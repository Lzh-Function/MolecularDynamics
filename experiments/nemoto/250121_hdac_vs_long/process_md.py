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

class MDPreparation():
    def __init__(self,args):
        self.args = args
    
    def prepare_protein(self):
        fixer = pdbfixer.PDBFixer(self.args.docked_file)
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
        fixer.addMissingHydrogens(self.args.ph)
        modeller = MyModeller(fixer.topology,fixer.positions)

        if hasattr(self.args,"metal_file"):
            metal_fixer = pdbfixer.PDBFixer(self.args.metal_file)
            metal_modeller = MyModeller(metal_fixer.topology,metal_fixer.positions)
            modeller.add(metal_modeller.topology,metal_modeller.positions)
        
        elif hasattr(self.args,"metal_name"):
            fixer = pdbfixer.PDBFixer(self.args.docked_file)
            mod_metal = MyModeller(fixer.topology,fixer.positions)
            toDelete = []
            for v in fixer.topology.residues():
                if v.name != self.args.metal_name:
                    toDelete.append(v)
            mod_metal.delete(toDelete)
            modeller.add(mod_metal.topology,mod_metal.positions)
        
        return modeller


    def prepare_ligand(self,pdb_file,smiles,resname="LG1"):
        rdkit_mol = Chem.MolFromPDBFile(pdb_file)
        ligand = Chem.rdmolops.SplitMolByPDBResidues(rdkit_mol)[resname]
        ligand = Chem.RemoveHs(ligand)
        reference = Chem.MolFromSmiles(smiles)
        prepared_ligand = AllChem.AssignBondOrdersFromTemplate(reference,ligand)
        prepared_ligand.AddConformer(ligand.GetConformer(0))
        prepared_ligand = Chem.rdmolops.AddHs(prepared_ligand,addCoords=True)
        prepared_ligand = Chem.MolFromMolBlock(Chem.MolToMolBlock(prepared_ligand))

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
    

    def prep_model(self,top,pos,off_mol,write_file="complex_model"):
        FF = app.ForceField(self.args.protein_ff,self.args.solvent_ff)
        smff = SMIRNOFFTemplateGenerator(molecules=off_mol)
        FF.registerTemplateGenerator(smff.generator)
        model = MyModeller(top,pos)
        model.addSolvent(FF,padding=self.args.padding,ionicStrength=self.args.ionicStrength * unit.molar)
        if len(write_file) > 0:
            app.PDBFile.writeFile(model.getTopology(),model.getPositions(),open("{}/{}.pdb".format(self.args.result_dir,write_file),"w"))
        return model, FF

    
    def process(self,smiles):
        mod_pro = self.prepare_protein()
        mod_mol, off_mol = self.prepare_ligand(self.args.docked_file,smiles)
        mod_pro.add(mod_mol.topology,mod_mol.positions)
        model, FF = self.prep_model(mod_pro.topology,mod_pro.positions,off_mol)
        return model, FF


    def run(self,model,FF,run_index):
        result_dir = os.path.join(self.args.result_dir,run_index)
        platform = mm.Platform.getPlatformByName("CUDA" if self.args.cuda else "cpu")
        platformProperties = {"Precision":"single"}

        system = FF.createSystem(model.topology,
                                nonbondedMethod=self.args.nonbondedMethod,
                                nonbondedCutoff=self.args.nonbondedCutoff,
                                constraints=self.args.constraints,
                                rigidWater=self.args.rigidWater,
                                ewaldErrorTolerance=self.args.ewaldErrorTolerance)
        system.addForce(mm.MonteCarloBarostat(self.args.pressure,self.args.temperature,self.args.barostatInterval))

        integrator = mm.LangevinMiddleIntegrator(self.args.temperature,self.args.friction,self.args.dt)
        integrator.setConstraintTolerance(self.args.constraintTolerance)
        simulation = app.Simulation(model.topology,system,integrator,platform,platformProperties)
        simulation.context.setPositions(model.positions)

        with open(os.path.join(result_dir,"sysmtem.xml"),mode="w") as f:
            f.write(mm.XmlSerializer.serialize(system))
        with open(os.path.join(result_dir,"integrator.xml"),mode="w") as f:
            f.write(mm.XmlSerializer.serialize(integrator))

        print("{} simulation start".format(run_index))
        simulation.minimizeEnergy()
        with open(os.path.join(result_dir,"minimized.pdb"),"w") as f:
            app.PDBFile.writeFile(
                simulation.topology,
                simulation.context.getState(getPositions=True,enforcePeriodicBox=False).getPositions(),
                file=f,
                keepIds=True
            )
        ts = time.perf_counter()
        simulation.context.setVelocitiesToTemperature(self.args.temperature)
        simulation.step(self.args.equilibrationSteps)

        dcdReporter = app.DCDReporter(os.path.join(result_dir,"trajectory.dcd"),self.args.reporterStep,enforcePeriodicBox=True)
        dataReporter = app.StateDataReporter(os.path.join(result_dir,"log.txt"),self.args.reporterStep,totalSteps=self.args.steps,step=True,speed=True,progress=True,potentialEnergy=True,temperature=True,separator="\t")
        checkpointReporter = app.CheckpointReporter(os.path.join(result_dir,"checkpoint.chk"),self.args.checkpointStep)
        simulation.reporters.append(dcdReporter)
        simulation.reporters.append(dataReporter)
        simulation.reporters.append(checkpointReporter)
        simulation.currentStep = 0

        simulation.step(self.args.steps)
        tg = time.perf_counter()
        simulation.saveState(os.path.join(result_dir,"final_state.xml"))
        state = simulation.context.getState(getPositions=True,enforcePeriodicBox=True)
        with open(os.path.join(result_dir,"final_state.cif"),mode="w") as file:
            app.PDBxFile.writeFile(simulation.topology,state.getPositions(),file)
        print("{} simulation finish".format(run_index))
        elapse = tg - ts
        h = elapse // 3600
        m = (elapse % 3600) // 60
        s = elapse % 60
        print(f"equilibration - simulation elapsed time: {h} h {m} min {s} sec")


def main():
    args = get_args()
    protein = glob.glob("input/vorinostat_docked/*.pdb")
    smiles = "O=C(Nc1ccccc1)CCCCCCC(=O)NO"
    for v in protein:
        p = v.split("/")[-1].split(".")[0]
        result_dir = os.path.join(args.result_dir,p)
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        args.docked_file = v
        print("{} process start".format(p))
        mdp = MDPreparation(args)
        model, FF = mdp.process(smiles)
        mdp.run(model,FF,p)

if __name__ == "__main__":
    main()