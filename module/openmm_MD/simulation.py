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
from parmed import openmm as pmm

import numpy as np
import time
import os
import gc
from tqdm.notebook import tqdm
import pickle

def simulate(complex_model, FF, lig_name, args):
    """
    input: complex_model, FF, lig_name, args
    output(to directory): config + results (config, trajectory, structure)
    """
    nonbondedMethod = PME
    nonbondedCutoff = args.nonbondedCutoff * unit.nanometers
    constraints = HBonds
    temperature = args.temperature * unit.kelvin
    pressure = args.pressure * unit.atmospheres
    friction = args.friction / unit.picoseconds
    dt = args.dt * unit.femtoseconds
    platform = Platform.getPlatformByName(args.platform)
    platformProperties = {"Precision":args.platformProperties}

    # convert to amber format(.prmtop) for mmpbsa
    system = FF.createSystem(complex_model.topology,
                             nonbondedMethod=nonbondedMethod,
                             nonbondedCutoff=nonbondedCutoff,
                             constraints=constraints,
                             rigidWater=args.rigidWater,
                             ewaldErrorTolerance=args.ewaldErrorTolerance,
                             flexibleConstraints=True)
    system.addForce(MonteCarloBarostat(pressure, temperature, args.barostatInterval))
    struct = pmm.load_topology(complex_model.topology, system)
    struct.save(str(args.pbsa_input_dir_path / "complex_solv.prmtop"), overwrite=True)
    complex_structure = struct["!(:HOH,NA,CL)"]
    complex_structure.box = None
    receptor_structure = struct["!(:UNK,HOH,NA,CL)"]
    receptor_structure.box = None
    ligand_structure = struct[":UNK"]
    ligand_structure.box = None
    complex_structure.save(str(args.pbsa_input_dir_path / "complex.prmtop"), overwrite=True)
    receptor_structure.save(str(args.pbsa_input_dir_path / "receptor.prmtop"), overwrite=True)
    ligand_structure.save(str(args.pbsa_input_dir_path / "ligand.prmtop"), overwrite=True)
    del system, struct, complex_structure, receptor_structure, ligand_structure
    gc.collect()
    
    # Simulation
    system = FF.createSystem(complex_model.topology,
                             nonbondedMethod=nonbondedMethod,
                             nonbondedCutoff=nonbondedCutoff,
                             constraints=constraints,
                             rigidWater=args.rigidWater,
                             ewaldErrorTolerance=args.ewaldErrorTolerance)
    system.addForce(MonteCarloBarostat(pressure, temperature, args.barostatInterval))

    integrator = LangevinMiddleIntegrator(temperature, friction, dt)
    integrator.setConstraintTolerance(args.constraintTolerance)

    simulation = Simulation(complex_model.topology, system, integrator, platform, platformProperties)
    simulation.context.setPositions(complex_model.positions)

    with open(args.config_dir_path / f"{lig_name}_system.xml", "w") as f:
        f.write(XmlSerializer.serialize(system))
    with open(args.config_dir_path / f"{lig_name}_integrator.xml", "w") as f:
        f.write(XmlSerializer.serialize(integrator))
    
    os.makedirs(args.results_dir_path / f"{lig_name}")

    print("quilibration time: {} ps".format(np.floor(dt*(args.equilibrationSteps)/unit.picoseconds)))
    print("simulation time: {} ns".format(np.floor(dt*(args.steps)/unit.nanoseconds)))

    print(f"{lig_name} minimization start...")
    simulation.minimizeEnergy()
    with open(args.results_dir_path / f"{lig_name}" / f"{lig_name}_minimized.pdb", "w") as f:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(),
            file=f,
            keepIds=True
        )
    
    dcdreporter = DCDReporter(str(args.results_dir_path / f"{lig_name}" / f"{lig_name}_trajectory.dcd"), args.reporterStep, enforcePeriodicBox=True)
    datareporter = StateDataReporter(str(args.results_dir_path / f"{lig_name}" / f"{lig_name}_log.txt"), args.reporterStep, totalSteps=args.steps, step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator="\t")
    checkpointreporter = CheckpointReporter(str(args.results_dir_path / f"{lig_name}" / f"{lig_name}_checkpoint.chk"), args.checkpointStep)

    print(f"{lig_name} equilibration start...")
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(args.equilibrationSteps)

    with open(args.results_dir_path / f"{lig_name}" / f"{lig_name}_equilibrated.pdb", "w") as f:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(),
            file=f,
            keepIds=True
        )
    
    simulation.reporters.append(dcdreporter)
    simulation.reporters.append(datareporter)
    simulation.reporters.append(checkpointreporter)
    simulation.currentStep = 0

    print(f"{lig_name} simulation start...")
    simulation.step(args.steps)

    simulation.saveState(str(args.results_dir_path / f"{lig_name}" / f"{lig_name}_final_state.xml"))
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(args.results_dir_path / f"{lig_name}" / f"{lig_name}_final_state.cif", "w") as f:
        PDBxFile.writeFile(simulation.topology, state.getPositions(), f)
    
    print(f"{lig_name} simualtion finished")