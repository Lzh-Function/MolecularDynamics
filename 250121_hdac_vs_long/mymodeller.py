from __future__ import division
from __future__ import absolute_import

__author__ = "Peter Eastman"
__version__ = "1.0"

from openmm import app
from openmm.app import Topology, PDBFile, ForceField
from openmm.app.forcefield import AllBonds, CutoffNonPeriodic, CutoffPeriodic, DrudeGenerator, _getDataDirectories
from openmm.app.internal import compiled
from openmm.vec3 import Vec3
from openmm import System, Context, NonbondedForce, CustomNonbondedForce, HarmonicBondForce, HarmonicAngleForce, VerletIntegrator, LangevinIntegrator, LocalEnergyMinimizer
from openmm.unit import nanometer, molar, elementary_charge, degree, acos, is_quantity, dot, norm, kilojoules_per_mole
import openmm.unit as unit
import openmm.app.element as elem
import gc
import os
import random
import sys
import xml.etree.ElementTree as etree
from copy import deepcopy
from math import ceil, floor, sqrt
from collections import defaultdict

class MyModeller(app.Modeller):
    def __init__(self,topology,positions):
        super().__init__(topology,positions)

    def addSolvent(self,forcefield,model="tip3p",boxSize=None,boxVectors=None,padding=None,numAdded=None,boxShape='cube',positiveIon='Na+',negativeIon='Cl-',ionicStrength=0*molar,neutralize=True,residueTemplates=dict()):
        if len([x for x in (boxSize, boxVectors, padding, numAdded) if x is not None]) > 1:
            raise ValueError('At most one of the following arguments may be specified: boxSize, boxVectors, padding, numAdded')

        # Load the pre-equilibrated water box.

        vdwRadiusPerSigma = 0.5612310241546864907
        if model == 'tip3p':
            waterRadius = 0.31507524065751241*vdwRadiusPerSigma
        elif model == 'spce':
            waterRadius = 0.31657195050398818*vdwRadiusPerSigma
        elif model == 'tip4pew':
            waterRadius = 0.315365*vdwRadiusPerSigma
        elif model == 'tip5p':
            waterRadius = 0.312*vdwRadiusPerSigma
        elif model == 'swm4ndp':
            waterRadius = 0.318395*vdwRadiusPerSigma
        else:
            raise ValueError('Unknown water model: %s' % model)
        pdb = PDBFile(os.path.join(os.path.dirname(__file__), model+'.pdb'))
        pdbTopology = pdb.getTopology()
        pdbPositions = pdb.getPositions().value_in_unit(nanometer)
        pdbResidues = list(pdbTopology.residues())
        pdbBoxSize = pdbTopology.getUnitCellDimensions().value_in_unit(nanometer)

        # Pick a unit cell size.

        if numAdded is not None:
            # Select a padding distance which is guaranteed to give more than the specified number of molecules.

            padding = 2.2*(numAdded/((len(pdbResidues)/pdbBoxSize[0]**3)*8))**(1.0/3.0)
            if padding < 0.5:
                padding = 0.5 # Ensure we have enough when adding very small numbers of molecules
        if boxVectors is not None:
            if is_quantity(boxVectors[0]):
                boxVectors = (boxVectors[0].value_in_unit(nanometer), boxVectors[1].value_in_unit(nanometer), boxVectors[2].value_in_unit(nanometer))
            box = Vec3(boxVectors[0][0], boxVectors[1][1], boxVectors[2][2])
            vectors = boxVectors
        elif boxSize is not None:
            if is_quantity(boxSize):
                boxSize = boxSize.value_in_unit(nanometer)
            box = Vec3(boxSize[0], boxSize[1], boxSize[2])
            vectors = (Vec3(boxSize[0], 0, 0), Vec3(0, boxSize[1], 0), Vec3(0, 0, boxSize[2]))
        elif padding is not None:
            if is_quantity(padding):
                padding = padding.value_in_unit(nanometer)
            if len(self.positions) == 0:
                radius = 0
            else:
                positions = self.positions.value_in_unit(nanometer)
                minRange = Vec3(*(min((pos[i] for pos in positions)) for i in range(3)))
                maxRange = Vec3(*(max((pos[i] for pos in positions)) for i in range(3)))
                center = 0.5*(minRange+maxRange)
                radius = max(unit.norm(center-pos) for pos in positions)
            width = max(2*radius+padding, 2*padding)
            vectors = self._computeBoxVectors(width, boxShape)
            box = Vec3(vectors[0][0], vectors[1][1], vectors[2][2])
        else:
            box = self.topology.getUnitCellDimensions().value_in_unit(nanometer)
            vectors = self.topology.getPeriodicBoxVectors().value_in_unit(nanometer)
            if box is None:
                raise ValueError('Neither the box size, box vectors, nor padding was specified, and the Topology does not define unit cell dimensions')

        # Have the ForceField build a System for the solute from which we can determine van der Waals radii.

        system = forcefield.createSystem(self.topology, residueTemplates=residueTemplates)
        nonbonded = None
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), NonbondedForce):
                nonbonded = system.getForce(i)
        if nonbonded is None:
            raise ValueError('The ForceField does not specify a NonbondedForce')
        cutoff = [waterRadius]*system.getNumParticles()
        for i in range(system.getNumParticles()):
            params = nonbonded.getParticleParameters(i)
            if params[2] != 0*kilojoules_per_mole:
                cutoff[i] += params[1].value_in_unit(nanometer)*vdwRadiusPerSigma
        waterCutoff = waterRadius
        if len(cutoff) == 0:
            maxCutoff = waterCutoff
        else:
            maxCutoff = max(waterCutoff, max(cutoff))

        # Copy the solute over.

        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(vectors*nanometer)
        newAtoms = {}
        newPositions = []*nanometer
        newResidueTemplates=dict()
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                if residue in residueTemplates:
                    newResidueTemplates[newResidue] = residueTemplates[residue]
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        # Sort the solute atoms into cells for fast lookup.

        if len(self.positions) == 0:
            positions = []
        else:
            positions = deepcopy(self.positions.value_in_unit(nanometer))
        cells = app.modeller._CellList(positions, maxCutoff, vectors, True)

        # Create a function to compute the distance between two points, taking periodic boundary conditions into account.

        periodicDistance = compiled.periodicDistance(vectors)

        # Find the list of water molecules to add.

        newChain = newTopology.addChain()
        if len(positions) == 0:
            center = Vec3(0, 0, 0)
        else:
            center = [(max((pos[i] for pos in positions))+min((pos[i] for pos in positions)))/2 for i in range(3)]
            center = Vec3(center[0], center[1], center[2])
        numBoxes = [int(ceil(box[i]/pdbBoxSize[i])) for i in range(3)]
        addedWaters = []
        for boxx in range(numBoxes[0]):
            for boxy in range(numBoxes[1]):
                for boxz in range(numBoxes[2]):
                    offset = Vec3(boxx*pdbBoxSize[0], boxy*pdbBoxSize[1], boxz*pdbBoxSize[2])
                    for residue in pdbResidues:
                        oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
                        atomPos = pdbPositions[oxygen.index]+offset
                        if not any((atomPos[i] > box[i] for i in range(3))):
                            # This molecule is inside the box, so see how close to it is to the solute.

                            atomPos += center-box/2
                            for i in cells.neighbors(atomPos):
                                if periodicDistance(atomPos, positions[i]) < cutoff[i]:
                                    break
                            else:
                                # Record this water molecule as one to add.

                                addedWaters.append((residue.index, atomPos))
        

        if numAdded is not None:
            # We added many more waters than we actually want.  Sort them based on distance to the nearest box edge and
            # only keep the ones in the middle.

            lowerBound = center-box/2
            upperBound = center+box/2
            distToEdge = (min(min(pos-lowerBound), min(upperBound-pos)) for index, pos in addedWaters)
            sortedIndex = [i[0] for i in sorted(enumerate(distToEdge), key=lambda x: -x[1])]
            addedWaters = [addedWaters[i] for i in sortedIndex[:numAdded]]

            # Compute a new periodic box size.

            maxSize = max(max((pos[i] for index, pos in addedWaters))-min((pos[i] for index, pos in addedWaters)) for i in range(3))
            maxSize += 0.1  # Add padding to reduce clashes at the edge.
            newTopology.setPeriodicBoxVectors(self._computeBoxVectors(maxSize, boxShape))
        else:
            # There could be clashes between water molecules at the box edges.  Find ones to remove.

            upperCutoff = center+box/2-Vec3(waterCutoff, waterCutoff, waterCutoff)
            lowerCutoff = center-box/2+Vec3(waterCutoff, waterCutoff, waterCutoff)
            lowerSkinPositions = [pos for index, pos in addedWaters if pos[0] < lowerCutoff[0] or pos[1] < lowerCutoff[1] or pos[2] < lowerCutoff[2]]
            filteredWaters = []
            cells.cells = {}
            for i in range(len(lowerSkinPositions)):
                cell = cells.cellForPosition(lowerSkinPositions[i])
                if cell in cells.cells:
                    cells.cells[cell].append(i)
                else:
                    cells.cells[cell] = [i]
            for entry in addedWaters:
                pos = entry[1]
                if pos[0] < upperCutoff[0] and pos[1] < upperCutoff[1] and pos[2] < upperCutoff[2]:
                    filteredWaters.append(entry)
                else:
                    if not any((periodicDistance(lowerSkinPositions[i], pos) < waterCutoff and norm(lowerSkinPositions[i]-pos) > waterCutoff for i in cells.neighbors(pos))):
                        filteredWaters.append(entry)
            addedWaters = filteredWaters
        
        # Add the water molecules.
        waterPos = {}
        for index, pos in addedWaters:
            newResidue = newTopology.addResidue(residue.name, newChain)
            residue = pdbResidues[index]
            oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
            oPos = pdbPositions[oxygen.index]
            molAtoms = []
            for atom in residue.atoms():
                molAtoms.append(newTopology.addAtom(atom.name, atom.element, newResidue))
                newPositions.append((pos+pdbPositions[atom.index]-oPos)*nanometer)
                if atom.element == elem.oxygen:
                    waterPos[newResidue] = newPositions[-1]
            for atom1 in molAtoms:
                if atom1.element == elem.oxygen:
                    for atom2 in molAtoms:
                        if atom2.element == elem.hydrogen:
                            newTopology.addBond(atom1, atom2)

        self.topology = newTopology
        self.positions = newPositions
        # Total number of waters in the box
        numTotalWaters = len(waterPos)

        # Add ions to neutralize the system.
        self._addIons(forcefield, numTotalWaters, waterPos, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength, neutralize=neutralize, residueTemplates=newResidueTemplates)

    def _addIons(self,forcefield,numWaters,replaceableMols,ionCutoff=0.05*nanometer,positiveIon='Na+',negativeIon='Cl-',ionicStrength=0*molar,neutralize=True,residueTemplates=dict()):
        posIonElements = {'Cs+': elem.cesium, 'K+': elem.potassium,
                          'Li+': elem.lithium, 'Na+': elem.sodium,
                          'Rb+': elem.rubidium}
        negIonElements = {'Cl-': elem.chlorine, 'Br-': elem.bromine,
                          'F-': elem.fluorine, 'I-': elem.iodine}

        ionPositions = []

        numReplaceableMols = len(replaceableMols)

        # Fetch ion elements from user input
        if positiveIon not in posIonElements:
            raise ValueError('Illegal value for positive ion: {}'.format(positiveIon))
        if negativeIon not in negIonElements:
            raise ValueError('Illegal value for negative ion: {}'.format(negativeIon))
        positiveElement = posIonElements[positiveIon]
        negativeElement = negIonElements[negativeIon]

        # Determine the total charge of the system
        system = forcefield.createSystem(self.topology, residueTemplates=residueTemplates)
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), NonbondedForce):
                nonbonded = system.getForce(i)
                break
        else:
            raise ValueError('The ForceField does not specify a NonbondedForce')

        totalCharge = 0.0
        for i in range(nonbonded.getNumParticles()):
            nb_i = nonbonded.getParticleParameters(i)
            totalCharge += nb_i[0].value_in_unit(elementary_charge)
        # Round to nearest integer
        totalCharge = int(floor(0.5 + totalCharge))

        # Figure out how many ions to add based on requested params/concentration
        numPositive, numNegative = 0, 0
        if neutralize:
            if abs(totalCharge) > numReplaceableMols:
                raise Exception('Cannot neutralize the system because the charge is greater than the number of available positions for ions')
            if totalCharge > 0:
                numNegative += totalCharge
            else:
                numPositive -= totalCharge

        if ionicStrength > 0 * molar:
            numIons = (numWaters - numPositive - numNegative) * ionicStrength / (55.4 * molar)  # Pure water is about 55.4 molar (depending on temperature)
            numPairs = int(floor(numIons + 0.5))
            numPositive += numPairs
            numNegative += numPairs
        totalIons = numPositive + numNegative

        if totalIons > 0:
            # Randomly select a set of waters
            # while ensuring ions are not placed too close to each other.
            modeller = MyModeller(self.topology, self.positions)

            replaceableList = list(replaceableMols.keys())
            numAddedIons = 0
            numTrials = 10  # Attempts to add ions N times before quitting
            toReplace = []  # list of molecules to be replaced
            while numAddedIons < totalIons:
                pickedMol = random.choice(replaceableList)
                replaceableList.remove(pickedMol)
                # Check distance to other ions
                for pos in ionPositions:
                    distance = norm(pos - replaceableMols[pickedMol])
                    if distance <= ionCutoff:
                        numTrials -= 1
                        break
                else:
                    toReplace.append(pickedMol)
                    ionPositions.append(replaceableMols[pickedMol])
                    numAddedIons += 1

                    n_trials = 10

                if n_trials == 0:
                    raise ValueError('Could not add more than {} ions to the system'.format(numAddedIons))
            
            # Replace waters/ions in the topology
            modeller.delete(toReplace)
            ionChain = modeller.topology.addChain()
            for i, water in enumerate(toReplace):
                element = (positiveElement if i < numPositive else negativeElement)
                newResidue = modeller.topology.addResidue(element.symbol.upper(), ionChain)
                ioncharge = 1 if i < numPositive else -1
                modeller.topology.addAtom(element.symbol, element, newResidue, formalCharge=ioncharge)
                modeller.positions.append(replaceableMols[water])

            # Update topology/positions
            self.topology = modeller.topology
            self.positions = modeller.positions
    
    def delete(self, toDelete):
        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
        newAtoms = {}
        newPositions = []*nanometer
        deleteSet = set(toDelete)
        for chain in self.topology.chains():
            if chain not in deleteSet:
                needNewChain = True;
                for residue in chain.residues():
                    if residue not in deleteSet:
                        needNewResidue = True
                        for atom in residue.atoms():
                            if atom not in deleteSet:
                                if needNewChain:
                                    newChain = newTopology.addChain(chain.id)
                                    needNewChain = False;
                                if needNewResidue:
                                    newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                                    needNewResidue = False;
                                newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                                newAtoms[atom] = newAtom
                                newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                if bond not in deleteSet and (bond[1], bond[0]) not in deleteSet:
                    newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)
        self.topology = newTopology
        self.positions = newPositions
