#!/bin/bash

pdb4amber -i /workspace/MD/data/test/protein/7zzs_ZINC.pdb -o 7zzs_ZINC_processed.pdb --reduce

#==================================================
#Summary of pdb4amber for: /workspace/MD/data/test/protein/7zzs_ZINC.pdb
#===================================================

#----------Chains
#The following (original) chains have been found:
#A

#---------- Alternate Locations (Original Residues!))

#The following residues had alternate locations:
#None
#-----------Non-standard-resnames


#---------- Missing heavy atom(s)

#None


antechamber -fi mol2 -i /workspace/MD/data/test/ligand/7zzs_vorinostat.mol2 \
    -fo mol2 -o /workspace/MD/experiments/250607_amber_testrun1/output/vorinostat_calculated.mol2

#Info: acdoctor mode is on: check and diagnose problems in the input file.
#Info: The atom type is set to gaff; the options available to the -at flag are
#      gaff, gaff2, amber, bcc, abcg2, and sybyl.

#-- Check Format for mol2 File --
#   Status: pass


tleap

#-I: Adding /opt/ambertools25/dat/leap/prep to search path.
#-I: Adding /opt/ambertools25/dat/leap/lib to search path.
#-I: Adding /opt/ambertools25/dat/leap/parm to search path.
#-I: Adding /opt/ambertools25/dat/leap/cmd to search path.

#Welcome to LEaP!
#(no leaprc in search path)