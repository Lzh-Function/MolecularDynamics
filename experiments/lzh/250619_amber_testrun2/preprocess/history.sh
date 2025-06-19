#!/bin/bash

pdb4amber -i 7zzs_ZINC.pdb -o hdac2_zinc_processed.pdb --reduce

antechamber -fi mol2 -i 7zzs_vorinostat.mol2 -fo mol2 -o vorinostat_processed.mol2 -at gaff2

# TODO modeling HDAC2 zinc ion