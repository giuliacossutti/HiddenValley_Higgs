#!/bin/bash

set -x

echo "test of signal EJ showering from a PowHeg-generated sample"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_Unique.tcl /gfsvol01/atlas/giuliac/condor/20250616_Task23/pythia8_cards/test_showering_signal_EJ.cmnd /gfsvol01/atlas/giuliac/condor/20250616_Task23/root_outputs/test_showering_signal_EJ.root > /gfsvol01/atlas/giuliac/condor/20250616_Task23/log_outputs/test_showering_signal_EJ.log 
