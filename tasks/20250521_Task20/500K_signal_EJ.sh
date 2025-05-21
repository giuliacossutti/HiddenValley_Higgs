#!/bin/bash

set -x

echo "500K Emerging Jet signal events, Z decay to charged leptons, DM decay to b bbar"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_Unique.tcl /gfsvol01/atlas/giuliac/condor/20250521_Task20/pythia8_cards/signal_EJ_${1}.cmnd /gfsvol01/atlas/giuliac/condor/20250521_Task20/root_outputs/50K_signal_EJ_${1}.root > /gfsvol01/atlas/giuliac/condor/20250521_Task20/log_outputs/50K_signal_EJ_${1}.log 
