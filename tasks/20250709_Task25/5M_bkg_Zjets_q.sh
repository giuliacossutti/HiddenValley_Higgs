#!/bin/bash

set -x

echo "5M background events, Z + jets from quark, Z decay to charged leptons"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

#DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_Unique.tcl /gfsvol01/atlas/giuliac/condor/20250709_Task25/pythia8_cards/bkg_Zjets_q_${1}.cmnd /eos/infnts/atlas/giuliac/condor/20250709_Task25/root_outputs/50K_bkg_Zjets_q_${1}.root > /eos/infnts/atlas/giuliac/condor/20250709_Task25/log_outputs/50K_bkg_Zjets_q_${1}.log 
#DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_Unique_Compact.tcl /gfsvol01/atlas/giuliac/condor/20250709_Task25/pythia8_cards/bkg_Zjets_q_${1}.cmnd /eos/infnts/atlas/giuliac/condor/20250709_Task25/root_outputs/test50K_bkg_Zjets_q_${1}.root > /eos/infnts/atlas/giuliac/condor/20250709_Task25/log_outputs/test50K_bkg_Zjets_q_${1}.log 

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_Unique_Compact.tcl /gfsvol01/atlas/giuliac/condor/20250709_Task25/pythia8_cards/bkg_Zjets_q_${1}.cmnd /eos/infnts/atlas/giuliac/condor/20250709_Task25/root_outputs/50K_bkg_Zjets_q_${1}.root > /eos/infnts/atlas/giuliac/condor/20250709_Task25/log_outputs/50K_bkg_Zjets_q_${1}.log 
