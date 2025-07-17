#!/bin/bash

set -x

echo "5M background events, Z + jets from gluon, Z decay to charged leptons"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/condor/20250717_Task27/delphes_cards/delphes_card_ATLAS_Reclustering_Unique_Compact_eb_${1}.tcl /gfsvol01/atlas/giuliac/condor/20250717_Task27/pythia8_cards/bkg_Zjets_g_${1}.cmnd /eos/infnts/atlas/giuliac/condor/20250717_Task27/root_outputs/50K_bkg_Zjets_g_eb_${1}.root > /eos/infnts/atlas/giuliac/condor/20250717_Task27/log_outputs/50K_bkg_Zjets_g_eb_${1}.log 
