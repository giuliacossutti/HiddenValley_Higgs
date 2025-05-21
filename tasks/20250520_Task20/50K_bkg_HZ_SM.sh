#!/bin/bash

set -x

echo "50K background events, Standard Model HZ -> b bbar l+ l-"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_Unique.tcl /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/bkg_HZ_SM.cmnd  /gfsvol01/atlas/giuliac/condor/20250520_Task20/50K_bkg_HZ_SM.root > /gfsvol01/atlas/giuliac/condor/20250520_Task20/50K_bkg_HZ_SM.log
