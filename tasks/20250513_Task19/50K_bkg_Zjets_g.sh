#!/bin/bash

set -x

echo "50K background events, Z + jets from gluon, Z decay to charged leptons"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_ReclpTmin40.tcl /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/bkg_Zjets_g.cmnd  /gfsvol01/atlas/giuliac/condor/20250513_Task19/50K_bkg_Zjets_g.root > /gfsvol01/atlas/giuliac/condor/20250513_Task19/50K_bkg_Zjets_g.log
