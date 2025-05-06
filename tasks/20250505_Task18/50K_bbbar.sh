#!/bin/bash

set -x

echo "50K events to test Higgs with HiddenValley, with pTmin of reclustered jets = 40 GeV, dark matter decay to b bbar"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_ReclpTmin40.tcl /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/Higgs_EJ_bbbar.cmnd  /gfsvol01/atlas/giuliac/condor/20250505_Task18/50K_bbbar.root > /gfsvol01/atlas/giuliac/condor/20250505_Task18/50K_bbbar.log
