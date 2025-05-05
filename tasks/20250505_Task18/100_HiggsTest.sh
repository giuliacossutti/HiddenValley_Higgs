#!/bin/bash

set -x

echo "100 events to test Higgs with HiddenValley"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering.tcl /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/Higgs_EJ.cmnd  /gfsvol01/atlas/giuliac/condor/20250505_Task18/100_HiggsTest.root > /gfsvol01/atlas/giuliac/condor/20250505_Task18/100_HiggsTest.log
