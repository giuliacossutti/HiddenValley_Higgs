#!/bin/bash

set -x

echo "100 test background events, Z + jets from quark, Z decay to charged leptons"

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

DelphesPythia8 /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/delphes_card_ATLAS_Reclustering_Unique.tcl /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/test_xsec_Zjets_q.cmnd  /gfsvol01/atlas/giuliac/condor/20250520_Task20/100_test_xsec_Zjets_q.root > /gfsvol01/atlas/giuliac/condor/20250520_Task20/100_test_xsec_Zjets_q.log
