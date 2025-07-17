#!/bin/bash

for seed in $(seq 1 10)
do
	declare -i myseed=$((129 + $seed))
	echo $myseed

	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_signal_EJ.cmnd > /gfsvol01/atlas/giuliac/condor/20250717_Task27/pythia8_cards/signal_EJ_$seed.cmnd
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_ZZ.cmnd > /gfsvol01/atlas/giuliac/condor/20250717_Task27/pythia8_cards/bkg_ZZ_$seed.cmnd
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_HZ_SM.cmnd > /gfsvol01/atlas/giuliac/condor/20250717_Task27/pythia8_cards/bkg_HZ_SM_$seed.cmnd
done

for seed in $(seq 1 100)
do
	declare -i myseed=$((129 + $seed))
	echo $myseed

	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_Zjets_q_pTHatMin.cmnd > /gfsvol01/atlas/giuliac/condor/20250717_Task27/pythia8_cards/bkg_Zjets_q_$seed.cmnd
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_Zjets_g_pTHatMin.cmnd > /gfsvol01/atlas/giuliac/condor/20250717_Task27/pythia8_cards/bkg_Zjets_g_$seed.cmnd

	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/condor_delphes_card_ATLAS_Reclustering_Unique_Compact.tcl > /gfsvol01/atlas/giuliac/condor/20250717_Task27/delphes_cards/delphes_card_ATLAS_Reclustering_Unique_Compact_$seed.tcl
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/condor_delphes_card_ATLAS_Reclustering_Unique_Compact_JES.tcl > /gfsvol01/atlas/giuliac/condor/20250717_Task27/delphes_cards/delphes_card_ATLAS_Reclustering_Unique_Compact_JES_$seed.tcl
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/condor_delphes_card_ATLAS_Reclustering_Unique_Compact_etrk.tcl > /gfsvol01/atlas/giuliac/condor/20250717_Task27/delphes_cards/delphes_card_ATLAS_Reclustering_Unique_Compact_etrk_$seed.tcl
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/delphes_cards/condor_delphes_card_ATLAS_Reclustering_Unique_Compact_eb.tcl > /gfsvol01/atlas/giuliac/condor/20250717_Task27/delphes_cards/delphes_card_ATLAS_Reclustering_Unique_Compact_eb_$seed.tcl
done
