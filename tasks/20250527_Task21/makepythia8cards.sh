#!/bin/bash

for seed in $(seq 11 50)
do
	declare -i myseed=$((129 + $seed))
	echo $myseed
#	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_signal_EJ.cmnd > /gfsvol01/atlas/giuliac/condor/20250527_Task21/pythia8_cards/signal_EJ_$seed.cmnd
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_Zjets_q.cmnd > /gfsvol01/atlas/giuliac/condor/20250527_Task21/pythia8_cards/bkg_Zjets_q_$seed.cmnd
	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_Zjets_g.cmnd > /gfsvol01/atlas/giuliac/condor/20250527_Task21/pythia8_cards/bkg_Zjets_g_$seed.cmnd
#	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_ZZ.cmnd > /gfsvol01/atlas/giuliac/condor/20250527_Task21/pythia8_cards/bkg_ZZ_$seed.cmnd
#	sed s/MYSEED/$myseed/  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/condor_bkg_HZ_SM.cmnd > /gfsvol01/atlas/giuliac/condor/20250527_Task21/pythia8_cards/bkg_HZ_SM_$seed.cmnd
done
