#!/bin/bash

export mylhef=/eos/infnts/atlas/giuliac/atlasmcsamples/TXT.33690942._005001.events

sed s@MYLHEF@$mylhef@  /gfsvol01/atlas/giuliac/HiddenValley_Higgs/pythia8_cards/showering_signal_EJ.cmnd > /gfsvol01/atlas/giuliac/condor/20250616_Task23/pythia8_cards/test_showering_signal_EJ.cmnd
