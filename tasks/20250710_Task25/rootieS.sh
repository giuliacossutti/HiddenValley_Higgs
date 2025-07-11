#!/bin/bash

set -x

cd /gfsvol01/atlas/giuliac/HiggsTutorial/MG5_aMC_v3_5_6/Delphes/

root -l -q /gfsvol01/atlas/giuliac/HiddenValley_Higgs/tasks/20250710_Task25/BuildHistos.C'("signal_EJ")'
