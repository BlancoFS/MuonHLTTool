#!/bin/bash

root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L2Muon")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L3IOFromL1")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L3Iter2FromL1")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L3Iter0FromL1")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L3OI")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L3MuonMerged")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L3MuonNoId")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "L3Muon")'

root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "hltL3fL1TkSingleMu22L3Filtered24Q")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70")'
root -l -b -q 'drawBDTEffFull_binary_origin.C("POG_v7", "DY PU 200", "Phase2_DYToLL_M50_PU200_Binary", "hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk")'
