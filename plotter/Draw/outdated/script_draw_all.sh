#!/bin/bash

###################
# -- TDR plots -- #
###################

VER="v1"
root -l -b -q 'drawPurityID.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "", false)'
# root -l -b -q 'drawPurityID.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep", "", false)'
root -l -b -q 'drawEfficiencyIso.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# root -l -b -q 'drawEfficiencyIso.C("'$VER'", "#scale[0.8]{QCD (#mu enriched) PU200}", "PU200-QCD_MuEn", "Iso / any L3")'

root -l -b -q 'drawEfficiencyIterL3L1.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
root -l -b -q 'drawEfficiencyIterL3L2.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "L1Tk")'
root -l -b -q 'drawEfficiencyL3.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "L1Tk")'
root -l -b -q 'drawEfficiencyID.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "L1Tk")'
root -l -b -q 'drawEfficiencyMu50.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
root -l -b -q 'drawEfficiencyIsoMu24.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
root -l -b -q 'drawEfficiencyDoubleMuon.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# root -l -b -q 'drawEfficiencyDoubleMuonRegions.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
root -l -b -q 'drawEfficiencyDoubleMuonRegions2D.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'

# VER="v2"
# root -l -b -q 'drawPurityID.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "", false)'
# # root -l -b -q 'drawPurityID.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep", "", false)'
# root -l -b -q 'drawEfficiencyIso.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# # root -l -b -q 'drawEfficiencyIso.C("'$VER'", "#scale[0.8]{QCD (#mu enriched) PU200}", "PU200-QCD_MuEn", "Iso / any L3")'

# root -l -b -q 'drawEfficiencyIterL3L1.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# root -l -b -q 'drawEfficiencyIterL3L2.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "L1Tk")'
# root -l -b -q 'drawEfficiencyL3.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "L1Tk")'
# root -l -b -q 'drawEfficiencyID.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "L1Tk")'
# root -l -b -q 'drawEfficiencyMu50.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# root -l -b -q 'drawEfficiencyIsoMu24.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# root -l -b -q 'drawEfficiencyDoubleMuon.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# # root -l -b -q 'drawEfficiencyDoubleMuonRegions.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'
# root -l -b -q 'drawEfficiencyDoubleMuonRegions2D.C("'$VER'", "DY PU 200", "Phase2_DYToLL_M50_PU200_default", "")'



# root -l -b -q 'drawRes.C("'$VER'", "qbpt", true)' >Res-$VER-qbpt.log
# root -l -b -q 'drawResSingleBinComp.C("'$VER'", "qbpt", true)'

# VER="v90"
# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 0, true )'  >Rate-v90-PU200-Stitching-SMu.log
# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 0, false )' >Rate-v90-PU200-Stitching-else.log

# VER="v91"
# root -l -b -q 'drawL1RateWeight.C( "'$VER'", "200", true, true )' >L1Rate-PU200-MB-offine.log
# root -l -b -q 'drawL1RateWeight.C( "'$VER'", "200", true, false )' >L1Rate-PU200-MB-online.log
# root -l -b -q 'drawL1RateWeight.C( "'$VER'", "200", false, true )' >L1Rate-PU200-Stitching-offine.log
# root -l -b -q 'drawL1RateWeight.C( "'$VER'", "200", false, false )' >L1Rate-PU200-Stitching-online.log

# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 0, true )'  >Rate-PU200-Stitching-SMu.log
# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 0, false )' >Rate-PU200-Stitching-else.log
# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 1, true )'  >Rate-PU200-QCDMuEn-SMu.log
# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 1, false )' >Rate-PU200-QCDMuEn-else.log
# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 2, true )'  >Rate-PU200-MB-SMu.log
# root -l -b -q 'drawRateWeight.C( "'$VER'", "200", 2, false )' >Rate-PU200-MB-else.log


# VER="v02"
# root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3Muon")'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt8", true)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt24", true)'

# TDRBDT timing
# root -l -b -q 'DRAW_DQM_Module.C("average", "t#bar{t}", "PU200_0p035x0p02")'
# root -l -b -q 'DRAW_DQM_Module.C("average", "t#bar{t}", "PU200_0p07x0p04")'


