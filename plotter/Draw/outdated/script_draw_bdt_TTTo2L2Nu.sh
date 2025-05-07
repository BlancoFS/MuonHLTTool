#!/bin/bash

#VER="v02"
VER="v14"

# root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3Iter2FromL1")'
# root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3IOFromL1")'
# root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3MuonNoId")'
# root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3Muon")'

root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Iter2FromL1")'
root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3IOFromL1")'
root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3MuonNoId")'
root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Muon")'

root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Iter2FromL1", true)'
root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3IOFromL1", true)'
root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3MuonNoId", true)'
root -l -b -q 'drawBDTEff.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Muon", true)'

# root -l -b -q 'drawBDTEff.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3Iter2FromL1")'
# root -l -b -q 'drawBDTEff.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3IOFromL1")'
# root -l -b -q 'drawBDTEff.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3MuonNoId")'
# root -l -b -q 'drawBDTEff.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3Muon")'


# root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3Iter2FromL1")'
# root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3IOFromL1")'
# root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3MuonNoId")'
# root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3Muon")'

root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Iter2FromL1")'
root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3IOFromL1")'
root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3MuonNoId")'
root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Muon")'

root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Iter2FromL1", true)'
root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3IOFromL1", true)'
root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3MuonNoId", true)'
root -l -b -q 'drawBDTEffFull.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3Muon", true)'
# root -l -b -q 'drawBDTEffFull.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3Iter2FromL1")'
# root -l -b -q 'drawBDTEffFull.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3IOFromL1")'
# root -l -b -q 'drawBDTEffFull.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3MuonNoId")'
# root -l -b -q 'drawBDTEffFull.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3Muon")'




# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt0", false)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt8", false)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt24", false)'

root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3pt8", false)'
root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3pt24", false)'

# root -l -b -q 'drawBDTFrac.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3pt0", false)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3pt8", false)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3pt24", false)'

# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt0", true)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt8", true)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "PU200-DYToLL_M50_TDRBDT", "L3pt24", true)'

root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3pt8", true)'
root -l -b -q 'drawBDTFrac.C("'$VER'", "DY PU 200", "Phase2_PU200_TTTo2L2Nu", "L3pt24", true)'

# root -l -b -q 'drawBDTFrac.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3pt0", true)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3pt8", true)'
# root -l -b -q 'drawBDTFrac.C("'$VER'", "t#bar{t} PU 200", "PU200-TTToSemiLep_TDRBDT", "L3pt24", true)'

