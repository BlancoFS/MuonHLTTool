#!/usr/bin/env python3
import sys
import os
import time
import glob
import subprocess
from shutil import copyfile
import gc

def jobSpritting( paths, nfiles ):
    str_dcap  = "dcap://cluster142.knu.ac.kr/"
    out = []

    lines = glob.glob(paths+"/ntuple*.root")
    lines.sort(key=os.path.getmtime)

    n = 0
    files = ''
    for i, line in enumerate(lines):
        if 'ntuple' not in line:
            continue
        if '.root' not in line:
            continue

        files += '\\"%s%s\\",' % (str_dcap, line)

        if i == nfiles*(n+1) - 1:
            filesout = '\'{'+files+'}\''
            filesout = filesout.replace(',}', '}')
            out.append( ( n, filesout) )
            n = n+1
            files = ''

        if i == len(lines)-1 and files != '':
            filesout = '\'{'+files+'}\''
            filesout = filesout.replace(',}', '}')
            out.append( ( n, filesout) )

    return out


if __name__ == '__main__':
    DCAP      = "dcap://cluster142.knu.ac.kr/"
    BASE      = "/pnfs/knu.ac.kr/data/cms/store/user/moh/PhaseII/20210512/"
    crabVER   = "_L3FromL1TkMuon_20210512"
    histoVER  = "v90"

    doHadd = False
    if len(sys.argv) > 1 and 'hadd' == sys.argv[1]:
        doHadd = True

    doRecover = False
    if len(sys.argv) > 1 and 'recover' == sys.argv[1]:
        doRecover = True

    samples = [

        # -- PU 200
        ( ("Eff", "EffFR", "EffDM", "Res"), "PU200-DYToLL_M50", "/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU200-DYToLL_M10to50", "/DYJetsToLL_M-10to50_TuneCP5_14TeV-madgraphMLM-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-WToLNu", "/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),

        ( ("Eff", "EffFR"), "PU200-TTToSemiLep", "/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU200-TTTo2L2Nu", "/TTTo2L2Nu_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU200-TT", "/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v2/FEVT"),

        # -- PU 140
        ( ("Eff", "EffFR", "EffDM", "Res"), "PU140-DYToLL_M50", "/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_pilot_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU140-DYToLL_M10to50", "/DYToLL_M-10To50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_pilot_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-WToLNu", "/WJetsToLNu_TuneCP5_14TeV-amcatnloFXFX-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),

        ( ("Eff", "EffFR"), "PU140-TTToSemiLep", "/TTToSemiLepton_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU140-TTTo2L2Nu", "/TTTo2L2Nu_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU140-TT", "/TT_TuneCP5_14TeV-powheg-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"),

        # -- QCDs
        # ( (), "PU140-MinBias", "/MinBias_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU140-MinBias-NewMB", "/MinBias_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT"),
        # ( (), "PU200-MinBias", "/MinBias_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU200-MinBias-NewMB", "/MinBias_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/FEVT"),

        ( ("Eff", "EffFR"), "PU140-QCD_Pt120to170_MuEn", "/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU200-QCD_Pt120to170_MuEn", "/QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU140-QCD_Pt15to20_MuEn", "/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU200-QCD_Pt15to20_MuEn", "/QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU140-QCD_Pt170to300_MuEn", "/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU200-QCD_Pt170to300_MuEn", "/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU140-QCD_Pt20to30_MuEn", "/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU200-QCD_Pt20to30_MuEn", "/QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt300toInf_MuEn", "/QCD_Pt-300toInf_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt300toInf_MuEn", "/QCD_Pt-300toInf_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU140-QCD_Pt30to50_MuEn", "/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU200-QCD_Pt30to50_MuEn", "/QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU140-QCD_Pt50to80_MuEn", "/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU200-QCD_Pt50to80_MuEn", "/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU140-QCD_Pt80to120_MuEn", "/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        ( ("Eff", "EffFR"), "PU200-QCD_Pt80to120_MuEn", "/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),

        # ( (), "PU140-QCD_Pt120to170", "/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt120to170", "/QCD_Pt_120to170_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt170to300", "/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt170to300", "/QCD_Pt_170to300_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt20to30-NewMB", "/QCD_Pt_20to30_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt20to30-NewMB", "/QCD_Pt_20to30_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1-v2/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt300to470", "/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt300to470", "/QCD_Pt_300to470_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt30to50", "/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt30to50-NewMB", "/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt30to50", "/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt30to50-NewMB", "/QCD_Pt_30to50_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt470to600", "/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt470to600", "/QCD_Pt_470to600_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt50to80", "/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt50to80-NewMB", "/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v2/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt50to80", "/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt50to80-NewMB", "/QCD_Pt_50to80_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_withNewMB_111X_mcRun4_realistic_T15_v1_ext1-v3/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt600oInf", "/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt600oInf", "/QCD_Pt_600oInf_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt80to120", "/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt80to120", "/QCD_Pt_80to120_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),

        # ( (), "PU140-QCD_Pt120to170_EMEn", "/QCD_Pt-120to170_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt120to170_EMEn", "/QCD_Pt-120to170_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt15to20_EMEn", "/QCD_Pt-15to20_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt15to20_EMEn", "/QCD_Pt-15to20_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt170to300_EMEn", "/QCD_Pt-170to300_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt170to300_EMEn", "/QCD_Pt-170to300_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt20to30_EMEn", "/QCD_Pt-20to30_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt20to30_EMEn", "/QCD_Pt-20to30_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt300toInf_EMEn", "/QCD_Pt-300toInf_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt300toInf_EMEn", "/QCD_Pt-300toInf_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt30to50_EMEn", "/QCD_Pt-30to50_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt30to50_EMEn", "/QCD_Pt-30to50_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt50to80_EMEn", "/QCD_Pt-50to80_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt50to80_EMEn", "/QCD_Pt-50to80_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-QCD_Pt80to120_EMEn", "/QCD_Pt-80to120_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU200-QCD_Pt80to120_EMEn", "/QCD_Pt-80to120_EMEnriched_TuneCP5_14TeV_pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),

        # Others
        # ( (), "PU140-Zprime_M6000", "/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( ("Eff", "EffDM"), "PU200-Zprime_M6000", "/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( (), "PU140-MuonGun", "/DoubleMuon_gun_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( ("Eff"), "PU200-MuonGun", "/DoubleMuon_gun_FlatPt-1To100/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1-v1/FEVT"),
        # ( (), "PU140-JPsi", "/JPsiToMuMu_Pt0to100-pythia8_TuneCP5-gun/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_111X_mcRun4_realistic_T15_v1-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
        # ( ("Eff"), "PU200-JPsi", "/JPsiToMuMu_Pt0to100-pythia8_TuneCP5-gun/Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext1-v3/GEN-SIM-DIGI-RAW-MINIAOD"),
    ]

    PWD = os.getcwd()

    # HERE
    nfilesDY    = 30
    nfilesTT    = 15
    nfilesOther = 30

    cmds = []
    for i, (analyzers, tag, sample) in enumerate(samples):
        pd   = sample.split('/')[1]
        name = ""

        if type(analyzers) != tuple:
            analyzers = (analyzers, )

        if 'SingleMuon' in sample:
            name = sample.split('/')[1]+'_'+sample.split('/')[2]
        else:
            name = tag

        crabName = "crab_"+name+crabVER

        crabDir = BASE + pd + "/" + crabName

        subdirs = [ x[0] for x in os.walk( crabDir ) ]

        paths = []
        for subdir in subdirs:
            if "/000" in subdir:
                paths.append( subdir )

        nfiles = nfilesOther
        if "DYToLL_M" in name:
            nfiles = nfilesDY
        elif "TT" in name:
            nfiles = nfilesTT

        for an in analyzers:

            if doHadd:
                if os.path.isfile('Outputs_%(histoVER)s/%(an)s/hist-%(histoVER)s-%(name)s-%(an)s.root' % locals()):
                    continue

                cmd = "hadd Outputs_%(histoVER)s/%(an)s/hist-%(histoVER)s-%(name)s-%(an)s.root Outputs_%(histoVER)s/%(an)s/%(name)s/hist-%(histoVER)s-%(name)s--Job*-%(an)s.root >Outputs_%(histoVER)s/%(an)s/hadd-%(name)s-%(an)s.log" % locals()

                print ""
                print cmd
                sys.stdout.flush()
                gc.collect()

                os.system(cmd)
                sys.stdout.flush()
                gc.collect()

                # if i % 5 == 4:
                #     time.sleep(60)

            else:
                for ipath, path in enumerate(paths):

                    jobid_files = jobSpritting(path, nfiles)

                    for jobid, files in jobid_files:

                        strJobId = 'Job'+str(ipath)+"0"+str(jobid)

                        doDimuon = "false"
                        if "DYToLL_M" in name or "Zprime_M6000" in name:
                            doDimuon = "true"

                        doGenMatchForIso = "false"
                        if "QCD_" not in name:
                            doGenMatchForIso = "true"

                        cmd = "source ./Phase2HLT-%s.sh %s %s %s %s %s %s %s" % (an, PWD, histoVER, name, files, strJobId, doDimuon, doGenMatchForIso)

                        if doRecover:
                            outpath = "Outputs_%s/%s/%s/hist-%s-%s--%s-%s.root" % (histoVER, an, name, histoVER, name, strJobId, an)
                            if not os.path.isfile(outpath):
                                print ""
                                print cmd
                                sys.stdout.flush()
                                gc.collect()

                                os.system(cmd)
                                sys.stdout.flush()
                                gc.collect()

                        else:
                            print ""
                            print cmd
                            sys.stdout.flush()
                            gc.collect()

                            os.system(cmd)
                            sys.stdout.flush()
                            gc.collect()

                        # HERE
                        # if jobid > 3:
                        #     sys.exit()

    print ""
    print "finished"
    sys.stdout.flush()
    gc.collect()





