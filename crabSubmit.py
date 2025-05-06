from CRABClient.UserUtilities import config, getUsername
import sys, os
import gc
import datetime
now = datetime.datetime.now()
date = now.strftime('%Y%m%d')

submitVersion = 'MuonHLTRun3_cmssw1409'
mainOutputDir = '/store/group/phys_muon/ec/HLT/%s/%s' % (submitVersion, date)
#mainOutputDir = '/store/user/%s/%s/%s' % (getUsername(), submitVersion, date)


# 'MultiCRAB' part
if __name__ == '__main__':
    # from CRABAPI.RawCommand import crabCommand

    crab_cfg = """
from CRABClient.UserUtilities import config, getUsername

config = config()

config.General.requestName = '%(datasetTag)s_%(menuTag)s_%(date)s'
config.General.workArea = 'crab_%(submitVersion)s_%(date)s'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '%(menu)s'

config.Data.inputDataset = '%(datasetPath)s'
config.Data.useParent = False                   ## Only MC
config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20    # 1(DY) # 20(Data)
#config.Data.totalUnits = 800 # 1000(DY)                                                                                                      ## Only MC
#config.Data.lumiMask = 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_Last5invfb.txt'                        ## Only Data
#config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json'   ## Only Data
#config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json'   ## Only Data
#config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_381594_Golden.json'   ## Only Data
config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/DCSOnly_JSONS/dailyDCSOnlyJSON/Collisions24_13p6TeV_378981_382344_DCSOnly_TkPx.json'
#config.Data.lumiMask = 'Run380237.json'

config.JobType.maxJobRuntimeMin = 2749
config.Data.outLFNDirBase = '%(mainOutputDir)s'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_KR_KNU'
    """

    datasets = [
        #("DYToLL_M50_133X", "/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/Run3Winter24Reco-KeepSi_133X_mcRun3_2024_realistic_v8-v2/AODSIM"),
        #("Bs_133X", "/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV-pythia8-evtgen/Run3Winter24Reco-133X_mcRun3_2024_realistic_v8-v2/AODSIM"),
        #("JPsi_133X", "/JPsiToMuMu_PT-0to100_pythia8-gun/Run3Winter24Reco-KeepSi_133X_mcRun3_2024_realistic_v8-v2/AODSIM"),

        #("DYToLL_M50_130X", "/DYto2L_M-50_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixDRPremix-KeepSi_130X_mcRun3_2023_realistic_postBPix_v2-v3/AODSIM"),
        #("Bs_126X", "/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/Run3Winter23Reco-KeepSi_126X_mcRun3_2023_forPU65_v1-v2/AODSIM"),
        #("JPsi_126X", "/JPsiTo2Mu_Pt-0To100_pythia8-gun/Run3Winter23Reco-KeepSi_126X_mcRun3_2023_forPU65_v1-v2/AODSIM"),
        #("DYToLL_M50_126X", "/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/AODSIM"),
        #("Zprime_126X", "/ZprimeToMuMu_M-6000_TuneCP5_13p6TeV_pythia8/Run3Winter23Reco-KeepSi_126X_mcRun3_2023_forPU65_v1-v2/AODSIM"),

        #("Muon0_Run2024B", "/Muon0/Run2024B-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon1_Run2024B", "/Muon1/Run2024B-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon0_Run2024C", "/Muon0/Run2024C-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon1_Run2024C", "/Muon1/Run2024C-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon0_Run2024D", "/Muon0/Run2024D-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon1_Run2024D", "/Muon1/Run2024D-ZMu-PromptReco-v1/RAW-RECO"),
        ("Muon0_Run2024Ev1", "/Muon0/Run2024E-ZMu-PromptReco-v1/RAW-RECO"),
        ("Muon1_Run2024Ev1", "/Muon1/Run2024E-ZMu-PromptReco-v1/RAW-RECO"),
        ("Muon0_Run2024Ev2", "/Muon0/Run2024E-ZMu-PromptReco-v2/RAW-RECO"),
        ("Muon1_Run2024Ev2", "/Muon1/Run2024E-ZMu-PromptReco-v2/RAW-RECO"),
        ("Muon0_Run2024F", "/Muon0/Run2024F-ZMu-PromptReco-v1/RAW-RECO"),
        ("Muon1_Run2024F", "/Muon1/Run2024F-ZMu-PromptReco-v1/RAW-RECO"),

        #("Muon0_Run2023B", "/Muon0/Run2023B-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon1_Run2023B", "/Muon1/Run2023B-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon0_Run2023C", "/Muon0/Run2023C-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon1_Run2023C", "/Muon1/Run2023C-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon0_Run2023Cv2", "/Muon0/Run2023C-ZMu-PromptReco-v2/RAW-RECO"),
        #("Muon1_Run2023Cv2", "/Muon1/Run2023C-ZMu-PromptReco-v2/RAW-RECO"),
        #("Muon0_Run2023Cv3", "/Muon0/Run2023C-ZMu-PromptReco-v3/RAW-RECO"),
        #("Muon1_Run2023Cv3", "/Muon1/Run2023C-ZMu-PromptReco-v3/RAW-RECO"),
        #("Muon0_Run2023Cv4", "/Muon0/Run2023C-ZMu-PromptReco-v4/RAW-RECO"),
        #("Muon1_Run2023Cv4", "/Muon1/Run2023C-ZMu-PromptReco-v4/RAW-RECO"),
        #("Muon0_Run2023Dv1", "/Muon0/Run2023D-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon1_Run2023Dv1", "/Muon1/Run2023D-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon0_Run2023Dv2", "/Muon0/Run2023D-ZMu-PromptReco-v2/RAW-RECO"),
        #("Muon1_Run2023Dv2", "/Muon1/Run2023D-ZMu-PromptReco-v2/RAW-RECO"),

        #("SingleMuon_Run2022B", "/SingleMuon/Run2022B-ZMu-10Dec2022-v2/RAW-RECO"),
        #("SingleMuon_Run2022C", "/SingleMuon/Run2022C-ZMu-10Dec2022-v2/RAW-RECO"),
        #("Muon_Run2022C", "/Muon/Run2022C-ZMu-27Jun2023-v1/RAW-RECO"),
        #("Muon_Run2022D", "/Muon/Run2022D-ZMu-27Jun2023-v2/RAW-RECO"),
        #("Muon_Run2022E", "/Muon/Run2022E-ZMu-27Jun2023-v1/RAW-RECO"),
        #("Muon_Run2022F", "/Muon/Run2022F-ZMu-PromptReco-v1/RAW-RECO"),
        #("Muon_Run2022G", "/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO"),

        #("SingleMuon_RunUL2018D", "/SingleMuon/Run2018D-ZMu-12Nov2019_UL2018-v4/RAW-RECO"),
    ]

    HLT_menus = [
        #"hlt_muon_mc_forTraining.py",
        "hlt_muon_data.py",
        "hlt_muon_data_CSC.py",
        "hlt_muon_data_CSC_RPC.py",
        #"hlt_muon_data_chaining.py",
        #"hlt_muon_data_chaining_v2.py",
        #"hlt_muon_data_BDT_PR.py",
        #"hlt_muon_data_CCC.py",
        #"hlt_muon_mc.py",
    ]

    # proxy = '"/tmp/x509up_u95096"'

    for menu in HLT_menus:
        #menuTag = 'HLTRun3'  # menu.replace(".py", "").replace("HLT_MC_", "")
        menuTag = menu.replace(".py", "").replace("HLT_", "")

        for datasetTag, datasetPath in datasets:

            Crab_Config = 'crabConfig_'+datasetTag+'.py'
            print( "\n\n", crab_cfg % locals() )
            sys.stdout.flush()
            gc.collect()

            open(Crab_Config, 'wt').write(crab_cfg % locals())

            cmd = 'crab submit -c '+Crab_Config  # +' --proxy='+proxy

            print( cmd )
            sys.stdout.flush()
            gc.collect()

            os.system(cmd)


