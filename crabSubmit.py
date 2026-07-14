from CRABClient.UserUtilities import config, getUsername
import sys, os
import gc
import datetime
now = datetime.datetime.now()
date = now.strftime('%Y%m%d')

submitVersion = 'MuonHLTPhase2_cmssw1600pre3'
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
config.Data.unitsPerJob = 1
config.JobType.maxMemoryMB = 8000
config.JobType.numCores = 4
config.JobType.maxJobRuntimeMin = 2749

config.Data.outLFNDirBase = '%(mainOutputDir)s'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
# config.Site.whitelist = ['T2_CH_CERN', '']
    """

    datasets = [
        ("DYToLL_M50_Phase2Spring24", "/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD"),
    ]

    HLT_menus = [
        "hlt_muon_mc.py",   
        "hlt_muon_mc_newTRKbaseline.py",   
        # "hlt_muon_mc_OIFirst_OIFromL1_dump.py",
        # "hlt_muon_mc_rerun_1500pre3_default.py",
        # "hlt_muon_mc_NoDupl_ALPAKA_Iter0.py",
        # "hlt_muon_mc_NoDupl_Iter02.py",
        # "hlt_muon_mc_NoDupl_Iter0.py"
        
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


