import FWCore.ParameterSet.Config as cms

process = cms.Process("testParticle")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
ev1 = 83
ev2 = ev1 + 1
process.source = cms.Source("PoolSource",
    eventsToProcess = cms.untracked.VEventRange("1:%d-1:%d"%(ev1, ev2)),
    #eventsToProcess = cms.untracked.VEventRange("1:%d-1:%d"%(4119, 4120)),
    #fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/a/abolshov/private/CMSSW_12_4_16/src/B2G-Run3Summer22DRPremix-00441.root')
    fileNames = cms.untracked.vstring('file:B2G-Run3Summer22DRPremix-00441.root')
    #fileNames = cms.untracked.vstring('file:/eos/user/a/abolshov/Analysis/HME_files/SingleLepton/MiniAOD/287CBF44-3E3D-E911-B46D-20040FE94274.root')
    #fileNames = cms.untracked.vstring('file:/eos/user/t/tahuang/GluGluToRadionToHHTo2B2WToLNu2J_M-250_Run2016/6C2B2F7D-DF3C-E911-B863-B083FED00117.root')
    #fileNames = cms.untracked.vstring('file:/eos/user/a/abolshov/Analysis/HME_files/SingleLepton/NanoAOD/NanoAOD_600GeV_1000Events.root')
   #fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_b2g/wprime/muahmad/ZgammaJetsToBB_Mass_13TeV_MadGraph5_pythia8_RECO_2018/ZgammaJetsToBB_13TeV_MadGraph5_pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1-RECO-ZgammaJetsToBB-v2/231117_220404/0000/SUS-RunIISummer20UL18RECO-00039_772.root')
)

process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    #src = cms.InputTag("prunedGenParticles"),
    maxEventsToPrint = cms.untracked.int32(1),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex = cms.untracked.bool(True),
    printStatus = cms.untracked.bool(True),
    printIndex  = cms.untracked.bool(True)
)

process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genParticles"),
    #src = cms.InputTag("prunedGenParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex  = cms.untracked.bool(True)
)

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

#process.p = cms.Path(process.printTree1*process.printTree2)
process.p = cms.Path(process.printTree2)
process.outpath = cms.EndPath(process.printEventNumber)


