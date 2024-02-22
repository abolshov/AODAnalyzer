import FWCore.ParameterSet.Config as cms

process = cms.Process("AODAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring('file:B2G-Run3Summer22DRPremix-00441.root')
)

process.demo = cms.EDAnalyzer('AODAnalyzer',
    genparticles = cms.untracked.InputTag('genParticles'),
    genjets = cms.untracked.InputTag('ak4GenJets')
)

process.p = cms.Path(process.demo)