import FWCore.ParameterSet.Config as cms

process = cms.Process( "TEST" )
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),allowUnscheduled=cms.untracked.bool(True))
#,
#				     SkipEvent = cms.untracked.vstring('ProductNotFound'))
filterMode = False # True                
 
######## Sequence settings ##########
corrJetsOnTheFly = False
runOnMC  = True
runOnSig = True
DOHLTFILTERS = True
#useJSON = not (runOnMC)
#JSONfile = 'Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#****************************************************************************************************#

#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
   #process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'#'MCRUN2_74_V9::All'
   process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v15'#'MCRUN2_74_V9::All'
elif not(runOnMC):
   process.GlobalTag.globaltag = '102X_dataRun2_v12'

hltFiltersProcessName = 'RECO'
if runOnMC:
   hltFiltersProcessName = 'PAT' #'RECO'

#if DOHLTFILTERS and not(runOnMC):
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
                                                    inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
                                                    reverseDecision = cms.bool(False)
                                                    )
process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
                                                       inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
                                                       reverseDecision = cms.bool(False)
                                                       )
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )

######### read JSON file for data ##########					                                                             
'''if not(runOnMC) and useJSON:

  import FWCore.PythonUtilities.LumiList as LumiList
  import FWCore.ParameterSet.Types as CfgTypes
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
  process.source.lumisToProcess.extend(myLumis) 
'''

# ---------------------------------------------------------
# DeepAK8: set up TransientTrackBuilder
process.load('Configuration.StandardSequences.MagneticField_cff')
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)
# ---------------------------------------------------------

####### Redo Jet clustering sequence ##########

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS, ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedMass, ak8PFJetsCHSSoftDropMass# , ak8PFJetsCSTrimmed, ak8PFJetsCSFiltered, ak8PFJetsCHSFilteredMass, ak8PFJetsCHSTrimmedMass

from CommonTools.PileupAlgos.Puppi_cff import puppi

process.puppi = puppi.clone()
process.puppi.useExistingWeights = True
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

process.ak8PFJetsCHS = ak8PFJetsCHS.clone( src = 'puppi', jetPtMin = 100.0 )
process.ak8PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = 'puppi', jetPtMin = 100.0 )
process.ak8PFJetsCHSPrunedMass = ak8PFJetsCHSPrunedMass.clone()
process.ak8PFJetsCHSSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'puppi', jetPtMin = 100.0 )
process.ak8PFJetsCHSSoftDropMass = ak8PFJetsCHSSoftDropMass.clone()


process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
                   src = cms.InputTag("ak8PFJetsCHS"),
                   Njets = cms.vuint32(1, 2, 3, 4),
                   # variables for measure definition : 
                   measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                   beta = cms.double(1.0),          # CMS default is 1
                   R0 = cms.double( 0.8 ),          # CMS default is jet cone size
                   Rcutoff = cms.double( 999.0),       # not used by default
                   # variables for axes definition :
                   axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                   nPass = cms.int32(0),         # not used by default
                   akAxesR0 = cms.double(-999.0)        # not used by default
                   )

process.substructureSequence = cms.Sequence()
process.substructureSequence+=process.puppi
process.substructureSequence+=process.ak8PFJetsCHS
process.substructureSequence+=process.NjettinessAK8
process.substructureSequence+=process.ak8PFJetsCHSPruned
process.substructureSequence+=process.ak8PFJetsCHSPrunedMass
process.substructureSequence+=process.ak8PFJetsCHSSoftDrop
process.substructureSequence+=process.ak8PFJetsCHSSoftDropMass

####### Redo pat jets sequence ##########
process.redoPatJets = cms.Sequence()
process.redoPrunedPatJets = cms.Sequence()
process.redoSoftDropPatJets = cms.Sequence()

from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

# Redo pat jets from ak8PFJetsCHS
process.patJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHS' )
process.patJetsAK8 = patJetsAK8.clone( jetSource = 'ak8PFJetsCHS' )
process.patJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8PFJetsCHSPrunedMass"), cms.InputTag("ak8PFJetsCHSSoftDropMass"), cms.InputTag("NjettinessAK8:tau1"), cms.InputTag("NjettinessAK8:tau2"), cms.InputTag("NjettinessAK8:tau3"),cms.InputTag("NjettinessAK8:tau4")]
process.patJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8") )
process.selectedPatJetsAK8 = selectedPatJetsAK8.clone( cut = cms.string('pt > 100') )

process.redoPatJets+=process.patJetCorrFactorsAK8
process.redoPatJets+=process.patJetsAK8
process.redoPatJets+=process.selectedPatJetsAK8

# Redo pat jets ak8PFJetsCHSPruned
process.patJetCorrFactorsAK8Pruned = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSPruned' )
process.patJetsAK8Pruned = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSPruned' )
process.patJetsAK8Pruned.userData.userFloats.src = [ "" ]
#process.patJetsAK8Pruned.userData.userFloats =cms.PSet(src = cms.VInputTag(""))
process.patJetsAK8Pruned.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Pruned") )
process.selectedPatJetsAK8Pruned = selectedPatJetsAK8.clone(cut = 'pt > 100', src = "patJetsAK8Pruned")

process.redoPrunedPatJets+=process.patJetCorrFactorsAK8Pruned
process.redoPrunedPatJets+=process.patJetsAK8Pruned
process.redoPrunedPatJets+=process.selectedPatJetsAK8Pruned

# Redo pat jets ak8PFJetsCHSSoftDrop
process.patJetCorrFactorsAK8Softdrop = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHSSoftDrop' )
process.patJetsAK8Softdrop = patJetsAK8.clone( jetSource = 'ak8PFJetsCHSSoftDrop' )
process.patJetsAK8Softdrop.userData.userFloats.src = [ "" ]
#process.patJetsAK8Softdrop.userData.userFloats =cms.PSet(src = cms.VInputTag(""))
process.patJetsAK8Softdrop.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("patJetCorrFactorsAK8Softdrop") )
process.selectedPatJetsAK8Softdrop = selectedPatJetsAK8.clone(cut = 'pt > 100', src = "patJetsAK8Softdrop")

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
## PATify soft drop subjets
addJetCollection(
    process,
    labelName = 'AK8SoftDropSubjets',
    jetSource = cms.InputTag('ak8PFJetsCHSSoftDrop','SubJets'),
    algo = 'ak',  # needed for subjet flavor clustering
    rParam = 0.8, # needed for subjet flavor clustering
    getJetMCFlavour = False,
    pvSource = cms.InputTag( 'offlineSlimmedPrimaryVertices' ),
    genJetCollection = cms.InputTag('slimmedGenJets'),
    genParticles = cms.InputTag( 'prunedGenParticles' ),
    btagDiscriminators = ['None'],
    jetCorrections = ('AK4PFPuppi', ['L2Relative', 'L3Absolute'], 'None'),
#    explicitJTA = True,  # needed for subjet b tagging
#    svClustering = True, # needed for subjet b tagging
#    fatJets=cms.InputTag('ak8PFJetsCHS'),             # needed for subjet flavor clustering
#    groomedFatJets=cms.InputTag('ak8PFJetsCHSSoftDrop') # needed for subjet flavor clustering
)
process.selectedPatJetsAK8SoftDropPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc = cms.InputTag("selectedPatJetsAK8Softdrop"),
    subjetSrc = cms.InputTag("selectedPatJetsAK8SoftDropSubjets")
)


process.redoSoftDropPatJets+=process.patJetCorrFactorsAK8Softdrop
process.redoSoftDropPatJets+=process.patJetsAK8Softdrop
process.redoSoftDropPatJets+=process.selectedPatJetsAK8Softdrop


option = 'RECO'

process.load("ExoDiBosonResonances.EDBRCommon.goodMuons_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodElectrons_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodJets_cff")
process.load("ExoDiBosonResonances.EDBRCommon.leptonicW_cff")
process.load("ExoDiBosonResonances.EDBRCommon.hadronicW_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodPuppi_cff")

if option == 'RECO':
    process.goodMuons.src = "slimmedMuons"
    process.goodElectrons.src = "slimmedElectrons"
    process.goodJets.src = "slimmedJetsAK8"
#    process.goodJets.src = "selectedPatJetsAK8"
    process.Wtoenu.MET  = "slimmedMETs"
    process.Wtomunu.MET = "slimmedMETs"
    #process.goodPuppi.src = "selectedPatJetsAK8"
    process.goodPuppi.src = "JetUserData"
    process.goodPuppiAK4.src = "JetUserDataak4"

process.goodOfflinePrimaryVertex = cms.EDFilter("VertexSelector",
                                       src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                       cut = cms.string("chi2!=0 && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0"),
                                       filter = cms.bool(True)
                                       )
if option == 'RECO':
    process.hadronicV.cut = ' '
if option == 'GEN':
    process.hadronicV.cut = ' '
WBOSONCUT = "pt > 200.0"

process.leptonicVSelector = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("leptonicV"),
                                       cut = cms.string( WBOSONCUT ),
                                       filter = cms.bool(True)
                                       )
process.leptonicVFilter = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("leptonicV"),
                                       minNumber = cms.uint32(1),
                                       #filter = cms.bool(True)
                                       )
process.hadronicVFilter = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("hadronicV"),
                                       minNumber = cms.uint32(1),
                                       #filter = cms.bool(True)
                                       )
process.graviton = cms.EDProducer("CandViewCombiner",
                                       decay = cms.string("leptonicV hadronicV"),
                                       checkCharge = cms.bool(False),
                                       cut = cms.string("mass > 180"),
                                       roles = cms.vstring('leptonicV', 'hadronicV'),
                                       )
process.gravitonFilter =  cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("graviton"),
                                       minNumber = cms.uint32(1),
                                       #filter = cms.bool(True)
                                       )
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.leptonSequence = cms.Sequence(process.muSequence +
                                      process.egmGsfElectronIDSequence*process.eleSequence +
                                      process.leptonicVSequence +
                                      process.leptonicVSelector +
                                      process.leptonicVFilter )

process.jetSequence = cms.Sequence(process.substructureSequence +
                                   process.redoPatJets + 
                                   process.redoPrunedPatJets+
                                   process.redoSoftDropPatJets+
                                   process.fatJetsSequence +
                                   process.fatPuppiSequence+
                                   process.hadronicV +
                                   process.hadronicVFilter)

process.gravitonSequence = cms.Sequence(process.graviton +
                                        process.gravitonFilter)


if filterMode == False:
    process.goodOfflinePrimaryVertex.filter = False
    process.Wtomunu.cut = ''
    process.Wtoenu.cut = ''
    process.leptonicVSelector.filter = False
    process.leptonicVSelector.cut = ''
    process.hadronicV.cut = ''
    process.graviton.cut = ''
    process.leptonicVFilter.minNumber = 0
    process.hadronicVFilter.minNumber = 0
    process.gravitonFilter.minNumber = 0

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.metfilterSequence = cms.Sequence(process.BadPFMuonFilter+process.BadChargedCandidateFilter)

######### JEC ########
METS = "slimmedMETs"
jetsAK8 = "slimmedJetsAK8"
jetsAK8pruned = "slimmedJetsAK8"
jetsAK8softdrop = "slimmedJetsAK8"
jetsAK8puppi = "cleanPuppi"
 
if runOnMC:
   jecLevelsAK8chs = [
                                   'Fall17_17Nov2017_V23_MC_L1FastJet_AK8PFchs.txt',
                                   'Fall17_17Nov2017_V23_MC_L2Relative_AK8PFchs.txt',
                                   'Fall17_17Nov2017_V23_MC_L3Absolute_AK8PFchs.txt'
     ]
   jecLevelsAK8chsGroomed = [
                                   'Fall17_17Nov2017_V23_MC_L2Relative_AK8PFchs.txt',
                                   'Fall17_17Nov2017_V23_MC_L3Absolute_AK8PFchs.txt'
     ]
   jecLevelsAK8puppi = [
                                   'Fall17_17Nov2017_V23_MC_L1FastJet_AK8PFPuppi.txt',
                                   'Fall17_17Nov2017_V23_MC_L2Relative_AK8PFPuppi.txt',
                                   'Fall17_17Nov2017_V23_MC_L3Absolute_AK8PFPuppi.txt'
     ]
   jecLevelsAK8puppiGroomed = [
                                   'Fall17_17Nov2017_V23_MC_L2Relative_AK8PFPuppi.txt',
                                   'Fall17_17Nov2017_V23_MC_L3Absolute_AK8PFPuppi.txt'
     ]
   BjecLevelsAK4chs = [
                                   'Fall17_17Nov2017_V23_MC_L1FastJet_AK4PFPuppi.txt',
                                   'Fall17_17Nov2017_V23_MC_L2Relative_AK4PFPuppi.txt',
                                   'Fall17_17Nov2017_V23_MC_L3Absolute_AK4PFPuppi.txt'
     ]
   jecLevelsAK4chs = [
                                   'Fall17_17Nov2017_V23_MC_L1FastJet_AK4PFchs.txt',
                                   'Fall17_17Nov2017_V23_MC_L2Relative_AK4PFchs.txt',
                                   'Fall17_17Nov2017_V23_MC_L3Absolute_AK4PFchs.txt'
    ]
else:
   jecLevelsAK8chs = [
                                   'Autumn18_RunA_V8_DATA_L1FastJet_AK8PFchs.txt',
                                   'Autumn18_RunA_V8_DATA_L2Relative_AK8PFchs.txt',
                                   'Autumn18_RunA_V8_DATA_L3Absolute_AK8PFchs.txt',
                                   'Autumn18_RunA_V8_DATA_L2L3Residual_AK8PFchs.txt'
     ]
   jecLevelsAK8chsGroomed = [
                                   'Autumn18_RunA_V8_DATA_L2Relative_AK8PFchs.txt',
                                   'Autumn18_RunA_V8_DATA_L3Absolute_AK8PFchs.txt',
                                   'Autumn18_RunA_V8_DATA_L2L3Residual_AK8PFchs.txt'
     ]
   jecLevelsAK8puppi = [
                                   'Autumn18_RunA_V8_DATA_L1FastJet_AK8PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L2Relative_AK8PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L3Absolute_AK8PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L2L3Residual_AK8PFPuppi.txt'
     ]
   jecLevelsAK8puppiGroomed = [
                                   'Autumn18_RunA_V8_DATA_L2Relative_AK8PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L3Absolute_AK8PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L2L3Residual_AK8PFPuppi.txt'
     ]
   BjecLevelsAK4chs = [
                                   'Autumn18_RunA_V8_DATA_L1FastJet_AK4PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L2Relative_AK4PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L3Absolute_AK4PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L2L3Residual_AK4PFPuppi.txt'

     ]
   jecLevelsAK4chs = [
                                   'Autumn18_RunA_V8_DATA_L1FastJet_AK4PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L2Relative_AK4PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L3Absolute_AK4PFPuppi.txt',
                                   'Autumn18_RunA_V8_DATA_L2L3Residual_AK4PFPuppi.txt'
     ]
#jLabel = jetsAK8puppi
jLabel = "slimmedJetsAK8"
jetAlgo    = 'AK8PFPuppi'
jer_era = "Fall17_17Nov2017_V23_MC"
triggerResultsLabel      = "TriggerResults"
triggerSummaryLabel      = "hltTriggerSummaryAOD"
hltProcess = "HLT"
process.JetUserData = cms.EDProducer(
                                     'JetUserData',
                                     jetLabel          = cms.InputTag(jLabel),
                                     rho               = cms.InputTag("fixedGridRhoFastjetAll"),
                                     coneSize          = cms.double(0.8),
                                     getJERFromTxt     = cms.bool(False),
                                     jetCorrLabel      = cms.string(jetAlgo),
                                     jerLabel          = cms.string(jetAlgo),
                                     resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgo+'.txt'),
                                     scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgo+'.txt'),
                                     ### TTRIGGER ###
                                     triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
                                     triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
                                     hltJetFilter       = cms.InputTag("hltPFHT"),
                                     hltPath            = cms.string("HLT_PFHT800"),
                                     hlt2reco_deltaRmax = cms.double(0.2),
                                     candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"),
                                     jecAk8chsPayloadNames_jetUserdata = cms.vstring( jecLevelsAK8puppi ),
                                     vertex_jetUserdata = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                     )
jLabelak4 = "slimmedJetsPuppi"
jetAlgoak4    = 'AK4PFPuppi'
jer_era = "Fall17_17Nov2017_V23_MC"
triggerResultsLabel      = "TriggerResults"
triggerSummaryLabel      = "hltTriggerSummaryAOD"
hltProcess = "HLT"
process.JetUserDataak4 = cms.EDProducer(
                                        'JetUserDataak4',
                                        jetLabel          = cms.InputTag(jLabelak4),
                                        rho               = cms.InputTag("fixedGridRhoFastjetAll"),
                                        coneSize          = cms.double(0.4),
                                        getJERFromTxt     = cms.bool(False),
                                        jetCorrLabel      = cms.string(jetAlgoak4),
                                        jerLabel          = cms.string(jetAlgoak4),
                                        resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgoak4+'.txt'),
                                        scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgoak4+'.txt'),
                                        ### TTRIGGER ###
                                        triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
                                        triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
                                        hltJetFilter       = cms.InputTag("hltPFHT"),
                                        hltPath            = cms.string("HLT_PFHT800"),
                                        hlt2reco_deltaRmax = cms.double(0.2),
                                        candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"),
                                        jecAK4chsPayloadNames_JetUserData = cms.vstring( BjecLevelsAK4chs ),
                                        vertex_JetUserData = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        )
# L1 prefiring
process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                 ThePhotons = cms.InputTag("slimmedPhotons"),
                                 TheJets = cms.InputTag("slimmedJets"),
#                                L1Maps = cms.string(relBase+"/src/L1Prefiring/EventWeightProducer/files/L1PrefiringMaps_new.root"),
                                # L1Maps = cms.string("L1PrefiringMaps_new.root"), # update this line with the location of this file
                                L1Maps = cms.string("L1PrefiringMaps_new.root"),
                                 DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
                                 UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                 PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                 )
process.treeDumper = cms.EDAnalyzer("EDBRTreeMaker",
                                    originalNEvents = cms.int32(1),
                                    crossSectionPb = cms.double(1),
                                    targetLumiInvPb = cms.double(1.0),
                                    EDBRChannel = cms.string("VW_CHANNEL"),
                                    lhe =  cms.InputTag("externalLHEProducer"),
                                    isGen = cms.bool(False),
                                    isJEC = cms.bool(corrJetsOnTheFly),
                                    RunOnMC  = cms.bool(runOnMC),
                                    RunOnSig = cms.bool(runOnSig),
                                    generator =  cms.InputTag("generator"),
                                    genSrc =  cms.InputTag("prunedGenParticles"),
                                    pileup  =   cms.InputTag("slimmedAddPileupInfo"),
                                    leptonicVSrc = cms.InputTag("leptonicV"),
                                    gravitonSrc = cms.InputTag("graviton"),
                                    t1jetSrc_userak4 = cms.InputTag("JetUserDataak4"),
                                    looseMuonSrc = cms.InputTag("looseMuons"),
                                    looseElectronSrc = cms.InputTag("looseElectrons"),
                                    goodMuSrc = cms.InputTag("goodMuons"),
                                    MuSrc = cms.InputTag("slimmedMuons"),
                                    EleSrc = cms.InputTag("slimmedElectrons"),
                                    t1muSrc = cms.InputTag("slimmedMuons"),
                                    metSrc = cms.InputTag("slimmedMETs"),
                                    mets = cms.InputTag(METS),
                                    #ak4jetsSrc = cms.InputTag("cleanAK4Jets"), 
                                    ak4jetsSrc = cms.InputTag("cleanPuppiAK4"), 
                                    hadronicVSrc = cms.InputTag("hadronicV"),
                                    hadronicVSrc_raw = cms.InputTag("slimmedJetsAK8"),
                                    hadronicVSoftDropSrc = cms.InputTag("selectedPatJetsAK8SoftDropPacked"),
                                    jets = cms.InputTag("slimmedJets"),
                                    fatjets = cms.InputTag(jetsAK8),
                                    ak8JetSrc = cms.InputTag(jetsAK8),
                                    prunedjets = cms.InputTag(jetsAK8pruned),
                                    softdropjets = cms.InputTag(jetsAK8softdrop),
                                    puppijets = cms.InputTag(jetsAK8puppi),
                                    jecAK8chsPayloadNames = cms.vstring( jecLevelsAK8chs ),
                                    jecAK8chsPayloadNamesGroomed = cms.vstring( jecLevelsAK8chsGroomed ),
                                    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
                                    BjecAK4chsPayloadNames = cms.vstring( BjecLevelsAK4chs ),
					
                                    jecAK8puppiPayloadNames = cms.vstring( jecLevelsAK8puppi ),
                                    jecAK8puppiPayloadNamesGroomed = cms.vstring( jecLevelsAK8puppiGroomed ),
                                    jecpath = cms.string(''),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    electronIDs = cms.InputTag("heepElectronID-HEEPV50-CSA14-25ns"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    hltToken    = cms.InputTag("TriggerResults","","HLT"),
                                    elPaths1     = cms.vstring("HLT_Ele32_WPTight_Gsf_v*"),#EXO-15-002#1
                                    elPaths2     = cms.vstring("HLT_Ele35_WPTight_Gsf_v*"), #B2G-15-005#2
                                    elPaths3     = cms.vstring("HLT_Ele45_WPLoose_Gsf_v*"),
                                    elPaths4     = cms.vstring("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"),#("HLT_Ele35_WPLoose_Gsf_v*"),#3
                                    elPaths5     = cms.vstring("HLT_Photon200_v*"),#("HLT_Ele35_WPLoose_Gsf_v*"),
                                    #elPaths5     = cms.vstring("HLT_Ele25_WPTight_Gsf_v*"),
                                    elPaths6     = cms.vstring("HLT_Ele38_WPTight_Gsf_v*"),#("HLT_Ele25_eta2p1_WPLoose_Gsf_v*"),
                                    elPaths7     = cms.vstring("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*"),#("HLT_Ele25_eta2p1_WPTight_Gsf_v*"),#4
                                    elPaths8     = cms.vstring("HLT_Ele27_WPTight_Gsf_v*"),
                                    muPaths1     = cms.vstring("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v*"),#EXO-15-002
                                    muPaths2     = cms.vstring("HLT_Mu50_v*"), #B2G-15-005
                                    muPaths3     = cms.vstring("HLT_TkMu50_v*"), #B2G-15-005
                                    muPaths4     = cms.vstring("HLT_OldMu100_v*"), #MET
                                    muPaths5     = cms.vstring("HLT_TkMu100_v*"), #MET
                                    muPaths6     = cms.vstring("HLT_PFHT900_v*"),
                                    muPaths7     = cms.vstring("HLT_PFJet450_v*"),
                                    muPaths8     = cms.vstring("HLT_PFJet500_v*"),
                                    muPaths9     = cms.vstring("HLT_AK8PFJet450_v*"),
                                    muPaths10     = cms.vstring("HLT_AK8PFJet500_v*"),
                                    muPaths11     = cms.vstring("HLT_AK8PFJet360_TrimMass30_v*"),#e1,MET110,e7,m1,m5-m11,HLT_hadronic
                                    muPaths12     = cms.vstring("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"),
                                    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
                                    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
                                    noiseFilterSelection_HBHENoiseIsoFilter = cms.string("Flag_HBHENoiseIsoFilter"),
                                    noiseFilterSelection_GlobalTightHaloFilter = cms.string('Flag_globalSuperTightHalo2016Filter'),
                                    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
                                    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
                                    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
                                    noiseFilterSelection_badMuon = cms.InputTag('BadPFMuonFilter'),
                                    noiseFilterSelection_badChargedHadron = cms.InputTag('BadChargedCandidateFilter'),
                                    )


if option=='GEN':
    process.treeDumper.metSrc = 'genMetTrue'
    process.treeDumper.isGen  = True
 

process.analysis = cms.Path(process.JetUserData +
                            process.JetUserDataak4 +
                            process.leptonSequence +
                            #process.substructureSequence+
                            #process.redoPatJets+
                            #process.redoPrunedPatJets+
                            #process.redoSoftDropPatJets+
                            process.HBHENoiseFilterResultProducer+
                            process.ApplyBaselineHBHENoiseFilter+
                            process.ApplyBaselineHBHEIsoNoiseFilter+
                            process.jetSequence +
                            process.metfilterSequence +
                            process.gravitonSequence +
                            process.ecalBadCalibReducedMINIAODFilter*process.prefiringweight*process.treeDumper)

if option=='RECO':
    process.analysis.replace(process.leptonSequence, process.goodOfflinePrimaryVertex + process.leptonSequence)

process.load("ExoDiBosonResonances.EDBRCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
#process.source.inputCommands = ['keep *','drop *_isolatedTracks_*_*']
process.source.fileNames = [
#'/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/70000/FC761907-BF54-E811-9076-0242AC130002.root',
#'/store/mc/RunIIFall17MiniAODv2/WZTo3LNu_3Jets_MLL-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/70000/F2EF406C-3F65-E811-AA59-0025905C43EC.root'
#'/store/mc/RunIIAutumn18MiniAOD/WW_TuneCP5_PSweights_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/110000/3133B703-B302-6A41-92F2-D546632164A8.root'
'/store/mc/RunIIAutumn18MiniAOD/WkkToWRadionToWWW_M3000-R0-5_TuneCP5_13TeV-madgraph/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/110000/C984DF43-BB3E-0348-A8BE-D0404B8C095A.root'
#'/store/mc/RunIIFall17MiniAODv2/WkkToWRadionToWWW_M5000-R0-3_TuneCP5_13TeV-madgraph/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/40000/C0C22A41-8A4C-E811-BD33-0025907D24F0.root'
]

process.maxEvents.input = 100
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000
process.MessageLogger.cerr.FwkReport.limit = 99999999
print "hh"
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("RStreeEDBR_pickup.root")
                                   )
