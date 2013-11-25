#ifndef HZZ4LeptonsCommonRootTree_h
#define HZZ4LeptonsCommonRootTree_h

/** \class  HZZ4LeptonsCommonRootTree
 *
 *  Root Tree for H->ZZ->4l analysis.
 *
 *  Author: N. De Filippis - Politecnico and INFN Bari
 *          
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/DataKeyTags.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

// data format
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"    
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

// Trigger
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Math/interface/deltaR.h"

// Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

// Maps
#include "DataFormats/Common/interface/ValueMap.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

// user include files
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
class MultiTrajectoryStateMode ;
class EgammaTowerIsolation ;


// Class to create TTree variables
#include <TFile.h> 
#include <TTree.h> 

// Namespaces
using namespace reco;
using namespace std;
using namespace pat;

class HZZ4LeptonsCommonRootTree : public edm::EDAnalyzer {
  
 public:
  
  HZZ4LeptonsCommonRootTree(const edm::ParameterSet& pset);
  
  ~HZZ4LeptonsCommonRootTree();
  
  void analyze(const edm::Event& e, const edm::EventSetup& c);
  void beginJob();
  void endJob();
  
  const EcalRecHit & getRecHit(DetId id, const EcalRecHitCollection *recHits);
  float E2overE9( const DetId , const EcalRecHitCollection & , float , float ,  bool , bool);
  float recHitE( const DetId , const EcalRecHitCollection & );
  float recHitE( const DetId , const EcalRecHitCollection & ,  int , int );
  float recHitApproxEt( const DetId , const EcalRecHitCollection &);
	
  void ReadParameters(const edm::ParameterSet& pset){ 
    std::cout << "Reading parameters from cfg" << std::endl;
    
    typedef std::vector<edm::InputTag> vtag;
    // Get the various input parameters
    decaychannel              = pset.getParameter<std::string>("decaychannel");
    // useBestCandidate          = pset.getParameter<bool> ("useBestCandidate" );
    // BestCandidatesLeptonsTag_ = pset.getParameter<edm::InputTag>("BestCandidatesLeptons");
    rootFileName              = pset.getUntrackedParameter<std::string>("rootFileName");
    useRECOformat             = pset.getUntrackedParameter<bool>("useRECOformat");

    // Get PU simulation info
    fillPUinfo                = pset.getUntrackedParameter<bool>("fillPUinfo");
    PileupSrc_                = pset.getParameter<edm::InputTag>("PileupSrc");

    // Get HLT flags
    fillHLTinfo               = pset.getUntrackedParameter<bool>("fillHLTinfo");
    HLTInfoFired              = pset.getParameter<edm::InputTag>("HLTInfoFired");
    HLTAnalysisinst           = pset.getParameter<string>("HLTAnalysisinst");
    flagHLTnames              = pset.getParameter<vtag>("flagHLTnames");
    // Get HLT matching
    triggerEvent              = pset.getParameter<edm::InputTag>("triggerEvent");  
    triggerFilter             = pset.getParameter<std::string>("triggerFilter");  
    triggerMatchObject        = pset.getParameter<edm::InputTag>("triggerMatchObject");
    triggerFilter_asym        = pset.getParameter<vector<std::string> >("triggerFilterAsym");  
    triggerMatchObject_asym   = pset.getParameter<edm::InputTag>("triggerMatchObject_asym");
    triggerMatchObjectEle     = pset.getParameter<edm::InputTag>("triggerMatchObjectEle");
    triggerHLTcollection      = pset.getParameter<std::string>("triggerHLTcollection");

     // Get SkimEarlyData flags
    useSkimEarlyData          = pset.getUntrackedParameter<bool>("useSkimEarlyData");
    cout << "early " << useSkimEarlyData << endl;
    SkimEarlyDataAnalysisinst = pset.getParameter<string>("SkimEarlyDataAnalysisinst");
    flagSkimEarlyDatanames    = pset.getParameter<vtag>("flagSkimEarlyDatanames");
    // Get flags 
    flaginst                  = pset.getParameter<std::string>("flaginst");
    flagtags                  = pset.getParameter<std::vector<std::string> >("flagtags");
    // MCtruth tags
    fillMCTruth               = pset.getUntrackedParameter<bool>("fillMCTruth");
    MCcollName                = pset.getParameter<vtag>("MCcollName");
    MCcollNameRestFrame       = pset.getParameter<vtag>("MCcollNameRestFrame");
    // RECO tags
    RECOcollNameBest2e2mu     = pset.getParameter<vtag>("RECOcollNameBest2e2mu");
    RECOcollNameBest4mu       = pset.getParameter<vtag>("RECOcollNameBest4mu");
    RECOcollNameBest4e        = pset.getParameter<vtag>("RECOcollNameBest4e");
    fillRF                    = pset.getUntrackedParameter<bool>("fillRF");
    RECOcollNameBestRestFrame2e2mu = pset.getParameter<vtag>("RECOcollNameBestRestFrame2e2mu");
    RECOcollNameBestRestFrame4mu   = pset.getParameter<vtag>("RECOcollNameBestRestFrame4mu");
    RECOcollNameBestRestFrame4e    = pset.getParameter<vtag>("RECOcollNameBestRestFrame4e");

    // RECO additional tags
    useAdditionalRECO         = pset.getUntrackedParameter<bool>("useAdditionalRECO");
    RECOcollNameZ             = pset.getParameter<vtag>("RECOcollNameZ");
    RECOcollNameZss           = pset.getParameter<vtag>("RECOcollNameZss");
    RECOcollNameDiLep         = pset.getParameter<edm::InputTag>("RECOcollNameDiLep");
    //RECOcollNameMMMM          = pset.getParameter<vtag>("RECOcollNameMMMM");
    //RECOcollNameEEEE          = pset.getParameter<vtag>("RECOcollNameEEEE");
    //RECOcollNameEEMM          = pset.getParameter<vtag>("RECOcollNameEEMM");
    RECOcollNameLLLLss        = pset.getParameter<vtag>("RECOcollNameLLLLss");
    RECOcollNameLLLLssos      = pset.getParameter<vtag>("RECOcollNameLLLLssos");
    RECOcollNameLLL           = pset.getParameter<vtag>("RECOcollNameLLL");
    RECOcollNameLLLl          = pset.getParameter<vtag>("RECOcollNameLLLl");
    RECOcollNameLLLL          = pset.getParameter<edm::InputTag>("RECOcollNameLLLL");
    RECOcollNameMMMM          = pset.getParameter<edm::InputTag>("RECOcollNameMMMM");
    RECOcollNameEEEE          = pset.getParameter<edm::InputTag>("RECOcollNameEEEE");
    RECOcollNameEEMM          = pset.getParameter<edm::InputTag>("RECOcollNameEEMM");
    RECOcollNameMMEE          = pset.getParameter<edm::InputTag>("RECOcollNameMMEE");
    RECOcollNameMMMT          = pset.getParameter<edm::InputTag>("RECOcollNameMMMT");
    RECOcollNameEEET          = pset.getParameter<edm::InputTag>("RECOcollNameEEET");
    RECOcollNameMMTT          = pset.getParameter<edm::InputTag>("RECOcollNameMMTT");
    RECOcollNameEETT          = pset.getParameter<edm::InputTag>("RECOcollNameEETT");
    RECOcollNameMMET          = pset.getParameter<edm::InputTag>("RECOcollNameMMET");
    RECOcollNameEEMT          = pset.getParameter<edm::InputTag>("RECOcollNameEEMT");
    RECOcollNameMMME          = pset.getParameter<edm::InputTag>("RECOcollNameMMME");
    RECOcollNameEEEM          = pset.getParameter<edm::InputTag>("RECOcollNameEEEM");
    RECOcollNameMT          = pset.getParameter<edm::InputTag>("RECOcollNameMT");
    RECOcollNameET          = pset.getParameter<edm::InputTag>("RECOcollNameET");
    RECOcollNameME          = pset.getParameter<edm::InputTag>("RECOcollNameME");

    // electrons and muons tags
    if (decaychannel=="2e2mu" || decaychannel=="4mu" || decaychannel=="2l2tau") {
      muonTag_              = pset.getParameter<edm::InputTag>("MuonsLabel");
      muonMapTag_           = pset.getParameter<edm::InputTag>("MuonsMapLabel");
      muonTkMapTag_         = pset.getParameter<edm::InputTag>("MuonsTkMapLabel");
      muonEcalMapTag_       = pset.getParameter<edm::InputTag>("MuonsEcalMapLabel");
      muonHcalMapTag_       = pset.getParameter<edm::InputTag>("MuonsHcalMapLabel");
      
      // Particle Flow Isolation
      muonPFIsoValueAllTag_              = pset.getParameter<edm::InputTag>("MuonsPFIsoValueAll");
      muonPFIsoValueChargedTag_          = pset.getParameter<edm::InputTag>("MuonsPFIsoValueCharged");
      muonPFIsoValueNeutralTag_          = pset.getParameter<edm::InputTag>("MuonsPFIsoValueNeutral");
      muonPFIsoValueGammaTag_            = pset.getParameter<edm::InputTag>("MuonsPFIsoValueGamma");
      muonPFIsoValuePUTag_               = pset.getParameter<edm::InputTag>("MuonsPFIsoValuePU");
      muonPFIsoValuePULowTag_            = pset.getParameter<edm::InputTag>("MuonsPFIsoValuePULow");

      electronPFIsoValueAllTag_              = pset.getParameter<edm::InputTag>("ElectronsPFIsoValueAll");
      electronPFIsoValueChargedTag_          = pset.getParameter<edm::InputTag>("ElectronsPFIsoValueCharged");
      electronPFIsoValueNeutralTag_          = pset.getParameter<edm::InputTag>("ElectronsPFIsoValueNeutral");
      electronPFIsoValueGammaTag_            = pset.getParameter<edm::InputTag>("ElectronsPFIsoValueGamma");
      electronPFIsoValuePUTag_               = pset.getParameter<edm::InputTag>("ElectronsPFIsoValuePU");
      electronPFIsoValuePULowTag_            = pset.getParameter<edm::InputTag>("ElectronsPFIsoValuePULow");
    }
    
    if (decaychannel=="2e2mu" || decaychannel=="4e" || decaychannel=="2l2tau"){ 
      clusterCollectionTag_    = pset.getParameter<edm::InputTag>("SuperClustersLabel");
      gsftrackCollection_      = pset.getParameter<edm::InputTag>("GsfTracksElectronsLabel");
      electronTag_             = pset.getParameter<edm::InputTag>("ElectronsLabel");
      electronMapTag_          = pset.getParameter<edm::InputTag>("ElectronsMapLabel");
      electronTkMapTag_        = pset.getParameter<edm::InputTag>("ElectronsTkMapLabel");
      electronEcalMapTag_      = pset.getParameter<edm::InputTag>("ElectronsEcalMapLabel");
      electronHcalMapTag_      = pset.getParameter<edm::InputTag>("ElectronsHcalMapLabel");
 
      electronEgmTag_          = pset.getParameter<edm::InputTag>("ElectronsEgmLabel");
      electronEgmTkMapTag_     = pset.getParameter<edm::InputTag>("ElectronsEgmTkMapLabel");
      electronEgmEcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmEcalMapLabel");
      electronEgmHcalMapTag_   = pset.getParameter<edm::InputTag>("ElectronsEgmHcalMapLabel");

      EleID_VeryLooseTag_      = pset.getParameter<edm::InputTag>("eleID_VeryLooseTag");
      EleID_LooseTag_          = pset.getParameter<edm::InputTag>("eleID_LooseTag");
      EleID_MediumTag_         = pset.getParameter<edm::InputTag>("eleID_MediumTag");
      EleID_TightTag_          = pset.getParameter<edm::InputTag>("eleID_TightTag");
      
      EleID_HZZVeryLooseTag_   = pset.getParameter<edm::InputTag>("eleID_HZZVeryLooseTag");
      EleID_HZZLooseTag_       = pset.getParameter<edm::InputTag>("eleID_HZZLooseTag");
      EleID_HZZMediumTag_      = pset.getParameter<edm::InputTag>("eleID_HZZMediumTag");
      EleID_HZZTightTag_       = pset.getParameter<edm::InputTag>("eleID_HZZTightTag");
      
    }
    
    // vertexing
    if (decaychannel=="2e2mu" || decaychannel=="4mu" || decaychannel=="2l2tau") {
      // 3D w.r.t primary vertex DA
      muonTag_Vert           = pset.getParameter<edm::InputTag>("MuonsLabelVert");
      muonMapTag_Vert        = pset.getParameter<edm::InputTag>("MuonsMapLabelVert");
      muonMapTag_VertValue   = pset.getParameter<edm::InputTag>("MuonsMapLabelVertValue");
      muonMapTag_VertError   = pset.getParameter<edm::InputTag>("MuonsMapLabelVertError");
      // KF
      muonMapTag_VertKF        = pset.getParameter<edm::InputTag>("MuonsMapLabelVertKF");
      muonMapTag_VertValueKF   = pset.getParameter<edm::InputTag>("MuonsMapLabelVertValueKF");
      muonMapTag_VertErrorKF   = pset.getParameter<edm::InputTag>("MuonsMapLabelVertErrorKF");


      // GD, Std and Kin vertex
      muonMapTag_VertGD      = pset.getParameter<edm::InputTag>("MuonsMapLabelVertGD");
      muonMapTag_VertGDMMMM  = pset.getParameter<edm::InputTag>("MuonsMapLabelVertGDMMMM");
      muonMapTag_VertStd     = pset.getParameter<edm::InputTag>("MuonsMapLabelVertStd");
      muonMapTag_VertStdMMMM = pset.getParameter<edm::InputTag>("MuonsMapLabelVertStdMMMM");
      muonMapTag_VertKin     = pset.getParameter<edm::InputTag>("MuonsMapLabelVertKin");
      muonMapTag_VertKinMMMM = pset.getParameter<edm::InputTag>("MuonsMapLabelVertKinMMMM");
      
      // STIP SLIP
      muonSTIPMapTag_Vert   = pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVert");
      muonSLIPMapTag_Vert   = pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVert");
      
      muonSTIPMapTag_VertValue   = pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertValue");
      muonSLIPMapTag_VertValue   = pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertValue");
      muonSTIPMapTag_VertError   = pset.getParameter<edm::InputTag>("MuonsSTIPMapLabelVertError");
      muonSLIPMapTag_VertError   = pset.getParameter<edm::InputTag>("MuonsSLIPMapLabelVertError");
      
    }
    
    if (decaychannel=="2e2mu" || decaychannel=="4e" || decaychannel=="2l2tau") {
      // 3D w.r.t primary vertex DA
      electronTag_Vert          = pset.getParameter<edm::InputTag>("ElectronsLabelVert");
      electronMapTag_Vert       = pset.getParameter<edm::InputTag>("ElectronsMapLabelVert");
      electronMapTag_VertValue  = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertValue");
      electronMapTag_VertError  = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertError");
      // KF
      electronMapTag_VertKF       = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertKF");
      electronMapTag_VertValueKF  = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertValueKF");
      electronMapTag_VertErrorKF  = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertErrorKF");


      // GD, Std and Kin vertex
      electronMapTag_VertGD      = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertGD");
      electronMapTag_VertGDEEEE  = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertGDEEEE");
      electronMapTag_VertStd     = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertStd");
      electronMapTag_VertStdEEEE = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertStdEEEE");
      electronMapTag_VertKin     = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertKin");
      electronMapTag_VertKinEEEE = pset.getParameter<edm::InputTag>("ElectronsMapLabelVertKinEEEE");

      // STIP SLIP
      electronSTIPMapTag_Vert   = pset.getParameter<edm::InputTag>("ElectronsSTIPMapLabelVert");
      electronSLIPMapTag_Vert   = pset.getParameter<edm::InputTag>("ElectronsSLIPMapLabelVert");
      
      electronSTIPMapTag_VertValue   = pset.getParameter<edm::InputTag>("ElectronsSTIPMapLabelVertValue");
      electronSLIPMapTag_VertValue   = pset.getParameter<edm::InputTag>("ElectronsSLIPMapLabelVertValue");
      electronSTIPMapTag_VertError   = pset.getParameter<edm::InputTag>("ElectronsSTIPMapLabelVertError");
      electronSLIPMapTag_VertError   = pset.getParameter<edm::InputTag>("ElectronsSLIPMapLabelVertError");
      
    }

    if(decaychannel=="2e2mu" || decaychannel=="4e" ){
    // geom. discr.
    ftsigma_Vert           = pset.getParameter<edm::InputTag>("ftsigmaVert");
    ftsigmalag_Vert        = pset.getParameter<edm::InputTag>("ftsigmalagVert");
    gdX_Vert               = pset.getParameter<edm::InputTag>("gdX_Vert");
    gdY_Vert               = pset.getParameter<edm::InputTag>("gdY_Vert");
    gdZ_Vert               = pset.getParameter<edm::InputTag>("gdZ_Vert");
    gdlagX_Vert            = pset.getParameter<edm::InputTag>("gdlagX_Vert");
    gdlagY_Vert            = pset.getParameter<edm::InputTag>("gdlagY_Vert");
    gdlagZ_Vert            = pset.getParameter<edm::InputTag>("gdlagZ_Vert");
    gdlagProb_Vert         = pset.getParameter<edm::InputTag>("gdlagProb_Vert");
    gdlagNdof_Vert         = pset.getParameter<edm::InputTag>("gdlagNdof_Vert");

    ftsigma_VertMMMM       = pset.getParameter<edm::InputTag>("ftsigmaVertMMMM");
    ftsigmalag_VertMMMM    = pset.getParameter<edm::InputTag>("ftsigmalagVertMMMM");
    gdX_VertMMMM           = pset.getParameter<edm::InputTag>("gdX_VertMMMM");
    gdY_VertMMMM           = pset.getParameter<edm::InputTag>("gdY_VertMMMM");
    gdZ_VertMMMM           = pset.getParameter<edm::InputTag>("gdZ_VertMMMM");
    gdlagX_VertMMMM        = pset.getParameter<edm::InputTag>("gdlagX_VertMMMM");
    gdlagY_VertMMMM        = pset.getParameter<edm::InputTag>("gdlagY_VertMMMM");
    gdlagZ_VertMMMM        = pset.getParameter<edm::InputTag>("gdlagZ_VertMMMM");
    gdlagProb_VertMMMM     = pset.getParameter<edm::InputTag>("gdlagProb_VertMMMM");
    gdlagNdof_VertMMMM     = pset.getParameter<edm::InputTag>("gdlagNdof_VertMMMM");

    ftsigma_VertEEEE       = pset.getParameter<edm::InputTag>("ftsigmaVertEEEE");
    ftsigmalag_VertEEEE    = pset.getParameter<edm::InputTag>("ftsigmalagVertEEEE");
    gdX_VertEEEE           = pset.getParameter<edm::InputTag>("gdX_VertEEEE");
    gdY_VertEEEE           = pset.getParameter<edm::InputTag>("gdY_VertEEEE");
    gdZ_VertEEEE           = pset.getParameter<edm::InputTag>("gdZ_VertEEEE");
    gdlagX_VertEEEE        = pset.getParameter<edm::InputTag>("gdlagX_VertEEEE");
    gdlagY_VertEEEE        = pset.getParameter<edm::InputTag>("gdlagY_VertEEEE");
    gdlagZ_VertEEEE        = pset.getParameter<edm::InputTag>("gdlagZ_VertEEEE");
    gdlagProb_VertEEEE     = pset.getParameter<edm::InputTag>("gdlagProb_VertEEEE");
    gdlagNdof_VertEEEE     = pset.getParameter<edm::InputTag>("gdlagNdof_VertEEEE");

    //ConstraintVertex 4l
    StandardFitVertex         = pset.getParameter<edm::InputTag>("StandardFitVertex");
    StandardFitVertexMMMM     = pset.getParameter<edm::InputTag>("StandardFitVertexMMMM");
    StandardFitVertexEEEE     = pset.getParameter<edm::InputTag>("StandardFitVertexEEEE");
    KinematicFitVertex      = pset.getParameter<edm::InputTag>("KinematicFitVertex");
    KinematicFitVertexMMMM  = pset.getParameter<edm::InputTag>("KinematicFitVertexMMMM");
    KinematicFitVertexEEEE  = pset.getParameter<edm::InputTag>("KinematicFitVertexEEEE");

    //ConstraintVertex 3l
    StandardFitVertexMMM     = pset.getParameter<edm::InputTag>("StandardFitVertexMMM");
    StandardFitVertexMME     = pset.getParameter<edm::InputTag>("StandardFitVertexMME");
    StandardFitVertexEEE     = pset.getParameter<edm::InputTag>("StandardFitVertexEEE");
    StandardFitVertexMEE     = pset.getParameter<edm::InputTag>("StandardFitVertexMEE");

    // ConstraintVertex Dileptons
    StandardFitVertexDiLep   = pset.getParameter<edm::InputTag>("StandardFitVertexDiLep");
    }

    //electronID 
    eleIDTag_                   = pset.getParameter<vtag>("eleIDLabel");

    // MCCP parameters
    MCCP_2e2mu_cosTheta1Tag_    = pset.getParameter<edm::InputTag>("MCCP_2e2mu_cosTheta1Label");
    MCCP_2e2mu_cosTheta2Tag_    = pset.getParameter<edm::InputTag>("MCCP_2e2mu_cosTheta2Label");
    MCCP_2e2mu_cosThetaStarTag_ = pset.getParameter<edm::InputTag>("MCCP_2e2mu_cosThetaStarLabel");
    MCCP_2e2mu_PhiTag_          = pset.getParameter<edm::InputTag>("MCCP_2e2mu_PhiLabel");
    MCCP_2e2mu_Phi1Tag_         = pset.getParameter<edm::InputTag>("MCCP_2e2mu_Phi1Label");
    MCCP_2e2mu_Phi2Tag_         = pset.getParameter<edm::InputTag>("MCCP_2e2mu_Phi2Label");
    MCCP_2e2mu_phi1RFTag_       = pset.getParameter<edm::InputTag>("MCCP_2e2mu_phi1RFLabel");
    MCCP_2e2mu_phi2RFTag_       = pset.getParameter<edm::InputTag>("MCCP_2e2mu_phi2RFLabel");
    
    MCCP_4mu_cosTheta1Tag_    = pset.getParameter<edm::InputTag>("MCCP_4mu_cosTheta1Label");
    MCCP_4mu_cosTheta2Tag_    = pset.getParameter<edm::InputTag>("MCCP_4mu_cosTheta2Label");
    MCCP_4mu_cosThetaStarTag_ = pset.getParameter<edm::InputTag>("MCCP_4mu_cosThetaStarLabel");
    MCCP_4mu_PhiTag_          = pset.getParameter<edm::InputTag>("MCCP_4mu_PhiLabel");
    MCCP_4mu_Phi1Tag_         = pset.getParameter<edm::InputTag>("MCCP_4mu_Phi1Label");
    MCCP_4mu_Phi2Tag_         = pset.getParameter<edm::InputTag>("MCCP_4mu_Phi2Label");
    MCCP_4mu_phi1RFTag_       = pset.getParameter<edm::InputTag>("MCCP_4mu_phi1RFLabel");
    MCCP_4mu_phi2RFTag_       = pset.getParameter<edm::InputTag>("MCCP_4mu_phi2RFLabel");
    
    MCCP_4e_cosTheta1Tag_    = pset.getParameter<edm::InputTag>("MCCP_4e_cosTheta1Label");
    MCCP_4e_cosTheta2Tag_    = pset.getParameter<edm::InputTag>("MCCP_4e_cosTheta2Label");
    MCCP_4e_cosThetaStarTag_ = pset.getParameter<edm::InputTag>("MCCP_4e_cosThetaStarLabel");
    MCCP_4e_PhiTag_          = pset.getParameter<edm::InputTag>("MCCP_4e_PhiLabel");
    MCCP_4e_Phi1Tag_         = pset.getParameter<edm::InputTag>("MCCP_4e_Phi1Label");
    MCCP_4e_Phi2Tag_         = pset.getParameter<edm::InputTag>("MCCP_4e_Phi2Label");
    MCCP_4e_phi1RFTag_       = pset.getParameter<edm::InputTag>("MCCP_4e_phi1RFLabel");
    MCCP_4e_phi2RFTag_       = pset.getParameter<edm::InputTag>("MCCP_4e_phi2RFLabel");
    


    // CP parameters
    CP2e2mu_cosTheta1Tag_    = pset.getParameter<edm::InputTag>("CP2e2mu_cosTheta1Label");
    CP2e2mu_cosTheta2Tag_    = pset.getParameter<edm::InputTag>("CP2e2mu_cosTheta2Label");
    CP2e2mu_cosThetaStarTag_ = pset.getParameter<edm::InputTag>("CP2e2mu_cosThetaStarLabel");
    CP2e2mu_PhiTag_       = pset.getParameter<edm::InputTag>("CP2e2mu_PhiLabel");
    CP2e2mu_Phi1Tag_      = pset.getParameter<edm::InputTag>("CP2e2mu_Phi1Label");
    CP2e2mu_Phi2Tag_      = pset.getParameter<edm::InputTag>("CP2e2mu_Phi2Label");
    CP2e2mu_phi1RFTag_      = pset.getParameter<edm::InputTag>("CP2e2mu_phi1RFLabel");
    CP2e2mu_phi2RFTag_      = pset.getParameter<edm::InputTag>("CP2e2mu_phi2RFLabel");

    CP4mu_cosTheta1Tag_    = pset.getParameter<edm::InputTag>("CP4mu_cosTheta1Label");
    CP4mu_cosTheta2Tag_    = pset.getParameter<edm::InputTag>("CP4mu_cosTheta2Label");
    CP4mu_cosThetaStarTag_ = pset.getParameter<edm::InputTag>("CP4mu_cosThetaStarLabel");
    CP4mu_PhiTag_       = pset.getParameter<edm::InputTag>("CP4mu_PhiLabel");
    CP4mu_Phi1Tag_      = pset.getParameter<edm::InputTag>("CP4mu_Phi1Label");
    CP4mu_Phi2Tag_      = pset.getParameter<edm::InputTag>("CP4mu_Phi2Label");
    CP4mu_phi1RFTag_      = pset.getParameter<edm::InputTag>("CP4mu_phi1RFLabel");
    CP4mu_phi2RFTag_      = pset.getParameter<edm::InputTag>("CP4mu_phi2RFLabel");
    
    CP4e_cosTheta1Tag_    = pset.getParameter<edm::InputTag>("CP4e_cosTheta1Label");
    CP4e_cosTheta2Tag_    = pset.getParameter<edm::InputTag>("CP4e_cosTheta2Label");
    CP4e_cosThetaStarTag_ = pset.getParameter<edm::InputTag>("CP4e_cosThetaStarLabel");
    CP4e_PhiTag_       = pset.getParameter<edm::InputTag>("CP4e_PhiLabel");
    CP4e_Phi1Tag_      = pset.getParameter<edm::InputTag>("CP4e_Phi1Label");
    CP4e_Phi2Tag_      = pset.getParameter<edm::InputTag>("CP4e_Phi2Label");
    CP4e_phi1RFTag_      = pset.getParameter<edm::InputTag>("CP4e_phi1RFLabel");
    CP4e_phi2RFTag_      = pset.getParameter<edm::InputTag>("CP4e_phi2RFLabel");


    // Other objets
    photonsTag_      = pset.getParameter<edm::InputTag>("PhotonsLabel");
    tracksTag_       = pset.getParameter<edm::InputTag>("TracksLabel");
    jetsTag_         = pset.getParameter<edm::InputTag>("JetsLabel");
    rhojetsTag_      = pset.getParameter<edm::InputTag>("RhoJetsLabel");
    verticesTag_     = pset.getParameter<edm::InputTag>("VerticesLabel");
    
    // MET reco    
    genmetTag_              = pset.getParameter<edm::InputTag>("GenMETLabel"); // GenMET
    trackermetTag_          = pset.getParameter<edm::InputTag>("TrackerMETLabel"); // TrackerMET
    // Calo MET
    calometTag_             = pset.getParameter<edm::InputTag>("CaloMETLabel"); 
    calometnohfTag_         = pset.getParameter<edm::InputTag>("CaloMET_NoHFLabel");
    useAdditionalMET_       = pset.getUntrackedParameter<bool>("useAdditionalMET"); 
    calomethoTag_           = pset.getParameter<edm::InputTag>("CaloMET_HOLabel");
    calometoptTag_          = pset.getParameter<edm::InputTag>("CaloMET_OptLabel");
    calometoptnohfTag_      = pset.getParameter<edm::InputTag>("CaloMET_OptNoHFLabel");
    calometoptnohfhoTag_    = pset.getParameter<edm::InputTag>("CaloMET_OptNoHFHOLabel");
    calometopthoTag_        = pset.getParameter<edm::InputTag>("CaloMET_OptHOLabel");   
    calometnohfhoTag_       = pset.getParameter<edm::InputTag>("CaloMET_NoHFHOLabel");
    
    // PF MET
    pfmetTag_               = pset.getParameter<edm::InputTag>("PfMETLabel");
    // htMET
    htmetic5Tag_            = pset.getParameter<edm::InputTag>("HtMET_IC5Label");
    htmetkt4Tag_            = pset.getParameter<edm::InputTag>("HtMET_KT4Label");
    htmetkt6Tag_            = pset.getParameter<edm::InputTag>("HtMET_KT6Label");
    htmetsc5Tag_            = pset.getParameter<edm::InputTag>("HtMET_SC5Label");
    htmetsc7Tag_            = pset.getParameter<edm::InputTag>("HtMET_SC7Label");
    // JES correction MET
    jescormetic5Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorIC5CaloJetLabel");
    jescormetkt4Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorKT4CaloJetLabel");
    jescormetkt6Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorKT6CaloJetLabel");
    jescormetsc5Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorSC5CaloJetLabel");
    jescormetsc7Tag_        = pset.getParameter<edm::InputTag>("MET_JESCorSC7CaloJetLabel");
    // Type I Muon Correction MET
    cormetMuTag_            = pset.getParameter<edm::InputTag>("CorMETGlobalMuLabel");
    
    // btagging
    tCHighEff_bTag_         = pset.getParameter<edm::InputTag>("tCHighEff_bTagLabel"); 
    tCHighPur_bTag_         = pset.getParameter<edm::InputTag>("tCHighPur_bTagLabel");
    jPHighEff_bTag_         = pset.getParameter<edm::InputTag>("jPHighEff_bTagLabel");
    jBP_bTag_               = pset.getParameter<edm::InputTag>("jBP_bTagLabel");
    sSVHighEff_bTag_        = pset.getParameter<edm::InputTag>("sSVHighEff_bTagLabel");
    sSVHighPur_bTag_        = pset.getParameter<edm::InputTag>("sSVHighPur_bTagLabel");
    cSV_bTag_               = pset.getParameter<edm::InputTag>("cSV_bTagLabel");
    cSVMVA_bTag_            = pset.getParameter<edm::InputTag>("cSVMVA_bTagLabel");
    sEByIP3d_bTag_          = pset.getParameter<edm::InputTag>("sEByIP3d_bTagLabel");
    sEByPt_bTag_            = pset.getParameter<edm::InputTag>("sEByPt_bTagLabel");
    sM_bTag_                = pset.getParameter<edm::InputTag>("sM_bTagLabel");
    sMByIP3d_bTag_          = pset.getParameter<edm::InputTag>("sMByIP3d_bTagLabel");
    sMByPt_bTag_            = pset.getParameter<edm::InputTag>("sMByPt_bTagLabel");

    // Conversion finder
    ConvMapDistTag_       = pset.getParameter<edm::InputTag>("ConvMapDist");
    ConvMapDcotTag_       = pset.getParameter<edm::InputTag>("ConvMapDcot");

  }
  
  
  void DefineBranches(TTree *Tree_){
    // Run event
    Tree_->Branch("Run",&irun,"irun/i");
    Tree_->Branch("Event",&ievt,"ievt/i");
    Tree_->Branch("LumiSection",&ils,"ils/i");
    
    Tree_->Branch("InternalXsection",&auto_cross_section,"auto_cross_section/D");
    Tree_->Branch("Xsection",&external_cross_section,"external_cross_section/D");
    //Tree_->Branch("GenFilterEff",&filter_eff, "filter_eff/D");
    
    Tree_->Branch("num_PU_vertices",&num_PU_vertices,"num_PU_vertices/I");
    Tree_->Branch("PU_BunchCrossing",&PU_BunchCrossing,"PU_BunchCrossing/I");
    Tree_->Branch("PU_zpos",PU_zpos,"PU_zpos[50]/F");
    Tree_->Branch("PU_sumpT_lowpT",PU_sumpT_lowpT,"PU_sumpT_lowpT[50]/F");
    Tree_->Branch("PU_sumpT_highpT",PU_sumpT_highpT,"PU_sumpT_highpT[50]/F");
      
    Tree_->Branch("PU_ntrks_lowpT",PU_ntrks_lowpT,"PU_ntrks_lowpT[50]/I");
    Tree_->Branch("PU_ntrks_highpT",PU_ntrks_highpT,"PU_ntrks_highpT[50]/I");

    // HLT flags
    // Tree_->Branch("boolMuonIso",&boolMuonIso,"boolMuonIso/b"); */
    // Tree_->Branch("boolMuonNonIso",&boolMuonNonIso,"boolMuonNonIso/b"); */
    // Tree_->Branch("bool2MuonNonIso",&bool2MuonNonIso,"bool2MuonNonIso/b"); */
    // Tree_->Branch("boolElectron",&boolElectron,"boolElectron/b"); */
    // Tree_->Branch("boolElectronRelaxed",&boolElectronRelaxed,"boolElectronRelaxed/b"); */
    // Tree_->Branch("bool2Electron",&bool2Electron,"bool2Electron/b"); */
    // Tree_->Branch("bool2ElectronRelaxed",&bool2ElectronRelaxed,"bool2ElectronRelaxed/b"); */
    // Tree_->Branch("boolCandHLT2ElectronStartup", &boolCandHLT2ElectronStartup,"boolCandHLT2ElectronStartup/b");
    // Tree_->Branch("boolXElectronMuon",&boolXElectronMuon,"boolXElectronMuon/b");  
    // Tree_->Branch("boolHLTaccept",&boolHLTaccept,"boolHLTaccept/b");

    Tree_->Branch("boolHLTpattern",&boolHLTpattern,"boolHLTpattern[50]/b");
    Tree_->Branch("RECO_nMuHLTMatch",&RECO_nMuHLTMatch,"RECO_nMuHLTMatch/I");
    Tree_->Branch("RECO_nMuHLTMatchPAT",&RECO_nMuHLTMatchPAT,"RECO_nMuHLTMatchPAT/I");
    Tree_->Branch("RECO_nMuHLTMatch_asym_PAT",&RECO_nMuHLTMatch_asym_PAT,"RECO_nMuHLTMatch_asym_PAT/I");
    Tree_->Branch("RECOMU_PT_MuHLTMatch",RECOMU_PT_MuHLTMatch,"RECOMU_PT_MuHLTMatch[100]/F");
    Tree_->Branch("RECOMU_PT_MuHLTMatchPAT",RECOMU_PT_MuHLTMatchPAT,"RECOMU_PT_MuHLTMatchPAT[100]/F");
    Tree_->Branch("RECOMU_PT_MuHLTMatch_asym_PAT",RECOMU_PT_MuHLTMatch_asym_PAT,"RECOMU_PT_MuHLTMatch_asym_PAT[100]/F");

    Tree_->Branch("RECO_nEleHLTMatchPAT",&RECO_nEleHLTMatchPAT,"RECO_nEleHLTMatchPAT/I");
    Tree_->Branch("RECOELE_PT_EleHLTMatchPAT",RECOELE_PT_EleHLTMatchPAT,"RECOELE_PT_EleHLTMatchPAT[100]/F");

    Tree_->Branch("HLTPathsFired",HLTPathsFired,"HLTPathsFired/C");

/*     // Skim Early Data */
/*     Tree_->Branch("boolSkim_highEnergyMuons",&boolSkim_highEnergyMuons,"boolSkim_highEnergyMuons/b"); */
/*     Tree_->Branch("boolSkim_highEnergyElectrons",&boolSkim_highEnergyElectrons,"boolSkim_highEnergyElectrons/b"); */
/*     Tree_->Branch("boolSkim_recoWMNfromPf",&boolSkim_recoWMNfromPf,"boolSkim_recoWMNfromPf/b"); */
/*     Tree_->Branch("boolSkim_recoWMNfromTc",&boolSkim_recoWMNfromTc,"boolSkim_recoWMNfromTc/b"); */
/*     Tree_->Branch("boolSkim_recoWENfromPf",&boolSkim_recoWENfromPf,"boolSkim_recoWENfromPf/b"); */
/*     Tree_->Branch("boolSkim_recoWENfromTc",&boolSkim_recoWENfromTc,"boolSkim_recoWENfromTc/b"); */
/*     Tree_->Branch("boolSkim_diMuonsJPsi",&boolSkim_diMuonsJPsi,"boolSkim_diMuonsJPsi/b"); */
/*     Tree_->Branch("boolSkim_diMuonsZ",&boolSkim_diMuonsZ,"boolSkim_diMuonsZ/b"); */
/*     Tree_->Branch("boolSkim_diElectronsZ",&boolSkim_diElectronsZ,"boolSkim_diElectronsZ/b"); */
/*     Tree_->Branch("boolSkim_triLeptonsMuMuMu",&boolSkim_triLeptonsMuMuMu,"boolSkim_triLeptonsMuMuMu/b"); */
/*     Tree_->Branch("boolSkim_triLeptonsMuMuEl",&boolSkim_triLeptonsMuMuEl,"boolSkim_triLeptonsMuMuEl/b"); */
/*     Tree_->Branch("boolSkim_triLeptonsMuElEl",&boolSkim_triLeptonsMuElEl,"boolSkim_triLeptonsMuElEl/b"); */
/*     Tree_->Branch("boolSkim_triLeptonsElElEl",&boolSkim_triLeptonsElElEl,"boolSkim_triLeptonsElElEl/b"); */
/*     Tree_->Branch("boolSkim_quadLeptons4Mu",  &boolSkim_quadLeptons4Mu,"boolSkim_quadLeptons4Mu/b"); */
/*     Tree_->Branch("boolSkim_quadLeptons2Mu2El",&boolSkim_quadLeptons2Mu2El,"boolSkim_quadLeptons2Mu2El/b"); */
/*     Tree_->Branch("boolSkim_quadLeptons4El",&boolSkim_quadLeptons4El,"boolSkim_quadLeptons4El/b"); */


    // Flags
    Tree_->Branch("flag",&flag,"flag[10]/i");
    
    /* // MC block 2e2mu */
/*     Tree_->Branch("MC_2e2mu_E", MC_2e2mu_E, "MC_2e2mu_E[7]/F");  */
/*     Tree_->Branch("MC_2e2mu_PT", MC_2e2mu_PT, "MC_2e2mu_PT[7]/F");  */
/*     Tree_->Branch("MC_2e2mu_ETA", MC_2e2mu_ETA, "MC_2e2mu_ETA[7]/F");  */
/*     Tree_->Branch("MC_2e2mu_THETA", MC_2e2mu_THETA, "MC_2e2mu_THETA[7]/F"); */
/*     Tree_->Branch("MC_2e2mu_PHI", MC_2e2mu_PHI, "MC_2e2mu_PHI[7]/F"); */
/*     Tree_->Branch("MC_2e2mu_MASS", MC_2e2mu_MASS, "MC_2e2mu_MASS[7]/F"); */
/*     Tree_->Branch("MC_2e2mu_PDGID", &MC_2e2mu_PDGID, "MC_2e2mu_PDGID[7]/F"); */
    
/*     Tree_->Branch("MCRF_2e2mu_E", MCRF_2e2mu_E, "MCRF_2e2mu_E[7]/F");  */
/*     Tree_->Branch("MCRF_2e2mu_PT_", MCRF_2e2mu_PT, "MCRF_2e2mu_PT[7]/F");  */
/*     Tree_->Branch("MCRF_2e2mu_ETA", MCRF_2e2mu_ETA, "MCRF_2e2mu_ETA[7]/F");  */
/*     Tree_->Branch("MCRF_2e2mu_THETA", MCRF_2e2mu_THETA, "MCRF_2e2mu_THETA[7]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_PHI", MCRF_2e2mu_PHI, "MCRF_2e2mu_PHI[7]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_MASS", MCRF_2e2mu_MASS, "MCRF_2e2mu_MASS[7]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_PDGID", &MCRF_2e2mu_PDGID, "MCRF_2e2mu_PDGID[7]/F"); */
    

/*     // MCCP variables */
/*     Tree_->Branch("MCRF_2e2mu_cosTheta1_spin",&MCRF_2e2mu_cosTheta1_spin,"MCRF_2e2mu_cosTheta1_spin[4]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_cosTheta2_spin",&MCRF_2e2mu_cosTheta2_spin,"MCRF_2e2mu_cosTheta2_spin[4]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_cosThetaStar_spin",&MCRF_2e2mu_cosThetaStar_spin,"MCRF_2e2mu_cosThetaStar_spin[4]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_Phi_spin",&MCRF_2e2mu_Phi_spin,"MCRF_2e2mu_Phi_spin[4]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_Phi1_spin",&MCRF_2e2mu_Phi1_spin,"MCRF_2e2mu_Phi1_spin[4]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_Phi2_spin",&MCRF_2e2mu_Phi2_spin,"MCRF_2e2mu_Phi2_spin[4]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_phi1RF_spin",&MCRF_2e2mu_phi1RF_spin,"MCRF_2e2mu_phi1RF_spin[4]/F"); */
/*     Tree_->Branch("MCRF_2e2mu_phi2RF_spin",&MCRF_2e2mu_phi2RF_spin,"MCRF_2e2mu_phi2RF_spin[4]/F"); */

/*     // MC block 4mu */
/*     Tree_->Branch("MC_4mu_E", MC_4mu_E, "MC_4mu_E[7]/F");  */
/*     Tree_->Branch("MC_4mu_PT", MC_4mu_PT, "MC_4mu_PT[7]/F");  */
/*     Tree_->Branch("MC_4mu_ETA", MC_4mu_ETA, "MC_4mu_ETA[7]/F");  */
/*     Tree_->Branch("MC_4mu_THETA", MC_4mu_THETA, "MC_4mu_THETA[7]/F"); */
/*     Tree_->Branch("MC_4mu_PHI", MC_4mu_PHI, "MC_4mu_PHI[7]/F"); */
/*     Tree_->Branch("MC_4mu_MASS", MC_4mu_MASS, "MC_4mu_MASS[7]/F"); */
/*     Tree_->Branch("MC_4mu_PDGID", &MC_4mu_PDGID, "MC_4mu_PDGID[7]/F"); */
    
/*     Tree_->Branch("MCRF_4mu_E", &MCRF_4mu_E, "MCRF_4mu_E[7]/F");  */
/*     Tree_->Branch("MCRF_4mu_PT_", MCRF_4mu_PT, "MCRF_4mu_PT[7]/F");  */
/*     Tree_->Branch("MCRF_4mu_ETA", MCRF_4mu_ETA, "MCRF_4mu_ETA[7]/F");  */
/*     Tree_->Branch("MCRF_4mu_THETA", MCRF_4mu_THETA, "MCRF_4mu_THETA[7]/F"); */
/*     Tree_->Branch("MCRF_4mu_PHI", MCRF_4mu_PHI, "MCRF_4mu_PHI[7]/F"); */
/*     Tree_->Branch("MCRF_4mu_MASS", MCRF_4mu_MASS, "MCRF_4mu_MASS[7]/F"); */
/*     Tree_->Branch("MCRF_4mu_PDGID", &MCRF_4mu_PDGID, "MCRF_4mu_PDGID[7]/F"); */
    

/*     // MCCP variables */
/*     Tree_->Branch("MCRF_4mu_cosTheta1_spin",&MCRF_4mu_cosTheta1_spin,"MCRF_4mu_cosTheta1_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4mu_cosTheta2_spin",&MCRF_4mu_cosTheta2_spin,"MCRF_4mu_cosTheta2_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4mu_cosThetaStar_spin",&MCRF_4mu_cosThetaStar_spin,"MCRF_4mu_cosThetaStar_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4mu_Phi_spin",&MCRF_4mu_Phi_spin,"MCRF_4mu_Phi_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4mu_Phi1_spin",&MCRF_4mu_Phi1_spin,"MCRF_4mu_Phi1_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4mu_Phi2_spin",&MCRF_4mu_Phi2_spin,"MCRF_4mu_Phi2_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4mu_phi1RF_spin",&MCRF_4mu_phi1RF_spin,"MCRF_4mu_phi1RF_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4mu_phi2RF_spin",&MCRF_4mu_phi2RF_spin,"MCRF_4mu_phi2RF_spin[4]/F"); */


/*     // MC block 4e */
/*     Tree_->Branch("MC_4e_E", MC_4e_E, "MC_4e_E[7]/F");  */
/*     Tree_->Branch("MC_4e_PT", MC_4e_PT, "MC_4e_PT[7]/F");  */
/*     Tree_->Branch("MC_4e_ETA", MC_4e_ETA, "MC_4e_ETA[7]/F");  */
/*     Tree_->Branch("MC_4e_THETA", MC_4e_THETA, "MC_4e_THETA[7]/F"); */
/*     Tree_->Branch("MC_4e_PHI", MC_4e_PHI, "MC_4e_PHI[7]/F"); */
/*     Tree_->Branch("MC_4e_MASS", MC_4e_MASS, "MC_4e_MASS[7]/F"); */
/*     Tree_->Branch("MC_4e_PDGID", &MC_4e_PDGID, "MC_4e_PDGID[7]/F"); */
    
/*     Tree_->Branch("MCRF_4e_E", &MCRF_4e_E, "MCRF_4e_E[7]/F");  */
/*     Tree_->Branch("MCRF_4e_PT_", MCRF_4e_PT, "MCRF_4e_PT[7]/F");  */
/*     Tree_->Branch("MCRF_4e_ETA", MCRF_4e_ETA, "MCRF_4e_ETA[7]/F");  */
/*     Tree_->Branch("MCRF_4e_THETA", MCRF_4e_THETA, "MCRF_4e_THETA[7]/F"); */
/*     Tree_->Branch("MCRF_4e_PHI", MCRF_4e_PHI, "MCRF_4e_PHI[7]/F"); */
/*     Tree_->Branch("MCRF_4e_MASS", MCRF_4e_MASS, "MCRF_4e_MASS[7]/F"); */
/*     Tree_->Branch("MCRF_4e_PDGID", &MCRF_4e_PDGID, "MCRF_4e_PDGID[7]/F"); */
    

/*     // MCCP variables */
/*     Tree_->Branch("MCRF_4e_cosTheta1_spin",&MCRF_4e_cosTheta1_spin,"MCRF_4e_cosTheta1_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4e_cosTheta2_spin",&MCRF_4e_cosTheta2_spin,"MCRF_4e_cosTheta2_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4e_cosThetaStar_spin",&MCRF_4e_cosThetaStar_spin,"MCRF_4e_cosThetaStar_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4e_Phi_spin",&MCRF_4e_Phi_spin,"MCRF_4e_Phi_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4e_Phi1_spin",&MCRF_4e_Phi1_spin,"MCRF_4e_Phi1_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4e_Phi2_spin",&MCRF_4e_Phi2_spin,"MCRF_4e_Phi2_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4e_phi1RF_spin",&MCRF_4e_phi1RF_spin,"MCRF_4e_phi1RF_spin[4]/F"); */
/*     Tree_->Branch("MCRF_4e_phi2RF_spin",&MCRF_4e_phi2RF_spin,"MCRF_4e_phi2RF_spin[4]/F"); */
    
/*     Tree_->Branch("MC_GENMET", &genmet, "MC_GENMET/I"); */
    
/*     // RECO block for reconstructed higgs, Z and their daughters */
/*     Tree_->Branch("RECOBEST_2e2mu_PT", RECOBEST_2e2mu_PT, "RECOBEST_2e2mu_PT[7]/F");  */
/*     Tree_->Branch("RECOBEST_2e2mu_ORDERED_PT", RECOBEST_2e2mu_ORDERED_PT, "RECOBEST_2e2mu_ORDERED_PT[7]/F");  */
/*     Tree_->Branch("RECOBEST_2e2mu_ETA", RECOBEST_2e2mu_ETA, "RECOBEST_2e2mu_ETA[7]/F");  */
/*     Tree_->Branch("RECOBEST_2e2mu_THETA", RECOBEST_2e2mu_THETA, "RECOBEST_2e2mu_THETA[7]/F"); */
/*     Tree_->Branch("RECOBEST_2e2mu_PHI", RECOBEST_2e2mu_PHI, "RECOBEST_2e2mu_PHI[7]/F"); */
/*     Tree_->Branch("RECOBEST_2e2mu_MASS", RECOBEST_2e2mu_MASS, "RECOBEST_2e2mu_MASS[7]/F");  */
/*     Tree_->Branch("RECOBEST_2e2mu_isGlobalMu", RECOBEST_2e2mu_isGlobalMu, "RECOBEST_2e2mu_isGlobalMu[7]/b"); */
/*     Tree_->Branch("RECOBEST_2e2mu_isStandAloneMu", RECOBEST_2e2mu_isStandAloneMu, "RECOBEST_2e2mu_isStandAloneMu[7]/b"); */
/*     Tree_->Branch("RECOBEST_2e2mu_isTrackerMu", RECOBEST_2e2mu_isTrackerMu, "RECOBEST_2e2mu_isTrackerMu[7]/b"); */
/*     Tree_->Branch("RECOBEST_2e2mu_isCaloMu", RECOBEST_2e2mu_isCaloMu, "RECOBEST_2e2mu_isCaloMu[7]/b"); */
/*     Tree_->Branch("RECOBEST_2e2mu_isElectron", RECOBEST_2e2mu_isElectron, "RECOBEST_2e2mu_isElectron[7]/b"); */

/*     Tree_->Branch("RECOBEST_4mu_PT", RECOBEST_4mu_PT, "RECOBEST_4mu_PT[7]/F");  */
/*     Tree_->Branch("RECOBEST_4mu_ORDERED_PT", RECOBEST_4mu_ORDERED_PT, "RECOBEST_4mu_ORDERED_PT[7]/F");  */
/*     Tree_->Branch("RECOBEST_4mu_ETA", RECOBEST_4mu_ETA, "RECOBEST_4mu_ETA[7]/F");  */
/*     Tree_->Branch("RECOBEST_4mu_THETA", RECOBEST_4mu_THETA, "RECOBEST_4mu_THETA[7]/F"); */
/*     Tree_->Branch("RECOBEST_4mu_PHI", RECOBEST_4mu_PHI, "RECOBEST_4mu_PHI[7]/F"); */
/*     Tree_->Branch("RECOBEST_4mu_MASS", RECOBEST_4mu_MASS, "RECOBEST_4mu_MASS[7]/F"); */
/*     Tree_->Branch("RECOBEST_4mu_isGlobalMu", RECOBEST_4mu_isGlobalMu, "RECOBEST_4mu_isGlobalMu[7]/b"); */
/*     Tree_->Branch("RECOBEST_4mu_isStandAloneMu", RECOBEST_4mu_isStandAloneMu, "RECOBEST_4mu_isStandAloneMu[7]/b"); */
/*     Tree_->Branch("RECOBEST_4mu_isTrackerMu", RECOBEST_4mu_isTrackerMu, "RECOBEST_4mu_isTrackerMu[7]/b"); */
/*     Tree_->Branch("RECOBEST_4mu_isCaloMu", RECOBEST_4mu_isCaloMu, "RECOBEST_4mu_isCaloMu[7]/b"); */


/*     Tree_->Branch("RECOBEST_4e_PT", RECOBEST_4e_PT, "RECOBEST_4e_PT[7]/F");  */
/*     Tree_->Branch("RECOBEST_4e_ORDERED_PT", RECOBEST_4e_ORDERED_PT, "RECOBEST_4e_ORDERED_PT[7]/F");  */
/*     Tree_->Branch("RECOBEST_4e_ETA", RECOBEST_4e_ETA, "RECOBEST_4e_ETA[7]/F");  */
/*     Tree_->Branch("RECOBEST_4e_THETA", RECOBEST_4e_THETA, "RECOBEST_4e_THETA[7]/F"); */
/*     Tree_->Branch("RECOBEST_4e_PHI", RECOBEST_4e_PHI, "RECOBEST_4e_PHI[7]/F"); */
/*     Tree_->Branch("RECOBEST_4e_MASS", RECOBEST_4e_MASS, "RECOBEST_4e_MASS[7]/F"); */
/*     Tree_->Branch("RECOBEST_4e_isElectron", RECOBEST_4e_isElectron, "RECOBEST_4e_isElectron[7]/b"); */

/*      // RECORF block 2e2mu */
/*     Tree_->Branch("RECORFBEST_2e2mu_PT", RECORFBEST_2e2mu_PT, "RECORFBEST_2e2mu_PT[7]/F");  */
/*     Tree_->Branch("RECORFBEST_2e2mu_ETA", RECORFBEST_2e2mu_ETA, "RECORFBEST_2e2mu_ETA[7]/F");  */
/*     Tree_->Branch("RECORFBEST_2e2mu_THETA", RECORFBEST_2e2mu_THETA, "RECORFBEST_2e2mu_THETA[7]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_PHI", RECORFBEST_2e2mu_PHI, "RECORFBEST_2e2mu_PHI[7]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_MASS", RECORFBEST_2e2mu_MASS, "RECORFBEST_2e2mu_MASS[7]/F"); */
    
/*     Tree_->Branch("RECORFBEST_2e2mu_cosTheta1_spin",&RECORFBEST_2e2mu_cosTheta1_spin,"RECORFBEST_2e2mu_cosTheta1_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_cosTheta2_spin",&RECORFBEST_2e2mu_cosTheta2_spin,"RECORFBEST_2e2mu_cosTheta2_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_cosThetaStar_spin",&RECORFBEST_2e2mu_cosThetaStar_spin,"RECORFBEST_2e2mu_cosThetaStar_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_Phi_spin",&RECORFBEST_2e2mu_Phi_spin,"RECORFBEST_2e2mu_Phi_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_Phi1_spin",&RECORFBEST_2e2mu_Phi1_spin,"RECORFBEST_2e2mu_Phi1_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_Phi2_spin",&RECORFBEST_2e2mu_Phi2_spin,"RECORFBEST_2e2mu_Phi2_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_phi1RF_spin",&RECORFBEST_2e2mu_phi1RF_spin,"RECORFBEST_2e2mu_phi1RF_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_2e2mu_phi2RF_spin",&RECORFBEST_2e2mu_phi2RF_spin,"RECORFBEST_2e2mu_phi2RF_spin[4]/F"); */

/*     Tree_->Branch("RECORFBEST_4e_PT", RECORFBEST_4e_PT, "RECORFBEST_4e_PT[7]/F");  */
/*     Tree_->Branch("RECORFBEST_4e_ETA", RECORFBEST_4e_ETA, "RECORFBEST_4e_ETA[7]/F");  */
/*     Tree_->Branch("RECORFBEST_4e_THETA", RECORFBEST_4e_THETA, "RECORFBEST_4e_THETA[7]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_PHI", RECORFBEST_4e_PHI, "RECORFBEST_4e_PHI[7]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_MASS", RECORFBEST_4e_MASS, "RECORFBEST_4e_MASS[7]/F"); */
    
/*     Tree_->Branch("RECORFBEST_4e_cosTheta1_spin",&RECORFBEST_4e_cosTheta1_spin,"RECORFBEST_4e_cosTheta1_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_cosTheta2_spin",&RECORFBEST_4e_cosTheta2_spin,"RECORFBEST_4e_cosTheta2_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_cosThetaStar_spin",&RECORFBEST_4e_cosThetaStar_spin,"RECORFBEST_4e_cosThetaStar_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_Phi_spin",&RECORFBEST_4e_Phi_spin,"RECORFBEST_4e_Phi_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_Phi1_spin",&RECORFBEST_4e_Phi1_spin,"RECORFBEST_4e_Phi1_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_Phi2_spin",&RECORFBEST_4e_Phi2_spin,"RECORFBEST_4e_Phi2_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_phi1RF_spin",&RECORFBEST_4e_phi1RF_spin,"RECORFBEST_4e_phi1RF_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4e_phi2RF_spin",&RECORFBEST_4e_phi2RF_spin,"RECORFBEST_4e_phi2RF_spin[4]/F"); */

    
/*     Tree_->Branch("RECORFBEST_4mu_PT", RECORFBEST_4mu_PT, "RECORFBEST_4mu_PT[7]/F");  */
/*     Tree_->Branch("RECORFBEST_4mu_ETA", RECORFBEST_4mu_ETA, "RECORFBEST_4mu_ETA[7]/F");  */
/*     Tree_->Branch("RECORFBEST_4mu_THETA", RECORFBEST_4mu_THETA, "RECORFBEST_4mu_THETA[7]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_PHI", RECORFBEST_4mu_PHI, "RECORFBEST_4mu_PHI[7]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_MASS", RECORFBEST_4mu_MASS, "RECORFBEST_4mu_MASS[7]/F"); */
    
/*     Tree_->Branch("RECORFBEST_4mu_cosTheta1_spin",&RECORFBEST_4mu_cosTheta1_spin,"RECORFBEST_4mu_cosTheta1_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_cosTheta2_spin",&RECORFBEST_4mu_cosTheta2_spin,"RECORFBEST_4mu_cosTheta2_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_cosThetaStar_spin",&RECORFBEST_4mu_cosThetaStar_spin,"RECORFBEST_4mu_cosThetaStar_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_Phi_spin",&RECORFBEST_4mu_Phi_spin,"RECORFBEST_4mu_Phi_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_Phi1_spin",&RECORFBEST_4mu_Phi1_spin,"RECORFBEST_4mu_Phi1_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_Phi2_spin",&RECORFBEST_4mu_Phi2_spin,"RECORFBEST_4mu_Phi2_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_phi1RF_spin",&RECORFBEST_4mu_phi1RF_spin,"RECORFBEST_4mu_phi1RF_spin[4]/F"); */
/*     Tree_->Branch("RECORFBEST_4mu_phi2RF_spin",&RECORFBEST_4mu_phi2RF_spin,"RECORFBEST_4mu_phi2RF_spin[4]/F"); */





    // RECO additional block for reconstructed higgs, Z and their daughters
    Tree_->Branch("RECO_ZMM_MASS", RECO_ZMM_MASS, "RECO_ZMM_MASS[50]/F");
    Tree_->Branch("RECO_ZEE_MASS", RECO_ZEE_MASS, "RECO_ZEE_MASS[50]/F");
    //Tree_->Branch("RECO_DiLep_MASS", RECO_DiLep_MASS, "RECO_DiLep_MASS[50]/F");
    Tree_->Branch("RECO_ZMM_PT", RECO_ZMM_PT, "RECO_ZMM_PT[3][50]/F");
    Tree_->Branch("RECO_ZEE_PT", RECO_ZEE_PT, "RECO_ZEE_PT[3][50]/F");
    //Tree_->Branch("RECO_DiLep_PT", RECO_DiLep_PT, "RECO_DiLep_PT[3][50]/F");
    Tree_->Branch("RECO_ZMM_ETA", RECO_ZMM_ETA, "RECO_ZMM_ETA[3][50]/F");
    Tree_->Branch("RECO_ZEE_ETA", RECO_ZEE_ETA, "RECO_ZEE_ETA[3][50]/F");
    //Tree_->Branch("RECO_DiLep_ETA", RECO_DiLep_ETA, "RECO_DiLep_ETA[3][50]/F");
    Tree_->Branch("RECO_ZMM_PHI", RECO_ZMM_PHI, "RECO_ZMM_PHI[3][50]/F");
    Tree_->Branch("RECO_ZEE_PHI", RECO_ZEE_PHI, "RECO_ZEE_PHI[3][50]/F");
    //Tree_->Branch("RECO_DiLep_PHI", RECO_DiLep_PHI, "RECO_DiLep_PHI[3][50]/F");
                              
/*     Tree_->Branch("RECO_ZMMss_MASS", RECO_ZMMss_MASS, "RECO_ZMMss_MASS[50]/F"); */
/*     Tree_->Branch("RECO_ZEEss_MASS", RECO_ZEEss_MASS, "RECO_ZEEss_MASS[50]/F"); */
/*     Tree_->Branch("RECO_ZEM_MASS", RECO_ZEM_MASS, "RECO_ZEM_MASS[50]/F"); */
/*     Tree_->Branch("RECO_ZMMss_PT", RECO_ZMMss_PT, "RECO_ZMMss_PT[3][50]/F"); */
/*     Tree_->Branch("RECO_ZEEss_PT", RECO_ZEEss_PT, "RECO_ZEEss_PT[3][50]/F"); */
/*     Tree_->Branch("RECO_ZEM_PT", RECO_ZEM_PT, "RECO_ZEM_PT[3][50]/F"); */
/*     Tree_->Branch("RECO_ZMMss_ETA", RECO_ZMMss_ETA, "RECO_ZMMss_ETA[3][50]/F"); */
/*     Tree_->Branch("RECO_ZEEss_ETA", RECO_ZEEss_ETA, "RECO_ZEEss_ETA[3][50]/F"); */
/*     Tree_->Branch("RECO_ZEM_ETA", RECO_ZEM_ETA, "RECO_ZEM_ETA[3][50]/F"); */
/*     Tree_->Branch("RECO_ZMMss_PHI", RECO_ZMMss_PHI, "RECO_ZMMss_PHI[3][50]/F"); */
/*     Tree_->Branch("RECO_ZEEss_PHI", RECO_ZEEss_PHI, "RECO_ZEEss_PHI[3][50]/F"); */
/*     Tree_->Branch("RECO_ZEM_PHI", RECO_ZEM_PHI, "RECO_ZEM_PHI[3][50]/F"); */

    
/*     Tree_->Branch("RECO_MMMM_MASS", RECO_MMMM_MASS, "RECO_MMMM_MASS[7][50]/F"); */
/*     Tree_->Branch("RECO_MMMM_PT", RECO_MMMM_PT, "RECO_MMMM_PT[7][50]/F"); */
/*     Tree_->Branch("RECO_MMMM_ETA", RECO_MMMM_ETA, "RECO_MMMM_ETA[7][50]/F"); */
/*     Tree_->Branch("RECO_MMMM_PHI", RECO_MMMM_PHI, "RECO_MMMM_PHI[7][50]/F"); */

/*     Tree_->Branch("RECO_EEEE_MASS", RECO_EEEE_MASS, "RECO_EEEE_MASS[7][50]/F"); */
/*     Tree_->Branch("RECO_EEEE_PT", RECO_EEEE_PT, "RECO_EEEE_PT[7][50]/F");  */
/*     Tree_->Branch("RECO_EEEE_ETA", RECO_EEEE_ETA, "RECO_EEEE_ETA[7][50]/F"); */
/*     Tree_->Branch("RECO_EEEE_PHI", RECO_EEEE_PHI, "RECO_EEEE_PHI[7][50]/F"); */

/*     Tree_->Branch("RECO_EEMM_MASS", RECO_EEMM_MASS, "RECO_EEMM_MASS[7][50]/F"); */
/*     Tree_->Branch("RECO_EEMM_PT", RECO_EEMM_PT, "RECO_EEMM_PT[7][50]/F"); */
/*     Tree_->Branch("RECO_EEMM_ETA", RECO_EEMM_ETA, "RECO_EEMM_ETA[7][50]/F"); */
/*     Tree_->Branch("RECO_EEMM_PHI", RECO_EEMM_PHI, "RECO_EEMM_PHI[7][50]/F"); */

//MMMT
    Tree_->Branch("RECO_MMMT_MASS", RECO_MMMT_MASS, "RECO_MMMT_MASS[9][100]/F");
    Tree_->Branch("RECO_MMMT_PT", RECO_MMMT_PT, "RECO_MMMT_PT[9][100]/F");
    Tree_->Branch("RECO_MMMT_ETA", RECO_MMMT_ETA, "RECO_MMMT_ETA[9][100]/F");
    Tree_->Branch("RECO_MMMT_PHI", RECO_MMMT_PHI, "RECO_MMMT_PHI[9][100]/F");
    Tree_->Branch("RECO_MMMT_CHARGE", RECO_MMMT_CHARGE, "RECO_MMMT_CHARGE[9][100]/F");

    //EEET
    Tree_->Branch("RECO_EEET_MASS", RECO_EEET_MASS, "RECO_EEET_MASS[9][100]/F");
    Tree_->Branch("RECO_EEET_PT", RECO_EEET_PT, "RECO_EEET_PT[9][100]/F");
    Tree_->Branch("RECO_EEET_ETA", RECO_EEET_ETA, "RECO_EEET_ETA[9][100]/F");
    Tree_->Branch("RECO_EEET_PHI", RECO_EEET_PHI, "RECO_EEET_PHI[9][100]/F");
    Tree_->Branch("RECO_EEET_CHARGE", RECO_EEET_CHARGE, "RECO_EEET_CHARGE[9][100]/F");

    //MMTT
    Tree_->Branch("RECO_MMTT_MASS", RECO_MMTT_MASS, "RECO_MMTT_MASS[9][200]/F");
    Tree_->Branch("RECO_MMTT_PT", RECO_MMTT_PT, "RECO_MMTT_PT[9][200]/F");
    Tree_->Branch("RECO_MMTT_ETA", RECO_MMTT_ETA, "RECO_MMTT_ETA[9][200]/F");
    Tree_->Branch("RECO_MMTT_PHI", RECO_MMTT_PHI, "RECO_MMTT_PHI[9][200]/F");
    Tree_->Branch("RECO_MMTT_CHARGE", RECO_MMTT_CHARGE, "RECO_MMTT_CHARGE[9][200]/F");

    //EETT
    Tree_->Branch("RECO_EETT_MASS", RECO_EETT_MASS, "RECO_EETT_MASS[9][200]/F");
    Tree_->Branch("RECO_EETT_PT", RECO_EETT_PT, "RECO_EETT_PT[9][200]/F");
    Tree_->Branch("RECO_EETT_ETA", RECO_EETT_ETA, "RECO_EETT_ETA[9][200]/F");
    Tree_->Branch("RECO_EETT_PHI", RECO_EETT_PHI, "RECO_EETT_PHI[9][200]/F");
    Tree_->Branch("RECO_EETT_CHARGE", RECO_EETT_CHARGE, "RECO_EETT_CHARGE[9][200]/F");

    //MMET
    Tree_->Branch("RECO_MMET_MASS", RECO_MMET_MASS, "RECO_MMET_MASS[9][100]/F");
    Tree_->Branch("RECO_MMET_PT", RECO_MMET_PT, "RECO_MMET_PT[9][100]/F");
    Tree_->Branch("RECO_MMET_ETA", RECO_MMET_ETA, "RECO_MMET_ETA[9][100]/F");
    Tree_->Branch("RECO_MMET_PHI", RECO_MMET_PHI, "RECO_MMET_PHI[9][100]/F");
    Tree_->Branch("RECO_MMET_CHARGE", RECO_MMET_CHARGE, "RECO_MMET_CHARGE[9][100]/F");

    //EEMT
    Tree_->Branch("RECO_EEMT_MASS", RECO_EEMT_MASS, "RECO_EEMT_MASS[9][100]/F");
    Tree_->Branch("RECO_EEMT_PT", RECO_EEMT_PT, "RECO_EEMT_PT[9][100]/F");
    Tree_->Branch("RECO_EEMT_ETA", RECO_EEMT_ETA, "RECO_EEMT_ETA[9][100]/F");
    Tree_->Branch("RECO_EEMT_PHI", RECO_EEMT_PHI, "RECO_EEMT_PHI[9][100]/F");
    Tree_->Branch("RECO_EEMT_CHARGE", RECO_EEMT_CHARGE, "RECO_EEMT_CHARGE[9][100]/F");

    //MMME
    Tree_->Branch("RECO_MMME_MASS", RECO_MMME_MASS, "RECO_MMME_MASS[9][100]/F");
    Tree_->Branch("RECO_MMME_PT", RECO_MMME_PT, "RECO_MMME_PT[9][100]/F");
    Tree_->Branch("RECO_MMME_ETA", RECO_MMME_ETA, "RECO_MMME_ETA[9][100]/F");
    Tree_->Branch("RECO_MMME_PHI", RECO_MMME_PHI, "RECO_MMME_PHI[9][100]/F");
    Tree_->Branch("RECO_MMME_CHARGE", RECO_MMME_CHARGE, "RECO_MMME_CHARGE[9][100]/F");

    //EEEM
    Tree_->Branch("RECO_EEEM_MASS", RECO_EEEM_MASS, "RECO_EEEM_MASS[9][100]/F");
    Tree_->Branch("RECO_EEEM_PT", RECO_EEEM_PT, "RECO_EEEM_PT[9][100]/F");
    Tree_->Branch("RECO_EEEM_ETA", RECO_EEEM_ETA, "RECO_EEEM_ETA[9][100]/F");
    Tree_->Branch("RECO_EEEM_PHI", RECO_EEEM_PHI, "RECO_EEEM_PHI[9][100]/F");
    Tree_->Branch("RECO_EEEM_CHARGE", RECO_EEEM_CHARGE, "RECO_EEEM_CHARGE[9][100]/F");

    //MMMM
    Tree_->Branch("RECO_MMMM_MASS", RECO_MMMM_MASS, "RECO_MMMM_MASS[9][100]/F");
    Tree_->Branch("RECO_MMMM_PT", RECO_MMMM_PT, "RECO_MMMM_PT[9][100]/F");
    Tree_->Branch("RECO_MMMM_ETA", RECO_MMMM_ETA, "RECO_MMMM_ETA[9][100]/F");
    Tree_->Branch("RECO_MMMM_PHI", RECO_MMMM_PHI, "RECO_MMMM_PHI[9][100]/F");
    Tree_->Branch("RECO_MMMM_CHARGE", RECO_MMMM_CHARGE, "RECO_MMMM_CHARGE[9][100]/F");

    //EEEE
    Tree_->Branch("RECO_EEEE_MASS", RECO_EEEE_MASS, "RECO_EEEE_MASS[9][100]/F");
    Tree_->Branch("RECO_EEEE_PT", RECO_EEEE_PT, "RECO_EEEE_PT[9][100]/F");
    Tree_->Branch("RECO_EEEE_ETA", RECO_EEEE_ETA, "RECO_EEEE_ETA[9][100]/F");
    Tree_->Branch("RECO_EEEE_PHI", RECO_EEEE_PHI, "RECO_EEEE_PHI[9][100]/F");
    Tree_->Branch("RECO_EEEE_CHARGE", RECO_EEEE_CHARGE, "RECO_EEEE_CHARGE[9][100]/F");

    //MMEE
    Tree_->Branch("RECO_MMEE_MASS", RECO_MMEE_MASS, "RECO_MMEE_MASS[9][100]/F");
    Tree_->Branch("RECO_MMEE_PT", RECO_MMEE_PT, "RECO_MMEE_PT[9][100]/F");
    Tree_->Branch("RECO_MMEE_ETA", RECO_MMEE_ETA, "RECO_MMEE_ETA[9][100]/F");
    Tree_->Branch("RECO_MMEE_PHI", RECO_MMEE_PHI, "RECO_MMEE_PHI[9][100]/F");
    Tree_->Branch("RECO_MMEE_CHARGE", RECO_MMEE_CHARGE, "RECO_MMEE_CHARGE[9][100]/F");

    //EEMM
    Tree_->Branch("RECO_EEMM_MASS", RECO_EEMM_MASS, "RECO_EEMM_MASS[9][100]/F");
    Tree_->Branch("RECO_EEMM_PT", RECO_EEMM_PT, "RECO_EEMM_PT[9][100]/F");
    Tree_->Branch("RECO_EEMM_ETA", RECO_EEMM_ETA, "RECO_EEMM_ETA[9][100]/F");
    Tree_->Branch("RECO_EEMM_PHI", RECO_EEMM_PHI, "RECO_EEMM_PHI[9][100]/F");
    Tree_->Branch("RECO_EEMM_CHARGE", RECO_EEMM_CHARGE, "RECO_EEMM_CHARGE[9][100]/F");

    //MT
    Tree_->Branch("RECO_MT_MASS", RECO_MT_MASS, "RECO_MT_MASS[5][100]/F");
    Tree_->Branch("RECO_MT_PT", RECO_MT_PT, "RECO_MT_PT[5][100]/F");
    Tree_->Branch("RECO_MT_ETA", RECO_MT_ETA, "RECO_MT_ETA[5][100]/F");
    Tree_->Branch("RECO_MT_PHI", RECO_MT_PHI, "RECO_MT_PHI[5][100]/F");
    Tree_->Branch("RECO_MT_CHARGE", RECO_MT_CHARGE, "RECO_MT_CHARGE[5][100]/F");

    //ET
    Tree_->Branch("RECO_ET_MASS", RECO_ET_MASS, "RECO_ET_MASS[5][100]/F");
    Tree_->Branch("RECO_ET_PT", RECO_ET_PT, "RECO_ET_PT[5][100]/F");
    Tree_->Branch("RECO_ET_ETA", RECO_ET_ETA, "RECO_ET_ETA[5][100]/F");
    Tree_->Branch("RECO_ET_PHI", RECO_ET_PHI, "RECO_ET_PHI[5][100]/F");
    Tree_->Branch("RECO_ET_CHARGE", RECO_ET_CHARGE, "RECO_ET_CHARGE[5][100]/F");

    //ME
    Tree_->Branch("RECO_ME_MASS", RECO_ME_MASS, "RECO_ME_MASS[5][100]/F");
    Tree_->Branch("RECO_ME_PT", RECO_ME_PT, "RECO_ME_PT[5][100]/F");
    Tree_->Branch("RECO_ME_ETA", RECO_ME_ETA, "RECO_ME_ETA[5][100]/F");
    Tree_->Branch("RECO_ME_PHI", RECO_ME_PHI, "RECO_ME_PHI[5][100]/F");
    Tree_->Branch("RECO_ME_CHARGE", RECO_ME_CHARGE, "RECO_ME_CHARGE[5][100]/F");

    /* Tree_->Branch("RECO_LLL0_MASS", RECO_LLL0_MASS, "RECO_LLL0_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLL1_MASS", RECO_LLL1_MASS, "RECO_LLL1_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLL2_MASS", RECO_LLL2_MASS, "RECO_LLL2_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLL3_MASS", RECO_LLL3_MASS, "RECO_LLL3_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLL0_PT", RECO_LLL0_PT, "RECO_LLL0_PT[4][50]/F");  */
/*     Tree_->Branch("RECO_LLL1_PT", RECO_LLL1_PT, "RECO_LLL1_PT[4][50]/F");  */
/*     Tree_->Branch("RECO_LLL2_PT", RECO_LLL2_PT, "RECO_LLL2_PT[4][50]/F");  */
/*     Tree_->Branch("RECO_LLL3_PT", RECO_LLL3_PT, "RECO_LLL3_PT[4][50]/F");  */

/*     Tree_->Branch("RECO_LLLl0_MASS", RECO_LLLl0_MASS, "RECO_LLLl0_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLLl1_MASS", RECO_LLLl1_MASS, "RECO_LLLl1_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLLl0_PT", RECO_LLLl0_PT, "RECO_LLLl0_PT[5][50]/F");  */
/*     Tree_->Branch("RECO_LLLl1_PT", RECO_LLLl1_PT, "RECO_LLLl1_PT[5][50]/F");  */

/*     Tree_->Branch("RECO_LLLL0ss_MASS", RECO_LLLL0ss_MASS, "RECO_LLLL0ss_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLLL0ss_PT", RECO_LLLL0ss_PT, "RECO_LLLL0ss_PT[5][50]/F");  */
/*     Tree_->Branch("RECO_LLLL1ss_MASS", RECO_LLLL1ss_MASS, "RECO_LLLL1ss_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLLL1ss_PT", RECO_LLLL1ss_PT, "RECO_LLLL1ss_PT[5][50]/F");  */
/*     Tree_->Branch("RECO_LLLL2ss_MASS", RECO_LLLL2ss_MASS, "RECO_LLLL2ss_MASS[50]/F");  */
/*     Tree_->Branch("RECO_LLLL2ss_PT", RECO_LLLL2ss_PT, "RECO_LLLL2ss_PT[5][50]/F");  */
       
/*     //Tree_->Branch("RECOcollNameLLLLssos_MASS",RECOcollNameLLLLssos_MASS,"RECOcollNameLLLLssos_MASS[50]/F"); */
/*     //Tree_->Branch("RECOcollNameLLLLssos_PT",RECOcollNameLLLLssos_PT,"RECOcollNameLLLLssos_PT[5][50]/F"); */

/*     Tree_->Branch("RECO_LLLL_MASS", RECO_LLLL_MASS, "RECO_LLLL_MASS[7][100]/F"); */
/*     Tree_->Branch("RECO_LLLL_PT", RECO_LLLL_PT, "RECO_LLLL_PT[7][100]/F"); */
/*     Tree_->Branch("RECO_LLLL_ETA", RECO_LLLL_ETA, "RECO_LLLL_ETA[7][100]/F"); */
/*     Tree_->Branch("RECO_LLLL_PHI", RECO_LLLL_PHI, "RECO_LLLL_PHI[7][100]/F"); */
 
    //HPSTau Block
    Tree_->Branch("HPSTAU_discByLER", HPSTAU_discByLER, "HPSTAU_discByLER[100]/D"); 
    Tree_->Branch("HPSTAU_discByMER", HPSTAU_discByMER, "HPSTAU_discByMER[100]/D"); 
    Tree_->Branch("HPSTAU_discByTER", HPSTAU_discByTER, "HPSTAU_discByTER[100]/D"); 
    Tree_->Branch("HPSTAU_discByLMR", HPSTAU_discByLMR, "HPSTAU_discByLMR[100]/D"); 
    Tree_->Branch("HPSTAU_discByTMR", HPSTAU_discByTMR, "HPSTAU_discByTMR[100]/D"); 
    Tree_->Branch("HPSTAU_discByDecayF", HPSTAU_discByDecayF, "HPSTAU_discByDecayF[100]/D"); 
    Tree_->Branch("HPSTAU_discLooseIso", HPSTAU_discLooseIso, "HPSTAU_discLooseIso[100]/D"); 
    Tree_->Branch("HPSTAU_discVLooseIso", HPSTAU_discVLooseIso, "HPSTAU_discVLooseIso[100]/D"); 
    Tree_->Branch("HPSTAU_discMediumIso", HPSTAU_discMediumIso, "HPSTAU_discMediumIso[100]/D"); 
    Tree_->Branch("HPSTAU_discTightIso", HPSTAU_discTightIso, "HPSTAU_discTightIso[100]/D"); 
    Tree_->Branch("HPSTAU_discLooseChargedIso", HPSTAU_discLooseChargedIso, "HPSTAU_discLooseChargedIso[100]/D"); 
    Tree_->Branch("HPSTAU_discVLooseChargedIso", HPSTAU_discVLooseChargedIso, "HPSTAU_discVLooseChargedIso[100]/D"); 
    Tree_->Branch("HPSTAU_discMediumChargedIso", HPSTAU_discMediumChargedIso, "HPSTAU_discMediumChargedIso[100]/D"); 
    Tree_->Branch("HPSTAU_discTightChargedIso", HPSTAU_discTightChargedIso, "HPSTAU_discTightChargedIso[100]/D"); 
    Tree_->Branch("HPSTAU_discLooseIsoDB", HPSTAU_discLooseIsoDB, "HPSTAU_discLooseIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_discVLooseIsoDB", HPSTAU_discVLooseIsoDB, "HPSTAU_discVLooseIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_discMediumIsoDB", HPSTAU_discMediumIsoDB, "HPSTAU_discMediumIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_discTightIsoDB", HPSTAU_discTightIsoDB, "HPSTAU_discTightIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_discLooseCombinedIsoDB", HPSTAU_discLooseCombinedIsoDB, "HPSTAU_discLooseCombinedIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_discVLooseCombinedIsoDB", HPSTAU_discVLooseCombinedIsoDB, "HPSTAU_discVLooseCombinedIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_discMediumCombinedIsoDB", HPSTAU_discMediumCombinedIsoDB, "HPSTAU_discMediumCombinedIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_discTightCombinedIsoDB", HPSTAU_discTightCombinedIsoDB, "HPSTAU_discTightCombinedIsoDB[100]/D"); 
    Tree_->Branch("HPSTAU_MASS", HPSTAU_MASS, "HPSTAU_MASS[100]/F"); 
    Tree_->Branch("HPSTAU_PT", HPSTAU_PT, "HPSTAU_PT[100]/F"); 
    Tree_->Branch("HPSTAU_ETA", HPSTAU_ETA, "HPSTAU_ETA[100]/F"); 
    Tree_->Branch("HPSTAU_PHI", HPSTAU_PHI, "HPSTAU_PHI[100]/F"); 
    Tree_->Branch("HPSTAU_CHARGE", HPSTAU_CHARGE, "HPSTAU_CHARGE[100]/F");   

    //ElectronDisc
    Tree_->Branch("EleId95relIso", EleId95relIso, "EleId95relIso[100]/D");
    Tree_->Branch("EleId90relIso", EleId90relIso, "EleId90relIso[100]/D");
    Tree_->Branch("EleId85relIso", EleId85relIso, "EleId85relIso[100]/D");
    Tree_->Branch("EleId80relIso", EleId80relIso, "EleId80relIso[100]/D");


    // Electron block
    Tree_->Branch("RECOELE_E", RECOELE_E, "RECOELE_E[100]/F"); 
    Tree_->Branch("RECOELE_PT",RECOELE_PT,"RECOELE_PT[100]/F");
    Tree_->Branch("RECOELE_P",RECOELE_P,"RECOELE_P[100]/F");
    Tree_->Branch("RECOELE_ETA",RECOELE_ETA,"RECOELE_ETA[100]/F"); 
    Tree_->Branch("RECOELE_THETA",RECOELE_THETA,"RECOELE_THETA[100]/F"); 
    Tree_->Branch("RECOELE_PHI",RECOELE_PHI,"RECOELE_PHI[100]/F"); 
    Tree_->Branch("RECOELE_MASS",RECOELE_MASS,"RECOELE_MASS[100]/F"); 
    Tree_->Branch("RECOELE_TMASS",RECOELE_TMASS,"RECOELE_TMASS[100]/F"); 
    Tree_->Branch("RECOELE_QUALITY",RECOELE_QUALITY,"RECOELE_QUALITY[100]/F"); 
    Tree_->Branch("RECOELE_CHARGE",RECOELE_CHARGE,"RECOELE_CHARGE[100]/F");   
    // Core attributes
    Tree_->Branch("RECOELE_isEcalDriven",   RECOELE_isEcalDriven,   "RECOELE_isEcalDriven[100]/b");   
    Tree_->Branch("RECOELE_isTrackerDriven",RECOELE_isTrackerDriven,"RECOELE_isTrackerDriven[100]/b");   
    Tree_->Branch("RECOELE_gsftrack_chi2",  RECOELE_gsftrack_chi2, "RECOELE_gsftrack_chi2[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dxyB",  RECOELE_gsftrack_dxyB,"RECOELE_gsftrack_dxyB[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dxy",   RECOELE_gsftrack_dxy,"RECOELE_gsftrack_dxy[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dxyError", RECOELE_gsftrack_dxyError,"RECOELE_gsftrack_dxyError[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dzB",   RECOELE_gsftrack_dzB,"RECOELE_gsftrack_dzB[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dz",    RECOELE_gsftrack_dz,"RECOELE_gsftrack_dz[100]/F");
    Tree_->Branch("RECOELE_gsftrack_dzError", RECOELE_gsftrack_dzError,"RECOELE_gsftrack_dzError[100]/F");
    Tree_->Branch("RECOELE_gsftrack_losthits", RECOELE_gsftrack_losthits,"RECOELE_gsftrack_losthits[100]/I");
    Tree_->Branch("RECOELE_gsftrack_validhits",RECOELE_gsftrack_validhits,"RECOELE_gsftrack_validhits[100]/I");
    Tree_->Branch("RECOELE_gsftrack_expected_inner_hits",RECOELE_gsftrack_expected_inner_hits,"RECOELE_gsftrack_expected_inner_hits[100]/I") ; 
    Tree_->Branch("RECOELE_scl_E",  RECOELE_scl_E,"RECOELE_scl_E[100]/F");
    Tree_->Branch("RECOELE_scl_Et", RECOELE_scl_Et,"RECOELE_scl_Et[100]/F");
    Tree_->Branch("RECOELE_scl_Eta",RECOELE_scl_Eta,"RECOELE_scl_Eta[100]/F");
    Tree_->Branch("RECOELE_scl_Phi",RECOELE_scl_Phi,"RECOELE_scl_Phi[100]/F");
    // Track-Cluster Matching
    Tree_->Branch("RECOELE_ep",          RECOELE_ep,"RECOELE_ep[100]/F");
    Tree_->Branch("RECOELE_eSeedp",      RECOELE_eSeedp,"RECOELE_eSeedp[100]/F");
    Tree_->Branch("RECOELE_eSeedpout",   RECOELE_eSeedpout,"RECOELE_eSeedpout[100]/F");
    Tree_->Branch("RECOELE_eElepout",    RECOELE_eElepout,"RECOELE_eElepout[100]/F");
    Tree_->Branch("RECOELE_deltaEtaIn",  RECOELE_deltaEtaIn,"RECOELE_deltaEtaIn[100]/F");
    Tree_->Branch("RECOELE_deltaEtaSeed",RECOELE_deltaEtaSeed,"RECOELE_deltaEtaSeed[100]/F");
    Tree_->Branch("RECOELE_deltaEtaEle", RECOELE_deltaEtaEle,"RECOELE_deltaEtaEle[100]/F");
    Tree_->Branch("RECOELE_deltaPhiIn",  RECOELE_deltaPhiIn,"RECOELE_deltaPhiIn[100]/F");
    Tree_->Branch("RECOELE_deltaPhiSeed",RECOELE_deltaPhiSeed,"RECOELE_deltaPhiSeed[100]/F");
    Tree_->Branch("RECOELE_deltaPhiEle", RECOELE_deltaPhiEle,"RECOELE_deltaPhiEle[100]/F");
    // Fiducial flags
    Tree_->Branch("RECOELE_isbarrel",    RECOELE_isbarrel,   "RECOELE_isbarrel[100]/I");   
    Tree_->Branch("RECOELE_isendcap",    RECOELE_isendcap,   "RECOELE_isendcap[100]/I");   
    Tree_->Branch("RECOELE_isEBetaGap",  RECOELE_isEBetaGap, "RECOELE_isEBetaGap[100]/I");   
    Tree_->Branch("RECOELE_isEBphiGap",  RECOELE_isEBphiGap, "RECOELE_isEBphiGap[100]/I");   
    Tree_->Branch("RECOELE_isEEdeeGap",  RECOELE_isEEdeeGap, "RECOELE_isEEdeeGap[100]/I");   
    Tree_->Branch("RECOELE_isEEringGap", RECOELE_isEEringGap,"RECOELE_isEEringGap[100]/I");   
    // Shower shape
    Tree_->Branch("RECOELE_sigmaIetaIeta", RECOELE_sigmaIetaIeta, "RECOELE_sigmaIetaIeta[100]/F");   
    Tree_->Branch("RECOELE_sigmaEtaEta",   RECOELE_sigmaEtaEta,   "RECOELE_sigmaEtaEta[100]/F");   
    Tree_->Branch("RECOELE_e15",           RECOELE_e15,           "RECOELE_e15[100]/F");   
    Tree_->Branch("RECOELE_e25max",        RECOELE_e25max,        "RECOELE_e25max[100]/F");   
    Tree_->Branch("RECOELE_e55",           RECOELE_e55,           "RECOELE_e55[100]/F");   
    Tree_->Branch("RECOELE_he",            RECOELE_he,            "RECOELE_he[100]/F");   
    // Particle flow
    Tree_->Branch("RECOELE_mva", RECOELE_mva,"RECOELE_mva[100]/F");   
    // Brem & Classifaction
    Tree_->Branch("RECOELE_fbrem",  RECOELE_fbrem,  "RECOELE_fbrem[100]/F");   
    Tree_->Branch("RECOELE_nbrems", RECOELE_nbrems, "RECOELE_nbrems[100]/I");   
    //  golden/bigbrem/(narrow)/showering/crack
    Tree_->Branch("RECOELE_Class",  RECOELE_Class,  "RECOELE_Class[100]/I");  
    //fBrem addition
    Tree_->Branch("RECOELE_fbrem_mode", &RECOELE_fbrem_mode,"RECOELE_fbrem_mode[100]/D");
    Tree_->Branch("RECOELE_fbrem_mean", &RECOELE_fbrem_mean,"RECOELE_fbrem_mean[100]/D");
    
    // Isolation
    Tree_->Branch("RECOELE_TRACKISO",RECOELE_TRACKISO,"RECOELE_TRACKISO[100]/F");  
    Tree_->Branch("RECOELE_HCALISO",RECOELE_HCALISO,"RECOELE_HCALISO[100]/F");  
    Tree_->Branch("RECOELE_ECALISO",RECOELE_ECALISO,"RECOELE_ECALISO[100]/F"); 
    Tree_->Branch("RECOELE_X",RECOELE_X,"RECOELE_X[100]/F"); 
    Tree_->Branch("RECOELE_EGMTRACKISO",RECOELE_EGMTRACKISO,"RECOELE_EGMTRACKISO[100]/F");  
    Tree_->Branch("RECOELE_EGMHCALISO",RECOELE_EGMHCALISO,"RECOELE_EGMHCALISO[100]/F");  
    Tree_->Branch("RECOELE_EGMECALISO",RECOELE_EGMECALISO,"RECOELE_EGMECALISO[100]/F"); 
    Tree_->Branch("RECOELE_EGMX",RECOELE_EGMX,"RECOELE_EGMX[100]/F"); 

    // Vertexing DA and KF
    Tree_->Branch("RECOELE_SIP",RECOELE_SIP,"RECOELE_SIP[100]/F"); 
    Tree_->Branch("RECOELE_IP",RECOELE_IP,"RECOELE_IP[100]/F"); 
    Tree_->Branch("RECOELE_IPERROR",RECOELE_IPERROR,"RECOELE_IPERROR[100]/F"); 
    Tree_->Branch("RECOELE_SIP_KF",RECOELE_SIP_KF,"RECOELE_SIP_KF[100]/F"); 
    Tree_->Branch("RECOELE_IP_KF",RECOELE_IP_KF,"RECOELE_IP_KF[100]/F"); 
    Tree_->Branch("RECOELE_IPERROR_KF",RECOELE_IPERROR_KF,"RECOELE_IPERROR_KF[100]/F"); 

     // GD vertex
    Tree_->Branch("RECOELE_SIP_GD",RECOELE_SIP_GD,"RECOELE_SIP_GD[100]/F"); //2e2mu
    Tree_->Branch("RECOELE_SIP_GDEEEE",RECOELE_SIP_GDEEEE,"RECOELE_SIP_GDEEEE[100]/F");  //4e
    // Std vertex
    Tree_->Branch("RECOELE_SIP_Std",RECOELE_SIP_Std,"RECOELE_SIP_Std[100]/F"); //2e2mu
    Tree_->Branch("RECOELE_SIP_StdEEEE",RECOELE_SIP_StdEEEE,"RECOELE_SIP_StdEEEE[100]/F");  //4e
    // Kin vertex
    Tree_->Branch("RECOELE_SIP_Kin",RECOELE_SIP_Kin,"RECOELE_SIP_Kin[100]/F"); //2e2mu
    Tree_->Branch("RECOELE_SIP_KinEEEE",RECOELE_SIP_KinEEEE,"RECOELE_SIP_KinEEEE[100]/F");  //4e


    Tree_->Branch("RECOELE_STIP",RECOELE_STIP,"RECOELE_STIP[100]/F"); 
    Tree_->Branch("RECOELE_SLIP",RECOELE_SLIP,"RECOELE_SLIP[100]/F"); 
    Tree_->Branch("RECOELE_TIP",RECOELE_TIP,"RECOELE_TIP[100]/F"); 
    Tree_->Branch("RECOELE_LIP",RECOELE_LIP,"RECOELE_LIP[100]/F"); 
    Tree_->Branch("RECOELE_TIPERROR",RECOELE_TIPERROR,"RECOELE_TIPERROR[100]/F"); 
    Tree_->Branch("RECOELE_LIPERROR",RECOELE_LIPERROR,"RECOELE_LIPERROR[100]/F"); 
    Tree_->Branch("RECOELE_eleID",RECOELE_eleID,"RECOELE_eleID[100][2]/b"); 

    Tree_->Branch("RECOELE_sclRawE",ele_sclRawE,"RECOELE_sclRawE[100]/D") ;
    Tree_->Branch("RECOELE_sclX",ele_sclX,"RECOELE_sclX[100]/D") ; 
    Tree_->Branch("RECOELE_sclY",ele_sclY,"RECOELE_sclY[100]/D") ; 
    Tree_->Branch("RECOELE_sclZ",ele_sclZ,"RECOELE_sclZ[100]/D") ;
    Tree_->Branch("RECOELE_seedSubdet1",ele_seedSubdet1,"RECOELE_seedSubdet1[100]/D") ;
    Tree_->Branch("RECOELE_seedDphi1",ele_seedDphi1,"RECOELE_seedDphi[100]/D") ; 
    Tree_->Branch("RECOELE_seedDrz1",ele_seedDrz1,"RECOELE_seedDrz1[100]/D") ;
    Tree_->Branch("RECOELE_seedSubdet2",ele_seedSubdet2,"RECOELE_seedSubdet2[100]/D") ;
    Tree_->Branch("RECOELE_seedDphi2",ele_seedDphi2,"RECOELE_seedDphi2[100]/D") ; 
    Tree_->Branch("RECOELE_seedDrz2",ele_seedDrz2,"RECOELE_seedDrz2[100]/D") ;
    Tree_->Branch("RECOELE_severityLevelSeed",ele_severityLevelSeed,"RECOELE_severityLevelSeed[100]/I") ; 
    Tree_->Branch("RECOELE_severityLevelClusters",ele_severityLevelClusters,"RECOELE_severityLevelClusters[100]/I") ; 
    Tree_->Branch("RECOELE_outOfTimeSeed",ele_outOfTimeSeed,"RECOELE_outOfTimeSeed[100]/I") ; 
    Tree_->Branch("RECOELE_outOfTimeClusters",ele_outOfTimeClusters,"RECOELE_outOfTimeClusters[100]/I") ;;
    Tree_->Branch("RECOELE_e2overe9",ele_e2overe9,"RECOELE_e2overe9[100]/D") ; 
    Tree_->Branch("RECOELE_eidVeryLoose",ele_eidVeryLoose,"RECOELE_eidVeryLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidLoose",ele_eidLoose,"RECOELE_eidLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidMedium",ele_eidMedium,"RECOELE_eidMedium[100]/D") ; 
    Tree_->Branch("RECOELE_eidTight",ele_eidTight,"RECOELE_eidTight[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZVeryLoose",ele_eidHZZVeryLoose,"RECOELE_eidHZZVeryLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZLoose",ele_eidHZZLoose,"RECOELE_eidHZZLoose[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZMedium",ele_eidHZZMedium,"RECOELE_eidHZZMedium[100]/D") ; 
    Tree_->Branch("RECOELE_eidHZZTight",ele_eidHZZTight,"RECOELE_eidHZZTight[100]/D") ; 

/*     Tree_->Branch("RECOELEBEST_2e2mu_MATCHED",RECOELEBEST_2e2mu_MATCHED,"RECOELEBEST_2e2mu_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELEBEST_4e_MATCHED",RECOELEBEST_4e_MATCHED,"RECOELEBEST_4e_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_ZEE_MATCHED",RECOELE_ZEE_MATCHED,"RECOELE_ZEE_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_ZssEE_MATCHED",RECOELE_ZssEE_MATCHED,"RECOELE_ZssEE_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_ZEM_MATCHED",RECOELE_ZEM_MATCHED,"RECOELE_ZEM_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_EEEE_MATCHED",RECOELE_EEEE_MATCHED,"RECOELE_EEEE_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_EEMM_MATCHED",RECOELE_EEMM_MATCHED,"RECOELE_EEMM_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLL0_MATCHED",RECOELE_LLL0_MATCHED,"RECOELE_LLL0_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLL1_MATCHED",RECOELE_LLL1_MATCHED,"RECOELE_LLL1_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLL2_MATCHED",RECOELE_LLL2_MATCHED,"RECOELE_LLL2_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLL3_MATCHED",RECOELE_LLL3_MATCHED,"RECOELE_LLL3_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLLLss0_MATCHED",RECOELE_LLLLss0_MATCHED,"RECOELE_LLLLss0_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLLLss1_MATCHED",RECOELE_LLLLss1_MATCHED,"RECOELE_LLLLss1_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLLLss2_MATCHED",RECOELE_LLLLss2_MATCHED,"RECOELE_LLLLss2_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLLl0_MATCHED",RECOELE_LLLl0_MATCHED,"RECOELE_LLLl0_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLLl1_MATCHED",RECOELE_LLLl1_MATCHED,"RECOELE_LLLl1_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOELE_LLLL_MATCHED",RECOELE_LLLL_MATCHED,"RECOELE_LLLL_MATCHED[100]/i"); */
    
    // Muon block
    Tree_->Branch("RECOMU_isGlobalMu", RECOMU_isGlobalMu, "RECOMU_isGlobalMu[100]/b");
    Tree_->Branch("RECOMU_isStandAloneMu", RECOMU_isStandAloneMu, "RECOMU_isStandAloneMu[100]/b");
    Tree_->Branch("RECOMU_isTrackerMu", RECOMU_isTrackerMu, "RECOMU_isTrackerMu[100]/b");
    Tree_->Branch("RECOMU_isCaloMu", RECOMU_isCaloMu, "RECOMU_isCaloMu[100]/b");
    Tree_->Branch("RECOMU_E",RECOMU_E,"RECOMU_E[100]/F"); 
    Tree_->Branch("RECOMU_PT",RECOMU_PT,"RECOMU_PT[100]/F"); 
    Tree_->Branch("RECOMU_P",RECOMU_P,"RECOMU_P[100]/F"); 
    Tree_->Branch("RECOMU_ETA",RECOMU_ETA,"RECOMU_ETA[100]/F"); 
    Tree_->Branch("RECOMU_THETA",RECOMU_THETA,"RECOMU_THETA[100]/F"); 
    Tree_->Branch("RECOMU_PHI",RECOMU_PHI,"RECOMU_PHI[100]/F"); 
    Tree_->Branch("RECOMU_MASS",RECOMU_MASS,"RECOMU_MASS[100]/F"); 
    Tree_->Branch("RECOMU_TMASS",RECOMU_TMASS,"RECOMU_TMASS[100]/F"); 
    Tree_->Branch("RECOMU_QUALITY",RECOMU_QUALITY,"RECOMU_QUALITY[100]/F"); 
    Tree_->Branch("RECOMU_CHARGE",RECOMU_CHARGE,"RECOMU_CHARGE[100]/F"); 
    Tree_->Branch("RECOMU_TRACKISO",RECOMU_TRACKISO,"RECOMU_TRACKISO[100]/F");  
    Tree_->Branch("RECOMU_HCALISO",RECOMU_HCALISO,"RECOMU_HCALISO[100]/F");  
    Tree_->Branch("RECOMU_ECALISO",RECOMU_ECALISO,"RECOMU_ECALISO[100]/F"); 
    Tree_->Branch("RECOMU_X",RECOMU_X,"RECOMU_X[100]/F"); 

    // Particle Flow Isolation
    //Tree_->Branch("RECOMU_ALLDEPOSITS",RECOMU_ALLDEPOSITS,"RECOMU_ALLDEPOSITS[100]/F");
    Tree_->Branch("RECOMU_ChargedDEPOSITS",RECOMU_ChargedDEPOSITS,"RECOMU_ChargedDEPOSITS[100]/F");
    Tree_->Branch("RECOMU_NeutralDEPOSITS",RECOMU_NeutralDEPOSITS,"RECOMU_NeutralDEPOSITS[100]/F");
    Tree_->Branch("RECOMU_GammaDEPOSITS",RECOMU_GammaDEPOSITS,"RECOMU_GammaDEPOSITS[100]/F");
    Tree_->Branch("RECOMU_PUDEPOSITS",RECOMU_PUDEPOSITS,"RECOMU_PUDEPOSITS[100]/F");
    //Tree_->Branch("RECOMU_PULowDEPOSITS",RECOMU_PULowDEPOSITS,"RECOMU_PULowDEPOSITS[100]/F");
    Tree_->Branch("RECOMU_PFX_DB",RECOMU_PFX_DB,"RECOMU_PFX_DB[100]/F");

    //Tree_->Branch("RECOELE_ALLDEPOSITS",RECOELE_ALLDEPOSITS,"RECOELE_ALLDEPOSITS[100]/F");
    Tree_->Branch("RECOELE_ChargedDEPOSITS",RECOELE_ChargedDEPOSITS,"RECOELE_ChargedDEPOSITS[100]/F");
    Tree_->Branch("RECOELE_NeutralDEPOSITS",RECOELE_NeutralDEPOSITS,"RECOELE_NeutralDEPOSITS[100]/F");
    Tree_->Branch("RECOELE_GammaDEPOSITS",RECOELE_GammaDEPOSITS,"RECOELE_GammaDEPOSITS[100]/F");
    Tree_->Branch("RECOELE_PUDEPOSITS",RECOELE_PUDEPOSITS,"RECOELE_PUDEPOSITS[100]/F");
    //Tree_->Branch("RECOELE_PULowDEPOSITS",RECOELE_PULowDEPOSITS,"RECOELE_PULowDEPOSITS[100]/F");
    Tree_->Branch("RECOELE_PFX_DB",RECOELE_PFX_DB,"RECOELE_PFX_DB[100]/F");
		  
    
    // vertexing DA and KF
    Tree_->Branch("RECOMU_SIP",RECOMU_SIP,"RECOMU_SIP[100]/F"); 
    Tree_->Branch("RECOMU_IP",RECOMU_IP,"RECOMU_IP[100]/F"); 
    Tree_->Branch("RECOMU_IPERROR",RECOMU_IPERROR,"RECOMU_IPERROR[100]/F"); 
    Tree_->Branch("RECOMU_SIP_KF",RECOMU_SIP_KF,"RECOMU_SIP_KF[100]/F"); 
    Tree_->Branch("RECOMU_IP_KF",RECOMU_IP_KF,"RECOMU_IP_KF[100]/F"); 
    Tree_->Branch("RECOMU_IPERROR_KF",RECOMU_IPERROR_KF,"RECOMU_IPERROR_KF[100]/F"); 

    // GD vertex
    Tree_->Branch("RECOMU_SIP_GD",RECOMU_SIP_GD,"RECOMU_SIP_GD[100]/F"); //2e2mu
    Tree_->Branch("RECOMU_SIP_GDMMMM",RECOMU_SIP_GDMMMM,"RECOMU_SIP_GDMMMM[100]/F");  //4mu
    // Std vertex
    Tree_->Branch("RECOMU_SIP_Std",RECOMU_SIP_Std,"RECOMU_SIP_Std[100]/F"); //2e2mu
    Tree_->Branch("RECOMU_SIP_StdMMMM",RECOMU_SIP_StdMMMM,"RECOMU_SIP_StdMMMM[100]/F");  //4mu
    // Kin vertex
    Tree_->Branch("RECOMU_SIP_Kin",RECOMU_SIP_Kin,"RECOMU_SIP_Kin[100]/F"); //2e2mu
    Tree_->Branch("RECOMU_SIP_KinMMMM",RECOMU_SIP_KinMMMM,"RECOMU_SIP_KinMMMM[100]/F");  //4mu



    Tree_->Branch("RECOMU_STIP",RECOMU_STIP,"RECOMU_STIP[100]/F"); 
    Tree_->Branch("RECOMU_SLIP",RECOMU_SLIP,"RECOMU_SLIP[100]/F"); 
    Tree_->Branch("RECOMU_TIP",RECOMU_TIP,"RECOMU_TIP[100]/F"); 
    Tree_->Branch("RECOMU_LIP",RECOMU_LIP,"RECOMU_LIP[100]/F"); 
    Tree_->Branch("RECOMU_TIPERROR",RECOMU_TIPERROR,"RECOMU_TIPERROR[100]/F"); 
    Tree_->Branch("RECOMU_LIPERROR",RECOMU_LIPERROR,"RECOMU_LIPERROR[100]/F"); 
    
/*     Tree_->Branch("RECOMUBEST_2e2mu_MATCHED",RECOMUBEST_2e2mu_MATCHED,"RECOMUBEST_2e2mu_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMUBEST_4mu_MATCHED",RECOMUBEST_4mu_MATCHED,"RECOMUBEST_4mu_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_ZMM_MATCHED",RECOMU_ZMM_MATCHED,"RECOMU_ZMM_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_ZssMM_MATCHED",RECOMU_ZssMM_MATCHED,"RECOMU_ZssMM_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_ZEM_MATCHED",RECOMU_ZEM_MATCHED,"RECOMU_ZEM_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_MMMM_MATCHED",RECOMU_MMMM_MATCHED,"RECOMU_MMMM_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_EEMM_MATCHED",RECOMU_EEMM_MATCHED,"RECOMU_EEMM_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLL0_MATCHED",RECOMU_LLL0_MATCHED,"RECOMU_LLL0_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLL1_MATCHED",RECOMU_LLL1_MATCHED,"RECOMU_LLL1_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLL2_MATCHED",RECOMU_LLL2_MATCHED,"RECOMU_LLL2_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLL3_MATCHED",RECOMU_LLL3_MATCHED,"RECOMU_LLL3_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLLLss0_MATCHED",RECOMU_LLLLss0_MATCHED,"RECOMU_LLLLss0_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLLLss1_MATCHED",RECOMU_LLLLss1_MATCHED,"RECOMU_LLLLss1_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLLLss2_MATCHED",RECOMU_LLLLss2_MATCHED,"RECOMU_LLLLss2_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLLl0_MATCHED",RECOMU_LLLl0_MATCHED,"RECOMU_LLLl0_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLLl1_MATCHED",RECOMU_LLLl1_MATCHED,"RECOMU_LLLl1_MATCHED[100]/i"); */
/*     Tree_->Branch("RECOMU_LLLL_MATCHED",RECOMU_LLLL_MATCHED,"RECOMU_LLLL_MATCHED[100]/i"); */
 

    Tree_->Branch("RECOMU_caloCompatibility",RECOMU_caloCompatibility,"RECOMU_caloCompatibility[100]/F");
    Tree_->Branch("RECOMU_segmentCompatibility",RECOMU_segmentCompatibility,"RECOMU_segmentCompatibility[100]/F"); 
    Tree_->Branch("RECOMU_numberOfMatches",RECOMU_numberOfMatches,"RECOMU_numberOfMatches[100]/i");
    Tree_->Branch("RECOMU_glbmuPromptTight",RECOMU_glbmuPromptTight,"RECOMU_glbmuPromptTight[100]/b");
 
    // track variables from muons:
    Tree_->Branch( "RECOMU_trkmuArbitration", RECOMU_trkmuArbitration, "RECOMU_trkmuArbitration[100]/b");
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityLoose", RECOMU_trkmu2DCompatibilityLoose, "RECOMU_trkmu2DCompatibilityLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmu2DCompatibilityTight", RECOMU_trkmu2DCompatibilityTight, "RECOMU_trkmu2DCompatibilityTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationLoose", RECOMU_trkmuOneStationLoose, "RECOMU_trkmuOneStationLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationTight", RECOMU_trkmuOneStationTight, "RECOMU_trkmuOneStationTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationLoose", RECOMU_trkmuLastStationLoose, "RECOMU_trkmuLastStationLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationTight", RECOMU_trkmuLastStationTight, "RECOMU_trkmuLastStationTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationAngLoose", RECOMU_trkmuOneStationAngLoose, "RECOMU_trkmuOneStationAngLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuOneStationAngTight", RECOMU_trkmuOneStationAngTight, "RECOMU_trkmuOneStationAngTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationAngLoose", RECOMU_trkmuLastStationAngLoose, "RECOMU_trkmuLastStationAngLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationAngTight", RECOMU_trkmuLastStationAngTight, "RECOMU_trkmuLastStationAngTight[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtLoose",RECOMU_trkmuLastStationOptimizedLowPtLoose , "RECOMU_trkmuLastStationOptimizedLowPtLoose[100]/b");
    Tree_->Branch( "RECOMU_trkmuLastStationOptimizedLowPtTight",RECOMU_trkmuLastStationOptimizedLowPtTight , "RECOMU_trkmuLastStationOptimizedLowPtTight[100]/b");

    Tree_->Branch( "RECOMU_mutrkDxy", RECOMU_mutrkDxy, "RECOMU_mutrkDxy[100]/F");
    Tree_->Branch( "RECOMU_mutrkDxyError", RECOMU_mutrkDxyError, "RECOMU_mutrkDxyError[100]/F");
    Tree_->Branch( "RECOMU_mutrkDxyB", RECOMU_mutrkDxyB, "RECOMU_mutrkDxyB[100]/F");
    Tree_->Branch( "RECOMU_mutrkDz", RECOMU_mutrkDz, "RECOMU_mutrkDz[100]/F");
    Tree_->Branch( "RECOMU_mutrkDzError", RECOMU_mutrkDzError, "RECOMU_mutrkDzError[100]/F");
    Tree_->Branch( "RECOMU_mutrkDzB", RECOMU_mutrkDzB, "RECOMU_mutrkDzB[100]/F");
    Tree_->Branch( "RECOMU_mutrkChi2PerNdof", RECOMU_mutrkChi2PerNdof, "RECOMU_mutrkChi2PerNdof[100]/F");
    Tree_->Branch( "RECOMU_mutrkCharge", RECOMU_mutrkCharge, "RECOMU_mutrkCharge[100]/F");
    Tree_->Branch( "RECOMU_mutrkNHits", RECOMU_mutrkNHits, "RECOMU_mutrkNHits[100]/F");
    Tree_->Branch( "RECOMU_mutrkNStripHits", RECOMU_mutrkNStripHits, "RECOMU_mutrkNStripHits[100]/F");
    Tree_->Branch( "RECOMU_mutrkNPixHits", RECOMU_mutrkNPixHits, "RECOMU_mutrkNPixHits[100]/F");
    Tree_->Branch( "RECOMU_mutrkNMuonHits", RECOMU_mutrkNMuonHits, "RECOMU_mutrkNMuonHits[100]/F");
    
    
    // Isolation
    Tree_->Branch("isolation",  &isolation,  "isolation[4]/F");
    Tree_->Branch("isoele",  &isoele,  "isoele[2]/F");
    Tree_->Branch("isomu",  &isomu,  "isomu[2]/F");
    
    // vertexing
    Tree_->Branch("vertexing",  &vertexing,  "vertexing[4]/F");
    /* // Geom. Discri. */
/*     Tree_->Branch("ftsigma",        &ftsigma,        "ftsigma[4]/D"); */
/*     Tree_->Branch("gdX",            &gdX,            "gdX[4]/D"); */
/*     Tree_->Branch("gdY",            &gdY,            "gdY[4]/D"); */
/*     Tree_->Branch("gdZ",            &gdZ,            "gdZ[4]/D"); */
/*     Tree_->Branch("ftsigmalag",     &ftsigmalag,     "ftsigmalag[4]/D"); */
/*     Tree_->Branch("gdlagX",         &gdlagX,         "gdlagX[4]/D"); */
/*     Tree_->Branch("gdlagY",         &gdlagY,         "gdlagY[4]/D"); */
/*     Tree_->Branch("gdlagZ",         &gdlagZ,         "gdlagZ[4]/D"); */
/*     Tree_->Branch("gdlagProb",      &gdlagProb,      "gdlagProb[4]/D"); */
/*     Tree_->Branch("gdlagNdof",      &gdlagNdof,      "gdlagNdof[4]/D"); */
/*     Tree_->Branch("ftsigmaMMMM",    &ftsigmaMMMM,    "ftsigmaMMMM[4]/D"); */
/*     Tree_->Branch("gdXMMMM",        &gdXMMMM,        "gdXMMMM[4]/D"); */
/*     Tree_->Branch("gdYMMMM",        &gdYMMMM,        "gdYMMMM[4]/D"); */
/*     Tree_->Branch("gdZMMMM",        &gdZMMMM,        "gdZMMMM[4]/D"); */
/*     Tree_->Branch("ftsigmalagMMMM", &ftsigmalagMMMM, "ftsigmalagMMMM[4]/D"); */
/*     Tree_->Branch("gdlagXMMMM",     &gdlagXMMMM,     "gdlagXMMMM[4]/D"); */
/*     Tree_->Branch("gdlagYMMMM",     &gdlagYMMMM,     "gdlagYMMMM[4]/D"); */
/*     Tree_->Branch("gdlagZMMMM",     &gdlagZMMMM,     "gdlagZMMMM[4]/D"); */
/*     Tree_->Branch("gdlagProbMMMM",  &gdlagProbMMMM,  "gdlagProbMMMM[4]/D"); */
/*     Tree_->Branch("gdlagNdofMMMM",  &gdlagNdofMMMM,  "gdlagNdofMMMM[4]/D"); */
/*     Tree_->Branch("ftsigmaEEEE",    &ftsigmaEEEE,    "ftsigmaEEEE[4]/D"); */
/*     Tree_->Branch("gdXEEEE",        &gdXEEEE,        "gdXEEEE[4]/D"); */
/*     Tree_->Branch("gdYEEEE",        &gdYEEEE,        "gdYEEEE[4]/D"); */
/*     Tree_->Branch("gdZEEEE",        &gdZEEEE,        "gdZEEEE[4]/D"); */
/*     Tree_->Branch("ftsigmalagEEEE", &ftsigmalagEEEE, "ftsigmalagEEEE[4]/D"); */
/*     Tree_->Branch("gdlagXEEEE",     &gdlagXEEEE,     "gdlagXEEEE[4]/D"); */
/*     Tree_->Branch("gdlagYEEEE",     &gdlagYEEEE,     "gdlagYEEEE[4]/D"); */
/*     Tree_->Branch("gdlagZEEEE",     &gdlagZEEEE,     "gdlagZEEEE[4]/D"); */
/*     Tree_->Branch("gdlagProbEEEE",  &gdlagProbEEEE,  "gdlagProbEEEE[4]/D"); */
/*     Tree_->Branch("gdlagNdofEEEE",  &gdlagNdofEEEE,  "gdlagNdofEEEE[4]/D"); */
    
/*     // ConstraintFit 4l */
/*     Tree_->Branch("StdFitVertexX",        &StdFitVertexX,        "StdFitVertexX[4]/D"); */
/*     Tree_->Branch("StdFitVertexY",        &StdFitVertexY,        "StdFitVertexY[4]/D"); */
/*     Tree_->Branch("StdFitVertexZ",        &StdFitVertexZ,        "StdFitVertexZ[4]/D"); */
/*     Tree_->Branch("StdFitVertexChi2r",    &StdFitVertexChi2r,    "StdFitVertexChi2r[4]/D"); */
/*     Tree_->Branch("StdFitVertexProb",     &StdFitVertexProb,     "StdFitVertexProb[4]/D"); */
/*     Tree_->Branch("KinFitVertexX",        &KinFitVertexX,        "KinFitVertexX[4]/D"); */
/*     Tree_->Branch("KinFitVertexY",        &KinFitVertexY,        "KinFitVertexY[4]/D"); */
/*     Tree_->Branch("KinFitVertexZ",        &KinFitVertexZ,        "KinFitVertexZ[4]/D"); */
/*     Tree_->Branch("KinFitVertexChi2r",    &KinFitVertexChi2r,    "KinFitVertexChi2r[4]/D"); */
/*     Tree_->Branch("KinFitVertexProb",     &KinFitVertexProb,     "KinFitVertexProb[4]/D"); */

/*     Tree_->Branch("StdFitVertexXMMMM",        &StdFitVertexXMMMM,        "StdFitVertexXMMMM[4]/D"); */
/*     Tree_->Branch("StdFitVertexYMMMM",        &StdFitVertexYMMMM,        "StdFitVertexYMMMM[4]/D"); */
/*     Tree_->Branch("StdFitVertexZMMMM",        &StdFitVertexZMMMM,        "StdFitVertexZMMMM[4]/D"); */
/*     Tree_->Branch("StdFitVertexChi2rMMMM",    &StdFitVertexChi2rMMMM,    "StdFitVertexChi2rMMMM[4]/D"); */
/*     Tree_->Branch("StdFitVertexProbMMMM",     &StdFitVertexProbMMMM,     "StdFitVertexProbMMMM[4]/D"); */
/*     Tree_->Branch("KinFitVertexXMMMM",        &KinFitVertexXMMMM,        "KinFitVertexXMMMM[4]/D"); */
/*     Tree_->Branch("KinFitVertexYMMMM",        &KinFitVertexYMMMM,        "KinFitVertexYMMMM[4]/D"); */
/*     Tree_->Branch("KinFitVertexZMMMM",        &KinFitVertexZMMMM,        "KinFitVertexZMMMM[4]/D"); */
/*     Tree_->Branch("KinFitVertexChi2rMMMM",    &KinFitVertexChi2rMMMM,    "KinFitVertexChi2rMMMM[4]/D"); */
/*     Tree_->Branch("KinFitVertexProbMMMM",     &KinFitVertexProbMMMM,     "KinFitVertexProbMMMM[4]/D"); */

/*     Tree_->Branch("StdFitVertexXEEEE",        &StdFitVertexXEEEE,        "StdFitVertexXEEEE[4]/D"); */
/*     Tree_->Branch("StdFitVertexYEEEE",        &StdFitVertexYEEEE,        "StdFitVertexYEEEE[4]/D"); */
/*     Tree_->Branch("StdFitVertexZEEEE",        &StdFitVertexZEEEE,        "StdFitVertexZEEEE[4]/D"); */
/*     Tree_->Branch("StdFitVertexChi2rEEEE",    &StdFitVertexChi2rEEEE,    "StdFitVertexChi2rEEEE[4]/D"); */
/*     Tree_->Branch("StdFitVertexProbEEEE",     &StdFitVertexProbEEEE,     "StdFitVertexProbEEEE[4]/D"); */
/*     Tree_->Branch("KinFitVertexXEEEE",        &KinFitVertexXEEEE,        "KinFitVertexXEEEE[4]/D"); */
/*     Tree_->Branch("KinFitVertexYEEEE",        &KinFitVertexYEEEE,        "KinFitVertexYEEEE[4]/D"); */
/*     Tree_->Branch("KinFitVertexZEEEE",        &KinFitVertexZEEEE,        "KinFitVertexZEEEE[4]/D"); */
/*     Tree_->Branch("KinFitVertexChi2rEEEE",    &KinFitVertexChi2rEEEE,    "KinFitVertexChi2rEEEE[4]/D"); */
/*     Tree_->Branch("KinFitVertexProbEEEE",     &KinFitVertexProbEEEE,     "KinFitVertexProbEEEE[4]/D"); */

/*     // constrintFit 3l */
/*     Tree_->Branch("StdFitVertexChi2rMMM",    &StdFitVertexChi2rMMM,    "StdFitVertexChi2rMMM[20]/D"); */
/*     Tree_->Branch("StdFitVertexProbMMM",     &StdFitVertexProbMMM,     "StdFitVertexProbMMM[20]/D"); */
/*     Tree_->Branch("StdFitVertexChi2rMME",    &StdFitVertexChi2rMME,    "StdFitVertexChi2rMME[20]/D"); */
/*     Tree_->Branch("StdFitVertexProbMME",     &StdFitVertexProbMME,     "StdFitVertexProbMME[20]/D"); */
/*     Tree_->Branch("StdFitVertexChi2rEEE",    &StdFitVertexChi2rEEE,    "StdFitVertexChi2rEEE[20]/D"); */
/*     Tree_->Branch("StdFitVertexProbEEE",     &StdFitVertexProbEEE,     "StdFitVertexProbEEE[20]/D"); */
/*     Tree_->Branch("StdFitVertexChi2rMEE",    &StdFitVertexChi2rMEE,    "StdFitVertexChi2rMEE[20]/D"); */
/*     Tree_->Branch("StdFitVertexProbMEE",     &StdFitVertexProbMEE,     "StdFitVertexProbMEE[20]/D"); */


/*      // constrintFit Dileptons */
/*     Tree_->Branch("StdFitVertexChi2rDiLep",    &StdFitVertexChi2rDiLep,    "StdFitVertexChi2rDiLep[40]/D"); */
/*     Tree_->Branch("StdFitVertexProbDiLep",     &StdFitVertexProbDiLep,     "StdFitVertexProbDiLep[40]/D"); */

    // Conversions
    Tree_->Branch("ConvMapDist",              ConvMapDist,              "ConvMapDist[100]/F");
    Tree_->Branch("ConvMapDcot",              ConvMapDcot,              "ConvMapDcot[100]/F");

    //Global Event 
    Tree_->Branch( "RECO_NMU", &RECO_NMU, "RECO_NMU/I"); 
    Tree_->Branch( "RECO_NELE", &RECO_NELE, "RECO_NELE/I"); 
    
    // Tracks
    Tree_->Branch( "RECO_NTRACK", &RECO_NTRACK, "RECO_NTRACK/I");
    Tree_->Branch( "RECO_TRACK_PT", &RECO_TRACK_PT, "RECO_TRACK_PT[200]/F");
    Tree_->Branch( "RECO_TRACK_ETA", &RECO_TRACK_ETA, "RECO_TRACK_ETA[200]/F");
    Tree_->Branch( "RECO_TRACK_PHI", &RECO_TRACK_PHI, "RECO_TRACK_PHI[200]/F");
    Tree_->Branch( "RECO_TRACK_CHI2", &RECO_TRACK_CHI2, "RECO_TRACK_CHI2[200]/F");
    Tree_->Branch( "RECO_TRACK_CHI2RED", &RECO_TRACK_CHI2RED, "RECO_TRACK_CHI2RED[200]/F");
    Tree_->Branch( "RECO_TRACK_CHI2PROB", &RECO_TRACK_CHI2PROB, "RECO_TRACK_CHI2PROB[200]/F");
    Tree_->Branch( "RECO_TRACK_NHITS", &RECO_TRACK_NHITS, "RECO_TRACK_NHITS[200]/I");
    Tree_->Branch( "RECO_TRACK_DXY", &RECO_TRACK_DXY, "RECO_TRACK_DXY[200]/F");
    Tree_->Branch( "RECO_TRACK_DXYERR", &RECO_TRACK_DXYERR, "RECO_TRACK_DXYERR[200]/F");
    Tree_->Branch( "RECO_TRACK_DZ", &RECO_TRACK_DZ, "RECO_TRACK_DZ[200]/F");
    Tree_->Branch( "RECO_TRACK_DZERR", &RECO_TRACK_DZERR, "RECO_TRACK_DZERR[200]/F");
    
    Tree_->Branch( "RECO_NPHOT", &RECO_NPHOT, "RECO_NPHOT/I");
    
    //Beam Spot position
    Tree_->Branch("BeamSpot_X",&BeamSpot_X,"BeamSpot_X/D");
    Tree_->Branch("BeamSpot_Y",&BeamSpot_Y,"BeamSpot_Y/D");
    Tree_->Branch("BeamSpot_Z",&BeamSpot_Z,"BeamSpot_Z/D");
    // Vertices
    Tree_->Branch( "RECO_NVTX", &RECO_NVTX, "RECO_NVTX/I");
    Tree_->Branch( "RECO_VERTEX_x", &RECO_VERTEX_x, "RECO_VERTEX_x[10]/F");
    Tree_->Branch( "RECO_VERTEX_y", &RECO_VERTEX_y, "RECO_VERTEX_y[10]/F");
    Tree_->Branch( "RECO_VERTEX_z", &RECO_VERTEX_z, "RECO_VERTEX_z[10]/F");
    Tree_->Branch( "RECO_VERTEX_ndof", &RECO_VERTEX_ndof, "RECO_VERTEX_ndof[10]/F");
    Tree_->Branch( "RECO_VERTEX_chi2", &RECO_VERTEX_chi2, "RECO_VERTEX_chi2[10]/F");
    Tree_->Branch( "RECO_VERTEX_ntracks", &RECO_VERTEX_ntracks, "RECO_VERTEX_ntracks[10]/I");
    Tree_->Branch( "RECO_VERTEXPROB", &RECO_VERTEXPROB, "RECO_VERTEXPROB[10]/F");
    Tree_->Branch( "RECO_VERTEX_isValid", &RECO_VERTEX_isValid, "RECO_VERTEX_isValid[10]/b");
    Tree_->Branch( "RECO_VERTEX_TRACK_PT",&RECO_VERTEX_TRACK_PT,"RECO_VERTEX_TRACK_PT[10][30]/F");
    
    // PFJets
    Tree_->Branch( "RECO_PFJET_N",   &RECO_PFJET_N,   "RECO_PFJET_N/I");
    Tree_->Branch( "RECO_PFJET_CHARGE",  &RECO_PFJET_CHARGE,  "RECO_PFJET_CHARGE[100]/I");
    Tree_->Branch( "RECO_PFJET_ET",  &RECO_PFJET_ET,  "RECO_PFJET_ET[100]/F");
    Tree_->Branch( "RECO_PFJET_PT",  &RECO_PFJET_PT,  "RECO_PFJET_PT[100]/F");
    Tree_->Branch( "RECO_PFJET_ETA", &RECO_PFJET_ETA, "RECO_PFJET_ETA[100]/F");
    Tree_->Branch( "RECO_PFJET_PHI", &RECO_PFJET_PHI, "RECO_PFJET_PHI[100]/F");
    Tree_->Branch( "RHO", &RHO, "RHO/D");
    
   /*  //CaloMET */
/*     Tree_->Branch( "RECO_CALOMET",          &calomet,          "RECO_CALOMET/F"); */
/*     Tree_->Branch( "RECO_CALOMETHO",        &calometho,        "RECO_CALOMETHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHFHO",    &calometnohfho,    "RECO_CALOMETNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETNOHF",      &calometnohf,      "RECO_CALOMETNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTHO",     &calometoptho,     "RECO_CALOMETOPTHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHFHO", &calometoptnohfho, "RECO_CALOMETOPTNOHFHO/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPTNOHF",   &calometoptnohf,   "RECO_CALOMETOPTNOHF/F"); */
/*     Tree_->Branch( "RECO_CALOMETOPT",       &calometopt,       "RECO_CALOMETOPT/F"); */
    //Particle Flow MET
    Tree_->Branch( "RECO_PFMET", &pfmet, "RECO_PFMET/F");
    Tree_->Branch( "RECO_PFMET_X", &pfmet_x, "RECO_PFMET_X/F");
    Tree_->Branch( "RECO_PFMET_Y", &pfmet_y, "RECO_PFMET_Y/F");
    Tree_->Branch( "RECO_PFMET_PHI", &pfmet_phi, "RECO_PFMET_PHI/F");
    Tree_->Branch( "RECO_PFMET_THETA", &pfmet_theta, "RECO_PFMET_THETA/F");
   /*  //Track Corrected MET */
/*     Tree_->Branch( "RECO_TCMET", &tcmet, "RECO_TCMET/F"); */
/*     //HTMET */
/*     Tree_->Branch( "RECO_HTMETIC5", &htmetic5, "RECO_HTMETIC5/F"); */
/*     Tree_->Branch( "RECO_HTMETKT4", &htmetkt4, "RECO_HTMETKT4/F"); */
/*     Tree_->Branch( "RECO_HTMETKT6", &htmetkt6, "RECO_HTMETKT6/F"); */
/*     Tree_->Branch( "RECO_HTMETSC5", &htmetsc5, "RECO_HTMETSC5/F"); */
/*     Tree_->Branch( "RECO_HTMETSC7", &htmetsc7, "RECO_HTMETSC7/F"); */
/*     //JES Correction MET */
/*     Tree_->Branch( "RECO_metJESCorIC5", &jescormetic5, "RECO_metJESCorIC5/F"); */
/*     Tree_->Branch( "RECO_metJESCorKT4", &jescormetkt4, "RECO_metJESCorKT4/F"); */
/*     Tree_->Branch( "RECO_metJESCorKT6", &jescormetkt6, "RECO_metJESCorKT6/F"); */
/*     Tree_->Branch( "RECO_metJESCorSC5", &jescormetsc5, "RECO_metJESCorsc4/F"); */
/*     Tree_->Branch( "RECO_metJESCorSC7", &jescormetsc7, "RECO_metJESCorsc7/F"); */
/*     //Type I correction MET */
/*     Tree_->Branch( "RECO_CORMETMUONS",  &cormetmuons,  "RECO_CORMETMUONS/F"); */
   

/*     // Btagging jets and discriminators */
/*     Tree_->Branch("tCHighEff_BTagJet_PT",tCHighEff_BTagJet_PT,"tCHighEff_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("tCHighPur_BTagJet_PT", tCHighPur_BTagJet_PT,"tCHighPur_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("jPHighEff_BTagJet_PT",jPHighEff_BTagJet_PT,"jPHighEff_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("jBP_BTagJet_PT",jBP_BTagJet_PT,"jBP_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("sSVHighEff_BTagJet_PT",sSVHighEff_BTagJet_PT,"sSVHighEff_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("sSVHighPur_BTagJet_PT",sSVHighPur_BTagJet_PT,"sSVHighPur_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("cSV_BTagJet_PT",cSV_BTagJet_PT,"cSV_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("cSVMVA_BTagJet_PT",cSVMVA_BTagJet_PT,"cSVMVA_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("sEByIP3d_BTagJet_PT",sEByIP3d_BTagJet_PT,"sEByIP3d_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("sEByPt_BTagJet_PT",sEByPt_BTagJet_PT,"sEByPt_BTagJet_PT[50],/F"); */
/*     Tree_->Branch("sM_BTagJet_PT",sM_BTagJet_PT,"sM_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("sMByIP3d_BTagJet_PT",sMByIP3d_BTagJet_PT,"sMByIP3d_BTagJet_PT[50]/F"); */
/*     Tree_->Branch("sMByPt_BTagJet_PT",sMByPt_BTagJet_PT,"sMByPt_BTagJet_PT[50]/F"); */

/*     Tree_->Branch("tCHighEff_BTagJet_ETA",tCHighEff_BTagJet_ETA,"tCHighEff_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("tCHighPur_BTagJet_ETA", tCHighPur_BTagJet_ETA,"tCHighPur_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("jPHighEff_BTagJet_ETA",jPHighEff_BTagJet_ETA,"jPHighEff_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("jBP_BTagJet_ETA",jBP_BTagJet_ETA,"jBP_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("sSVHighEff_BTagJet_ETA",sSVHighEff_BTagJet_ETA,"sSVHighEff_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("sSVHighPur_BTagJet_ETA",sSVHighPur_BTagJet_ETA,"sSVHighPur_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("cSV_BTagJet_ETA",cSV_BTagJet_ETA,"cSV_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("cSVMVA_BTagJet_ETA",cSVMVA_BTagJet_ETA,"cSVMVA_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("sEByIP3d_BTagJet_ETA",sEByIP3d_BTagJet_ETA,"sEByIP3d_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("sEByPt_BTagJet_ETA",sEByPt_BTagJet_ETA,"sEByPt_BTagJet_ETA[50],/F"); */
/*     Tree_->Branch("sM_BTagJet_ETA",sM_BTagJet_ETA,"sM_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("sMByIP3d_BTagJet_ETA",sMByIP3d_BTagJet_ETA,"sMByIP3d_BTagJet_ETA[50]/F"); */
/*     Tree_->Branch("sMByPt_BTagJet_ETA",sMByPt_BTagJet_ETA,"sMByPt_BTagJet_ETA[50]/F"); */

/*     Tree_->Branch("tCHighEff_BTagJet_PHI",tCHighEff_BTagJet_PHI,"tCHighEff_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("tCHighPur_BTagJet_PHI", tCHighPur_BTagJet_PHI,"tCHighPur_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("jPHighEff_BTagJet_PHI",jPHighEff_BTagJet_PHI,"jPHighEff_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("jBP_BTagJet_PHI",jBP_BTagJet_PHI,"jBP_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("sSVHighEff_BTagJet_PHI",sSVHighEff_BTagJet_PHI,"sSVHighEff_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("sSVHighPur_BTagJet_PHI",sSVHighPur_BTagJet_PHI,"sSVHighPur_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("cSV_BTagJet_PHI",cSV_BTagJet_PHI,"cSV_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("cSVMVA_BTagJet_PHI",cSVMVA_BTagJet_PHI,"cSVMVA_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("sEByIP3d_BTagJet_PHI",sEByIP3d_BTagJet_PHI,"sEByIP3d_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("sEByPt_BTagJet_PHI",sEByPt_BTagJet_PHI,"sEByPt_BTagJet_PHI[50],/F"); */
/*     Tree_->Branch("sM_BTagJet_PHI",sM_BTagJet_PHI,"sM_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("sMByIP3d_BTagJet_PHI",sMByIP3d_BTagJet_PHI,"sMByIP3d_BTagJet_PHI[50]/F"); */
/*     Tree_->Branch("sMByPt_BTagJet_PHI",sMByPt_BTagJet_PHI,"sMByPt_BTagJet_PHI[50]/F"); */
 
/*     Tree_->Branch("tCHighEff_BTagJet_DISCR",tCHighEff_BTagJet_DISCR,"tCHighEff_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("tCHighPur_BTagJet_DISCR", tCHighPur_BTagJet_DISCR,"tCHighPur_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("jPHighEff_BTagJet_DISCR",jPHighEff_BTagJet_DISCR,"jPHighEff_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("jBP_BTagJet_DISCR",jBP_BTagJet_DISCR,"jBP_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("sSVHighEff_BTagJet_DISCR",sSVHighEff_BTagJet_DISCR,"sSVHighEff_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("sSVHighPur_BTagJet_DISCR",sSVHighPur_BTagJet_DISCR,"sSVHighPur_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("cSV_BTagJet_DISCR",cSV_BTagJet_DISCR,"cSV_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("cSVMVA_BTagJet_DISCR",cSVMVA_BTagJet_DISCR,"cSVMVA_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("sEByIP3d_BTagJet_DISCR",sEByIP3d_BTagJet_DISCR,"sEByIP3d_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("sEByPt_BTagJet_DISCR",sEByPt_BTagJet_DISCR,"sEByPt_BTagJet_DISCR[50],/F"); */
/*     Tree_->Branch("sM_BTagJet_DISCR",sM_BTagJet_DISCR,"sM_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("sMByIP3d_BTagJet_DISCR",sMByIP3d_BTagJet_DISCR,"sMByIP3d_BTagJet_DISCR[50]/F"); */
/*     Tree_->Branch("sMByPt_BTagJet_DISCR",sMByPt_BTagJet_DISCR,"sMByPt_BTagJet_DISCR[50]/F"); */
     
  }
  
  
  void Initialize(){
    
    irun=-999,ievt=-999,ils=-999;
    auto_cross_section=-999.,external_cross_section=-999., filter_eff=-999.;
    weight=-888.;
    RHO=-999.;

    // PU
    num_PU_vertices=-999;
    PU_BunchCrossing=-999;
    for (int j=0;j<50;j++){
      PU_zpos[j]=-999.;
      PU_sumpT_lowpT[j]=-999.; 
      PU_sumpT_highpT[j]=-999.; 
      PU_ntrks_lowpT[j]=-999; 
      PU_ntrks_highpT[j]=-999; 
    }
    
    // HLT flags
    boolMuonIso=false;
    boolMuonNonIso=false;
    bool2MuonNonIso=false;
    boolElectron=false;
    boolElectronRelaxed=false;
    bool2Electron=false;
    bool2ElectronRelaxed=false;
    boolHLTaccept=false;
    for (int i=0; i<50;i++){
      boolHLTpattern[i]=false;  
    }

    RECO_nMuHLTMatch=0;
    RECO_nMuHLTMatchPAT=0;
    RECO_nMuHLTMatch_asym_PAT=0;
    RECO_nEleHLTMatchPAT=0;
    
    //for (int ii=0;ii<200;ii++){
    //  HLTPathsFired[ii]="";
    //}

    for (int i=0; i<10;i++){
      flag[i]=0;  
    }

    leptonscands2e2mu_= new (CandidateCollection);
    leptonscands2e2murf_= new (CandidateCollection);
    leptonscands4mu_= new (CandidateCollection);
    leptonscands4murf_= new (CandidateCollection);
    leptonscands4e_= new (CandidateCollection);
    leptonscands4erf_= new (CandidateCollection);

    leptonscands_Z0= new (CandidateCollection);
    leptonscands_Z1= new (CandidateCollection);
    leptonscands_Zss0= new (CandidateCollection);
    leptonscands_Zss1= new (CandidateCollection);
    leptonscands_Zcross= new (CandidateCollection);
    leptonscands_DiLep= new (CandidateCollection);
    //leptonscands_MMMM= new (CandidateCollection);
    //leptonscands_EEEE= new (CandidateCollection);
    //leptonscands_EEMM= new (CandidateCollection);
    leptonscands_LLL0= new (CandidateCollection);
    leptonscands_LLL1= new (CandidateCollection);
    leptonscands_LLL2= new (CandidateCollection);
    leptonscands_LLL3= new (CandidateCollection);
    leptonscands_LLLLss0= new (CandidateCollection);
    leptonscands_LLLLss1= new (CandidateCollection);
    leptonscands_LLLLss2= new (CandidateCollection);
    leptonscands_LLLl0= new (CandidateCollection);
    leptonscands_LLLl1= new (CandidateCollection);
    leptonscands_LLLL= new (CandidateCollection);
    leptonscands_MMMM= new (CandidateCollection);
    leptonscands_EEEE= new (CandidateCollection);
    leptonscands_EEMM= new (CandidateCollection);
    leptonscands_MMEE= new (CandidateCollection);
    leptonscands_MMMT= new (CandidateCollection);
    leptonscands_EEET= new (CandidateCollection);
    leptonscands_MMTT= new (CandidateCollection);
    leptonscands_EETT= new (CandidateCollection);
    leptonscands_MMET= new (CandidateCollection);
    leptonscands_EEMT= new (CandidateCollection);
    leptonscands_MMME= new (CandidateCollection);
    leptonscands_EEEM= new (CandidateCollection);
    leptonscands_MT= new (CandidateCollection);
    leptonscands_ET= new (CandidateCollection);
    leptonscands_ME= new (CandidateCollection);

      
    for (int i=0; i<7;i++){
      MC_2e2mu_E[i]=-999.;
      MC_2e2mu_PT[i]=-999.;
      MC_2e2mu_ETA[i]=-999.;
      MC_2e2mu_THETA[i]=-999.;
      MC_2e2mu_PHI[i]=-999.;
      MC_2e2mu_MASS[i]=-999.;
      MC_2e2mu_PDGID[i]=-999.;
      MC_4mu_E[i]=-999.;
      MC_4mu_PT[i]=-999.;
      MC_4mu_ETA[i]=-999.;
      MC_4mu_THETA[i]=-999.;
      MC_4mu_PHI[i]=-999.;
      MC_4mu_MASS[i]=-999.;
      MC_4mu_PDGID[i]=-999.;
      MC_4e_E[i]=-999.;
      MC_4e_PT[i]=-999.;
      MC_4e_ETA[i]=-999.;
      MC_4e_THETA[i]=-999.;
      MC_4e_PHI[i]=-999.;
      MC_4e_MASS[i]=-999.;
      MC_4e_PDGID[i]=-999.;


      MCRF_2e2mu_E[i]=-999.;
      MCRF_2e2mu_PT[i]=-999.;
      MCRF_2e2mu_ETA[i]=-999.;
      MCRF_2e2mu_THETA[i]=-999.;
      MCRF_2e2mu_PHI[i]=-999.;
      MCRF_2e2mu_MASS[i]=-999.;
      MCRF_2e2mu_PDGID[i]=-999.;

      MCRF_4mu_E[i]=-999.;
      MCRF_4mu_PT[i]=-999.;
      MCRF_4mu_ETA[i]=-999.;
      MCRF_4mu_THETA[i]=-999.;
      MCRF_4mu_PHI[i]=-999.;
      MCRF_4mu_MASS[i]=-999.;
      MCRF_4mu_PDGID[i]=-999.;

      MCRF_4e_E[i]=-999.;
      MCRF_4e_PT[i]=-999.;
      MCRF_4e_ETA[i]=-999.;
      MCRF_4e_THETA[i]=-999.;
      MCRF_4e_PHI[i]=-999.;
      MCRF_4e_MASS[i]=-999.;
      MCRF_4e_PDGID[i]=-999.;


      RECOBEST_2e2mu_PT[i]=-999.;
      RECOBEST_2e2mu_ORDERED_PT[i]=-999.;
      RECOBEST_2e2mu_ETA[i]=-999.;
      RECOBEST_2e2mu_THETA[i]=-999.;
      RECOBEST_2e2mu_PHI[i]=-999.;
      RECOBEST_2e2mu_MASS[i]=-999.;
      RECOBEST_2e2mu_isGlobalMu[i]=false;
      RECOBEST_2e2mu_isStandAloneMu[i]=false;
      RECOBEST_2e2mu_isTrackerMu[i]=false;
      RECOBEST_2e2mu_isCaloMu[i]=false;
      RECOBEST_2e2mu_isElectron[i]=false;
 
      RECOBEST_4e_PT[i]=-999.;
      RECOBEST_4e_ORDERED_PT[i]=-999.;
      RECOBEST_4e_ETA[i]=-999.;
      RECOBEST_4e_THETA[i]=-999.;
      RECOBEST_4e_PHI[i]=-999.;
      RECOBEST_4e_MASS[i]=-999.;
      RECOBEST_4e_isElectron[i]=false;

      RECOBEST_4mu_PT[i]=-999.;
      RECOBEST_4mu_ORDERED_PT[i]=-999.;
      RECOBEST_4mu_ETA[i]=-999.;
      RECOBEST_4mu_THETA[i]=-999.;
      RECOBEST_4mu_PHI[i]=-999.;
      RECOBEST_4mu_MASS[i]=-999.;
      RECOBEST_4mu_isGlobalMu[i]=false;
      RECOBEST_4mu_isStandAloneMu[i]=false;
      RECOBEST_4mu_isTrackerMu[i]=false;
      RECOBEST_4mu_isCaloMu[i]=false;

      RECORFBEST_2e2mu_PT[i]=-999.;
      RECORFBEST_2e2mu_ETA[i]=-999.;
      RECORFBEST_2e2mu_THETA[i]=-999.;
      RECORFBEST_2e2mu_PHI[i]=-999.;
      RECORFBEST_2e2mu_MASS[i]=-999.;

      RECORFBEST_4mu_PT[i]=-999.;
      RECORFBEST_4mu_ETA[i]=-999.;
      RECORFBEST_4mu_THETA[i]=-999.;
      RECORFBEST_4mu_PHI[i]=-999.;
      RECORFBEST_4mu_MASS[i]=-999.;

      RECORFBEST_4e_PT[i]=-999.;
      RECORFBEST_4e_ETA[i]=-999.;
      RECORFBEST_4e_THETA[i]=-999.;
      RECORFBEST_4e_PHI[i]=-999.;
      RECORFBEST_4e_MASS[i]=-999.;

    }

   
    for (int i=0; i<4;i++){
      MCRF_2e2mu_cosTheta1_spin[i]=-999.;
      MCRF_2e2mu_cosTheta2_spin[i]=-999.;
      MCRF_2e2mu_cosThetaStar_spin[i]=-999.;
      MCRF_2e2mu_Phi_spin[i]=-999.;
      MCRF_2e2mu_Phi1_spin[i]=-999.;
      MCRF_2e2mu_Phi2_spin[i]=-999.;
      MCRF_2e2mu_phi1RF_spin[i]=-999.;
      MCRF_2e2mu_phi2RF_spin[i]=-999.;			
      
      RECORFBEST_2e2mu_cosTheta1_spin[i]=-999.;
      RECORFBEST_2e2mu_cosTheta2_spin[i]=-999.;
      RECORFBEST_2e2mu_cosThetaStar_spin[i]=-999.;
      RECORFBEST_2e2mu_Phi_spin[i]=-999.;
      RECORFBEST_2e2mu_Phi1_spin[i]=-999.;
      RECORFBEST_2e2mu_Phi2_spin[i]=-999.;
      RECORFBEST_2e2mu_phi1RF_spin[i]=-999.;
      RECORFBEST_2e2mu_phi2RF_spin[i]=-999.;

      MCRF_4mu_cosTheta1_spin[i]=-999.;
      MCRF_4mu_cosTheta2_spin[i]=-999.;
      MCRF_4mu_cosThetaStar_spin[i]=-999.;
      MCRF_4mu_Phi_spin[i]=-999.;
      MCRF_4mu_Phi1_spin[i]=-999.;
      MCRF_4mu_Phi2_spin[i]=-999.;
      MCRF_4mu_phi1RF_spin[i]=-999.;
      MCRF_4mu_phi2RF_spin[i]=-999.;			
      
      RECORFBEST_4mu_cosTheta1_spin[i]=-999.;
      RECORFBEST_4mu_cosTheta2_spin[i]=-999.;
      RECORFBEST_4mu_cosThetaStar_spin[i]=-999.;
      RECORFBEST_4mu_Phi_spin[i]=-999.;
      RECORFBEST_4mu_Phi1_spin[i]=-999.;
      RECORFBEST_4mu_Phi2_spin[i]=-999.;
      RECORFBEST_4mu_phi1RF_spin[i]=-999.;
      RECORFBEST_4mu_phi2RF_spin[i]=-999.;

      MCRF_4e_cosTheta1_spin[i]=-999.;
      MCRF_4e_cosTheta2_spin[i]=-999.;
      MCRF_4e_cosThetaStar_spin[i]=-999.;
      MCRF_4e_Phi_spin[i]=-999.;
      MCRF_4e_Phi1_spin[i]=-999.;
      MCRF_4e_Phi2_spin[i]=-999.;
      MCRF_4e_phi1RF_spin[i]=-999.;
      MCRF_4e_phi2RF_spin[i]=-999.;			
      
      RECORFBEST_4e_cosTheta1_spin[i]=-999.;
      RECORFBEST_4e_cosTheta2_spin[i]=-999.;
      RECORFBEST_4e_cosThetaStar_spin[i]=-999.;
      RECORFBEST_4e_Phi_spin[i]=-999.;
      RECORFBEST_4e_Phi1_spin[i]=-999.;
      RECORFBEST_4e_Phi2_spin[i]=-999.;
      RECORFBEST_4e_phi1RF_spin[i]=-999.;
      RECORFBEST_4e_phi2RF_spin[i]=-999.;


    }
    
    genmet=-999.,calomet=-999.;  
    RECO_NMU=0,RECO_NELE=0;
    RECO_NTRACK=0;
    
    
    RECO_NPHOT=0,RECO_NJET=0,RECO_NVTX=0;
    
		
    for(int ivtx=0;ivtx<10;ivtx++) {
      RECO_VERTEX_x[ivtx] = -999;
      RECO_VERTEX_y[ivtx] = -999;
      RECO_VERTEX_z[ivtx] = -999;
      RECO_VERTEXPROB[ivtx]=-999.;
      RECO_VERTEX_ndof[ivtx] = -999;
      RECO_VERTEX_chi2[ivtx] = -999;
      RECO_VERTEX_ntracks[ivtx]=-999;
      for (int i=0;i<30;i++){
	RECO_VERTEX_TRACK_PT[ivtx][i]=-999;
      }
    }
    
    RECO_PFJET_N = 0;
    for (int ijets=0;ijets<30;ijets++) {
      RECO_PFJET_CHARGE[ijets]   = -999.; 
      RECO_PFJET_ET[ijets]   = -999.; 
      RECO_PFJET_PT[ijets]  = -999.; 
      RECO_PFJET_ETA[ijets] = -999.; 
      RECO_PFJET_PHI[ijets] = -999.;
    }
    
    calomet=-999.;
    calometopt=-999.;
    calometoptnohf=-999.;
    calometoptnohfho=-999.;
    calometoptho=-999.;
    calometnohf=-999.;
    calometnohfho=-999.;
    calometho=-999.;
    pfmet=-999.;
    pfmet_x=-999.;
    pfmet_y=-999.;
    pfmet_phi=-999.;
    pfmet_theta=-999.;
    htmetic5=-999.;
    htmetkt4=-999.;
    htmetkt6=-999.;
    htmetsc5=-999.;
    htmetsc7=-999.;
    tcmet=-999.;
    jescormetic5=-999.;
    jescormetkt4=-999.;
    jescormetkt6=-999.;
    jescormetsc5=-999.;
    jescormetsc7=-999.;
    cormetmuons=-999.;
    
    BeamSpot_X=-999.;
    BeamSpot_Y=-999.;
    BeamSpot_Z=-999.;
    
    
    for (int i=0; i<4;i++){
      isolation[i]=-999.;
      vertexing[i]=-999.;
      ftsigma[i]=-999.;
      ftsigmalag[i]=-999.;
      gdX[i]=-999.;
      gdY[i]=-999.;
      gdZ[i]=-999.;
      gdlagX[i]=-999.;
      gdlagY[i]=-999.;
      gdlagZ[i]=-999.;
      gdlagProb[i]=-999.;
      gdlagNdof[i]=-999.;
      ftsigmaMMMM[i]=-999.;
      ftsigmalagMMMM[i]=-999.;
      gdXMMMM[i]=-999.;
      gdYMMMM[i]=-999.;
      gdZMMMM[i]=-999.;
      gdlagXMMMM[i]=-999.;
      gdlagYMMMM[i]=-999.;
      gdlagZMMMM[i]=-999.;
      gdlagProbMMMM[i]=-999.;
      gdlagNdofMMMM[i]=-999.;
      ftsigmaEEEE[i]=-999.;
      ftsigmalagEEEE[i]=-999.;
      gdXEEEE[i]=-999.;
      gdYEEEE[i]=-999.;
      gdZEEEE[i]=-999.;
      gdlagXEEEE[i]=-999.;
      gdlagYEEEE[i]=-999.;
      gdlagZEEEE[i]=-999.;
      gdlagProbEEEE[i]=-999.;
      gdlagNdofEEEE[i]=-999.;

      StdFitVertexX[i]=-999.;
      StdFitVertexY[i]=-999.;
      StdFitVertexZ[i]=-999.;
      StdFitVertexChi2r[i]=-999.;
      StdFitVertexProb[i]=-999.;
      KinFitVertexX[i]=-999.;
      KinFitVertexY[i]=-999.;
      KinFitVertexZ[i]=-999.;
      KinFitVertexChi2r[i]=-999.;
      KinFitVertexProb[i]=-999.;

      StdFitVertexXMMMM[i]=-999.;
      StdFitVertexYMMMM[i]=-999.;
      StdFitVertexZMMMM[i]=-999.;
      StdFitVertexChi2rMMMM[i]=-999.;
      StdFitVertexProbMMMM[i]=-999.;
      KinFitVertexXMMMM[i]=-999.;
      KinFitVertexYMMMM[i]=-999.;
      KinFitVertexZMMMM[i]=-999.;
      KinFitVertexChi2rMMMM[i]=-999.;
      KinFitVertexProbMMMM[i]=-999.;

      StdFitVertexXEEEE[i]=-999.;
      StdFitVertexYEEEE[i]=-999.;
      StdFitVertexZEEEE[i]=-999.;
      StdFitVertexChi2rEEEE[i]=-999.;
      StdFitVertexProbEEEE[i]=-999.;
      KinFitVertexXEEEE[i]=-999.;
      KinFitVertexYEEEE[i]=-999.;
      KinFitVertexZEEEE[i]=-999.;
      KinFitVertexChi2rEEEE[i]=-999.;
      KinFitVertexProbEEEE[i]=-999.;

    }
    
    for (int i=0; i<20;i++){
      StdFitVertexChi2rMMM[i]=-999.;
      StdFitVertexProbMMM[i]=-999.;
      StdFitVertexChi2rMME[i]=-999.;
      StdFitVertexProbMME[i]=-999.;
      StdFitVertexChi2rEEE[i]=-999.;
      StdFitVertexProbEEE[i]=-999.;
      StdFitVertexChi2rMEE[i]=-999.;
      StdFitVertexProbMEE[i]=-999.;
    }
    

    for (int i=0; i<40;i++){
      StdFitVertexChi2rDiLep[i]=-999.;
      StdFitVertexProbDiLep[i]=-999.;
    }

    for (int i=0; i<2;i++){
      isoele[i]=-999.;
      isomu[i]=-999.;
    }
    
    /* for (int j=0; j<7;j++){ */
/*       for (int i=0; i<50;i++){ */
/* 	RECO_MMMM_MASS[j][i]=-999.; */
/* 	RECO_MMMM_PT[j][i]=-999.; */
/* 	RECO_EEMM_MASS[j][i]=-999.; */
/* 	RECO_EEMM_PT[j][i]=-999.; */
/* 	RECO_EEEE_MASS[j][i]=-999.; */
/* 	RECO_EEEE_PT[j][i]=-999.; */

/* 	RECO_MMMM_ETA[j][i]=-999.; */
/* 	RECO_MMMM_PHI[j][i]=-999.; */
/* 	RECO_EEMM_ETA[j][i]=-999.; */
/* 	RECO_EEMM_PHI[j][i]=-999.; */
/* 	RECO_EEEE_ETA[j][i]=-999.; */
/* 	RECO_EEEE_PHI[j][i]=-999.; */

/*       } */
/*     } */

    for (int j=0; j<9;j++){
      for (int i=0; i<200;i++){
	RECO_MMTT_MASS[j][i]=-999.;
	RECO_MMTT_PT[j][i]=-999.;
	RECO_MMTT_ETA[j][i]=-999.;
	RECO_MMTT_PHI[j][i]=-999.;
	RECO_MMTT_CHARGE[j][i]=-999.;
	
	RECO_EETT_MASS[j][i]=-999.;
	RECO_EETT_PT[j][i]=-999.;
	RECO_EETT_ETA[j][i]=-999.;
	RECO_EETT_PHI[j][i]=-999.;
	RECO_EETT_CHARGE[j][i]=-999.;
      }
    }
    for (int j=0; j<9;j++){
      for (int i=0; i<100;i++){
	RECO_EEET_MASS[j][i]=-999.;
	RECO_EEET_PT[j][i]=-999.;
	RECO_EEET_ETA[j][i]=-999.;
	RECO_EEET_PHI[j][i]=-999.;
	RECO_EEET_CHARGE[j][i]=-999.;
      }
    }
    
    for (int j=0; j<9;j++){
      for (int i=0; i<100;i++){
	RECO_MMMM_MASS[j][i]=-999.;
	RECO_MMMM_PT[j][i]=-999.;
	RECO_MMMM_ETA[j][i]=-999.;
	RECO_MMMM_PHI[j][i]=-999.;
	RECO_MMMM_CHARGE[j][i]=-999.;
	
	RECO_EEMM_MASS[j][i]=-999.;
	RECO_EEMM_PT[j][i]=-999.;
	RECO_EEMM_ETA[j][i]=-999.;
	RECO_EEMM_PHI[j][i]=-999.;
	RECO_EEMM_CHARGE[j][i]=-999.;

	RECO_MMEE_MASS[j][i]=-999.;
	RECO_MMEE_PT[j][i]=-999.;
	RECO_MMEE_ETA[j][i]=-999.;
	RECO_MMEE_PHI[j][i]=-999.;
	RECO_MMEE_CHARGE[j][i]=-999.;

	RECO_EEEE_MASS[j][i]=-999.;
	RECO_EEEE_PT[j][i]=-999.;
	RECO_EEEE_ETA[j][i]=-999.;
	RECO_EEEE_PHI[j][i]=-999.;
	RECO_EEEE_CHARGE[j][i]=-999.;

	RECO_MMMT_MASS[j][i]=-999.;
	RECO_MMMT_PT[j][i]=-999.;
	RECO_MMMT_ETA[j][i]=-999.;
	RECO_MMMT_PHI[j][i]=-999.;
	RECO_MMMT_CHARGE[j][i]=-999.;

	

	RECO_MMET_MASS[j][i]=-999.;
	RECO_MMET_PT[j][i]=-999.;
	RECO_MMET_ETA[j][i]=-999.;
	RECO_MMET_PHI[j][i]=-999.;
	RECO_MMET_CHARGE[j][i]=-999.;

	RECO_EEMT_MASS[j][i]=-999.;
	RECO_EEMT_PT[j][i]=-999.;
	RECO_EEMT_ETA[j][i]=-999.;
	RECO_EEMT_PHI[j][i]=-999.;
	RECO_EEMT_CHARGE[j][i]=-999.;

	RECO_MMME_MASS[j][i]=-999.;
	RECO_MMME_PT[j][i]=-999.;
	RECO_MMME_ETA[j][i]=-999.;
	RECO_MMME_PHI[j][i]=-999.;
	RECO_MMME_CHARGE[j][i]=-999.;

	RECO_EEEM_MASS[j][i]=-999.;
	RECO_EEEM_PT[j][i]=-999.;
	RECO_EEEM_ETA[j][i]=-999.;
	RECO_EEEM_PHI[j][i]=-999.;
	RECO_EEEM_CHARGE[j][i]=-999.;
      }
    }
    
    for (int i=0; i<100;i++){
      HPSTAU_discByLER[i] = -999.;
      HPSTAU_discByMER[i] = -999.;
      HPSTAU_discByTER[i] = -999.;
      HPSTAU_discByLMR[i] = -999.;
      HPSTAU_discByTMR[i] = -999.;
      HPSTAU_discByDecayF[i] = -999.;
      HPSTAU_discLooseIso[i] = -999.;
      HPSTAU_discVLooseIso[i] = -999.;
      HPSTAU_discMediumIso[i] = -999.;
      HPSTAU_discTightIso[i] = -999.;
      HPSTAU_discLooseChargedIso[i] = -999.;
      HPSTAU_discVLooseChargedIso[i] = -999.;
      HPSTAU_discMediumChargedIso[i] = -999.;
      HPSTAU_discTightChargedIso[i] = -999.;
      HPSTAU_discLooseIsoDB[i] = -999.;
      HPSTAU_discVLooseIsoDB[i] = -999.;
      HPSTAU_discMediumIsoDB[i] = -999.;
      HPSTAU_discTightIsoDB[i] = -999.;
      HPSTAU_discLooseCombinedIsoDB[i] = -999.;
      HPSTAU_discVLooseCombinedIsoDB[i] = -999.;
      HPSTAU_discMediumCombinedIsoDB[i] = -999.;
      HPSTAU_discTightCombinedIsoDB[i] = -999.;
      HPSTAU_MASS[i] = -999.;
      HPSTAU_PT[i] = -999.;
      HPSTAU_ETA[i] = -999.;
      HPSTAU_PHI[i] = -999.;
      HPSTAU_CHARGE[i] = -999.;
    }

    for (int i=0; i<100;i++){
      EleId95relIso[i]=-999.;
      EleId90relIso[i]=-999.;
      EleId85relIso[i]=-999.;      
      EleId80relIso[i]=-999.;
    }


    for (int i=0; i<50;i++){

      RECO_ZMM_MASS[i]=-999.;
      RECO_ZMM_PT[0][i]=-999.;
      RECO_ZMM_PT[1][i]=-999.;
      RECO_ZMM_PT[2][i]=-999.;
      RECO_ZMM_ETA[0][i]=-999.;
      RECO_ZMM_ETA[1][i]=-999.;
      RECO_ZMM_ETA[2][i]=-999.;
      RECO_ZMM_PHI[0][i]=-999.;
      RECO_ZMM_PHI[1][i]=-999.;
      RECO_ZMM_PHI[2][i]=-999.;
      
      RECO_ZEE_MASS[i]=-999.;
      RECO_ZEE_PT[0][i]=-999.;
      RECO_ZEE_PT[1][i]=-999.;
      RECO_ZEE_PT[2][i]=-999.;
      RECO_ZEE_ETA[0][i]=-999.;
      RECO_ZEE_ETA[1][i]=-999.;
      RECO_ZEE_ETA[2][i]=-999.;
      RECO_ZEE_PHI[0][i]=-999.;
      RECO_ZEE_PHI[1][i]=-999.;
      RECO_ZEE_PHI[2][i]=-999.;
      
      RECO_ZMMss_MASS[i]=-999.;
      RECO_ZMMss_PT[0][i]=-999.;
      RECO_ZMMss_PT[1][i]=-999.;
      RECO_ZMMss_PT[2][i]=-999.;
      RECO_ZMMss_ETA[0][i]=-999.;
      RECO_ZMMss_ETA[1][i]=-999.;
      RECO_ZMMss_ETA[2][i]=-999.;
      RECO_ZMMss_PHI[0][i]=-999.;
      RECO_ZMMss_PHI[1][i]=-999.;
      RECO_ZMMss_PHI[2][i]=-999.;
      

      RECO_ZEEss_MASS[i]=-999.;
      RECO_ZEEss_PT[0][i]=-999.;
      RECO_ZEEss_PT[1][i]=-999.;
      RECO_ZEEss_PT[2][i]=-999.;
      RECO_ZEEss_ETA[0][i]=-999.;
      RECO_ZEEss_ETA[1][i]=-999.;
      RECO_ZEEss_ETA[2][i]=-999.;
      RECO_ZEEss_PHI[0][i]=-999.;
      RECO_ZEEss_PHI[1][i]=-999.;
      RECO_ZEEss_PHI[2][i]=-999.;


      RECO_ZEM_MASS[i]=-999.;
      RECO_ZEM_PT[0][i]=-999.;
      RECO_ZEM_PT[1][i]=-999.;
      RECO_ZEM_PT[2][i]=-999.;
      RECO_ZEM_ETA[0][i]=-999.;
      RECO_ZEM_ETA[1][i]=-999.;
      RECO_ZEM_ETA[2][i]=-999.;
      RECO_ZEM_PHI[0][i]=-999.;
      RECO_ZEM_PHI[1][i]=-999.;
      RECO_ZEM_PHI[2][i]=-999.;

      RECO_DiLep_MASS[i]=-999.;
      RECO_DiLep_PT[0][i]=-999.;
      RECO_DiLep_PT[1][i]=-999.;
      RECO_DiLep_PT[2][i]=-999.;
      RECO_DiLep_ETA[0][i]=-999.;
      RECO_DiLep_ETA[1][i]=-999.;
      RECO_DiLep_ETA[2][i]=-999.;
      RECO_DiLep_PHI[0][i]=-999.;
      RECO_DiLep_PHI[1][i]=-999.;
      RECO_DiLep_PHI[2][i]=-999.;
     
      RECO_LLL0_MASS[i]=-999;
      RECO_LLL1_MASS[i]=-999;
      RECO_LLL2_MASS[i]=-999;
      RECO_LLL3_MASS[i]=-999;
      for (int j=0; j<4;j++){
	RECO_LLL0_PT[j][i]=-999;
	RECO_LLL1_PT[j][i]=-999;
	RECO_LLL2_PT[j][i]=-999;
	RECO_LLL3_PT[j][i]=-999;
      }
      
      RECO_LLLL0ss_MASS[i]=-999;
      RECO_LLLL1ss_MASS[i]=-999;
      RECO_LLLL2ss_MASS[i]=-999;

      RECO_LLLl0_MASS[i]=-999;
      RECO_LLLl1_MASS[i]=-999;
      
      for (int j=0; j<4;j++){
	RECO_LLLL1ss_PT[j][i]=-999;
	RECO_LLLL1ss_PT[j][i]=-999;
	RECO_LLLL2ss_PT[j][i]=-999;
	RECO_LLLl0_PT[j][i]=-999;
	RECO_LLLl1_PT[j][i]=-999;
      }
    }
    
    
    for (int j=0; j<7;j++){
      for (int i=0; i<100;i++){
	RECO_LLLL_MASS[j][i]=-999.;
	RECO_LLLL_PT[j][i]=-999.;
	RECO_LLLL_ETA[j][i]=-999.;
	RECO_LLLL_PHI[j][i]=-999.;
      }
    }

    for (int i=0; i<100;i++){
      RECOELE_E[i]     = -999.;
      RECOELE_PT[i]=-999.;
      RECOELE_P[i]=-999.;
      RECOELE_PHI[i]=-999.;
      RECOELE_ETA[i]=-999.;
      RECOELE_THETA[i]=-999.;
      RECOELE_MASS[i]=-999.;
      RECOELE_TMASS[i]=-999.;
      RECOELE_QUALITY[i]=-999.;			
      RECOELE_CHARGE[i]=-999.;
      
      // Core attributes
      RECOELE_isEcalDriven[i]    = false;
      RECOELE_isTrackerDriven[i] = false;
      RECOELE_gsftrack_chi2[i]   = -999;
      RECOELE_gsftrack_dxyB[i]   = -999;
      RECOELE_gsftrack_dxy[i]    = -999;
      RECOELE_gsftrack_dxyError[i]    = -999;
      RECOELE_gsftrack_dzB[i]    = -999;
      RECOELE_gsftrack_dz[i]     = -999;
      RECOELE_gsftrack_dzError[i]   = -999;
      RECOELE_gsftrack_losthits[i]  = -999;
      RECOELE_gsftrack_validhits[i] = -999;
      RECOELE_gsftrack_expected_inner_hits[i]=-999;
      RECOELE_scl_E[i]   = -999;
      RECOELE_scl_Et[i]  = -999;
      RECOELE_scl_Eta[i] = -999;
      RECOELE_scl_Phi[i] = -999;
      // Track-Cluster matching attributes
      RECOELE_ep[i]             = -999.;
      RECOELE_eSeedp[i]         = -999.;
      RECOELE_eSeedpout[i]      = -999.;
      RECOELE_eElepout[i]       = -999.;
      //
      RECOELE_deltaEtaIn[i]     = -999.;
      RECOELE_deltaEtaSeed[i]   = -999.;
      RECOELE_deltaEtaEle[i]    = -999.;
      RECOELE_deltaPhiIn[i]     = -999.;
      RECOELE_deltaPhiSeed[i]   = -999.;
      RECOELE_deltaPhiEle[i]    = -999.;
      // Fiducial flags 
      RECOELE_isbarrel[i]    = 0;
      RECOELE_isendcap[i]    = 0;
      RECOELE_isEBetaGap[i]  = 0;
      RECOELE_isEBphiGap[i]  = 0;
      RECOELE_isEEdeeGap[i]  = 0;
      RECOELE_isEEringGap[i] = 0;
      // Shower shape
      RECOELE_sigmaIetaIeta[i] = -999.;
      RECOELE_sigmaEtaEta[i]   = -999.;
      RECOELE_e15[i]           = -999.;
      RECOELE_e25max[i]        = -999.;
      RECOELE_e55[i]           = -999.;
      RECOELE_he[i]            = -999.;
      // Particle flow
      RECOELE_mva[i] = -999.; 
      // Brem & Classifaction
      RECOELE_fbrem[i]  = -999.; 
      RECOELE_nbrems[i] = -999;
      RECOELE_Class[i]  = -999;
      RECOELE_fbrem_mode[i]=-999.;
      RECOELE_fbrem_mean[i]=-999.;
      // isolation
      RECOELE_TRACKISO[i]=-999.; 
      RECOELE_HCALISO[i]=-999.;
      RECOELE_ECALISO[i]=-999.;
      RECOELE_X[i]=-999.;
      RECOELE_EGMTRACKISO[i]=-999.; 
      RECOELE_EGMHCALISO[i]=-999.;
      RECOELE_EGMECALISO[i]=-999.;
      RECOELE_EGMX[i]=-999.;
      RECOELE_IP[i]=-9999.;
      RECOELE_SIP[i]=-9999.;
      RECOELE_IPERROR[i]=-9999.;
      RECOELE_IP_KF[i]=-999.;
      RECOELE_SIP_KF[i]=-999.;
      RECOELE_IPERROR_KF[i]=-999.;
      RECOELE_SIP_GD[i]=-999.;
      RECOELE_SIP_GDEEEE[i]=-999.;
      RECOELE_SIP_Std[i]=-999.;
      RECOELE_SIP_StdEEEE[i]=-999.;
      RECOELE_SIP_Kin[i]=-999.;
      RECOELE_SIP_KinEEEE[i]=-999.;
      RECOELE_STIP[i]=-999.;
      RECOELE_TIP[i]=-999.;
      RECOELE_TIPERROR[i]=-999.;
      RECOELE_SLIP[i]=-999.;
      RECOELE_LIP[i]=-999.;
      RECOELE_LIPERROR[i]=-999.;
      RECOELE_eleID[i][0]=-999.;
      RECOELE_eleID[i][1]=-999.;
      ele_sclRawE[i]=-999. ;
      ele_sclX[i]=-999.; 
      ele_sclY[i]=-999.; 
      ele_sclZ[i]=-999.;
      ele_seedSubdet1[i]=-999.;
      ele_seedDphi1[i]=-999.; 
      ele_seedDrz1[i]=-999.;
      ele_seedSubdet2[i]=-999.;
      ele_seedDphi2[i]=-999.; 
      ele_seedDrz2[i]=-999.;
      ele_severityLevelSeed[i]=-999.;
      ele_severityLevelClusters[i]=-999.; 
      ele_outOfTimeSeed[i]=-999.; 
      ele_outOfTimeClusters[i]=-999.;
      ele_e2overe9[i]=-999. ;
      ele_eidVeryLoose[i]=-999.; 
      ele_eidLoose[i]=-999.; 
      ele_eidMedium[i]=-999.; 
      ele_eidTight[i]=-999. ;
      ele_eidHZZVeryLoose[i]=-999.; 
      ele_eidHZZLoose[i]=-999.; 
      ele_eidHZZMedium[i]=-999.; 
      ele_eidHZZTight[i]=-999. ;

      RECOELEBEST_2e2mu_MATCHED[i]=0;
      RECOELEBEST_4e_MATCHED[i]=0;
      
      RECOELE_ZEE_MATCHED[i]=0;
      RECOELE_ZssEE_MATCHED[i]=0;
      RECOELE_ZEM_MATCHED[i]=0;
      RECOELE_EEEE_MATCHED[i]=0; 
      RECOELE_EEMM_MATCHED[i]=0; 
      RECOELE_LLL0_MATCHED[i]=0;
      RECOELE_LLL1_MATCHED[i]=0;
      RECOELE_LLL2_MATCHED[i]=0;
      RECOELE_LLL3_MATCHED[i]=0;
      RECOELE_LLLLss0_MATCHED[i]=0;
      RECOELE_LLLLss1_MATCHED[i]=0;
      RECOELE_LLLLss2_MATCHED[i]=0;
      RECOELE_LLLl0_MATCHED[i]=0;
      RECOELE_LLLl1_MATCHED[i]=0;
      RECOELE_LLLL_MATCHED[i]=0;

      // Conversion
      ConvMapDist[i]=-999.;
      ConvMapDcot[i]=-999.;
      
      // Muon block
      RECOMU_PT_MuHLTMatch[i]=-999.;
      RECOMU_PT_MuHLTMatchPAT[i]=-999.;
      RECOMU_PT_MuHLTMatch_asym_PAT[i]=-999.;
      RECOELE_PT_EleHLTMatchPAT[i]=-999.;

      RECOMU_isGlobalMu[i]=false;
      RECOMU_isStandAloneMu[i]=false;
      RECOMU_isTrackerMu[i]=false;
      RECOMU_isCaloMu[i]=false;

      RECOMU_E[i]=-999.;
      RECOMU_PT[i]=-999.;
      RECOMU_P[i]=-999.;
      RECOMU_PHI[i]=-999.;
      RECOMU_ETA[i]=-999.;
      RECOMU_THETA[i]=-999.;
      RECOMU_MASS[i]=-999.;
      RECOMU_TMASS[i]=-999.;
      RECOMU_QUALITY[i]=-999.;
      
      RECOMU_CHARGE[i]=-999.;
      RECOMU_TRACKISO[i]=-999.; 
      RECOMU_ECALISO[i]=-999.;
      RECOMU_HCALISO[i]=-999.;
      RECOMU_X[i]=-999.;
      
      // Particle Flow Isolation
      RECOMU_ALLDEPOSITS[i]=-999.;
      RECOMU_ChargedDEPOSITS[i]=-999.;
      RECOMU_NeutralDEPOSITS[i]=-999.;
      RECOMU_GammaDEPOSITS[i]=-999.;
      RECOMU_PUDEPOSITS[i]=-999.;
      RECOMU_PULowDEPOSITS[i]=-999.;
      RECOMU_PFX_DB[i]=-999.;
      
      RECOELE_ALLDEPOSITS[i]=-999.;
      RECOELE_ChargedDEPOSITS[i]=-999.;
      RECOELE_NeutralDEPOSITS[i]=-999.;
      RECOELE_GammaDEPOSITS[i]=-999.;
      RECOELE_PUDEPOSITS[i]=-999.;
      RECOELE_PULowDEPOSITS[i]=-999.;
      RECOELE_PFX_DB[i]=-999.;
      
      RECOMU_IP[i]=-9999.;
      RECOMU_SIP[i]=-9999.;
      RECOMU_IPERROR[i]=-9999.;
      RECOMU_IP_KF[i]=-999.;
      RECOMU_SIP_KF[i]=-999.;
      RECOMU_IPERROR_KF[i]=-999.;

      RECOMU_SIP_GD[i]=-999.;
      RECOMU_SIP_GDMMMM[i]=-999.;
      RECOMU_SIP_Std[i]=-999.;
      RECOMU_SIP_StdMMMM[i]=-999.;
      RECOMU_SIP_Kin[i]=-999.;
      RECOMU_SIP_KinMMMM[i]=-999.;
      RECOMU_STIP[i]=-999.;
      RECOMU_TIP[i]=-999.;
      RECOMU_TIPERROR[i]=-999.;
      RECOMU_SLIP[i]=-999.;
      RECOMU_LIP[i]=-999.;
      RECOMU_LIPERROR[i]=-999.;
      
      RECOMUBEST_2e2mu_MATCHED[i]=0;
      RECOMUBEST_4mu_MATCHED[i]=0;
      RECOMU_ZMM_MATCHED[i]=0;
      RECOMU_ZssMM_MATCHED[i]=0;
      RECOMU_ZEM_MATCHED[i]=0;
      RECOMU_MMMM_MATCHED[i]=0;
      RECOMU_EEMM_MATCHED[i]=0; 
      RECOMU_LLL0_MATCHED[i]=0;
      RECOMU_LLL1_MATCHED[i]=0;
      RECOMU_LLL2_MATCHED[i]=0;
      RECOMU_LLL3_MATCHED[i]=0;      
      RECOMU_LLLLss0_MATCHED[i]=0;
      RECOMU_LLLLss1_MATCHED[i]=0;
      RECOMU_LLLLss2_MATCHED[i]=0;
      RECOMU_LLLl0_MATCHED[i]=0;
      RECOMU_LLLl1_MATCHED[i]=0;
      RECOMU_LLLL_MATCHED[i]=0;

      RECOMU_numberOfMatches[i]=-999;
      RECOMU_caloCompatibility[i]=-999.;
      RECOMU_segmentCompatibility[i]=-999.;
      RECOMU_glbmuPromptTight[i]=false;
      
      RECOMU_mutrkDxy[i]=-999.;
      RECOMU_mutrkDxyError[i]=-999.;
      RECOMU_mutrkDxyB[i]=-999.;
      RECOMU_mutrkDz[i]=-999.;
      RECOMU_mutrkDzError[i]=-999.;
      RECOMU_mutrkDzB[i]=-999.;
      RECOMU_mutrkChi2PerNdof[i]=-999.;
      RECOMU_mutrkCharge[i]=-999.;
      RECOMU_mutrkNHits[i]=-999.;
      RECOMU_mutrkNPixHits[i]=-999.;
      RECOMU_mutrkNMuonHits[i]=-999.;
      RECOMU_mutrkNStripHits[i]=-999.;
      RECOMU_trkmuArbitration[i]=false;
      RECOMU_trkmu2DCompatibilityLoose[i]=false;
      RECOMU_trkmu2DCompatibilityTight[i]=false;
      RECOMU_trkmuOneStationLoose[i]=false;
      RECOMU_trkmuOneStationTight[i]=false;
      RECOMU_trkmuLastStationLoose[i]=false;
      RECOMU_trkmuLastStationTight[i]=false;
      RECOMU_trkmuOneStationAngLoose[i]=false;
      RECOMU_trkmuOneStationAngTight[i]=false;
      RECOMU_trkmuLastStationAngLoose[i]=false;
      RECOMU_trkmuLastStationAngTight[i]=false;
      RECOMU_trkmuLastStationOptimizedLowPtLoose[i]=false;
      RECOMU_trkmuLastStationOptimizedLowPtTight[i]=false;


    }
    
    for (int i=0; i<200;i++){
      RECO_TRACK_PT[i]=-999.;
      RECO_TRACK_ETA[i]=-999.;
      RECO_TRACK_PHI[i]=-999.;
      RECO_TRACK_CHI2[i]=-999.;
      RECO_TRACK_CHI2RED[i]=-999.;
      RECO_TRACK_CHI2PROB[i]=-999.;
      RECO_TRACK_NHITS[i]=0;
      RECO_TRACK_DXY[i]=-999.;
      RECO_TRACK_DXYERR[i]=-999.;
      RECO_TRACK_DZ[i]=-999.;
      RECO_TRACK_DZERR[i]=-999.;   			       					
    }

    for (int i=0; i<50;i++){
      tCHighEff_BTagJet_PT[i]=-999.;
      tCHighPur_BTagJet_PT[i]=-999.;
      jPHighEff_BTagJet_PT[i]=-999.;
      jBP_BTagJet_PT[i]=-999.;
      sSVHighEff_BTagJet_PT[i]=-999.;
      sSVHighPur_BTagJet_PT[i]=-999.;
      cSV_BTagJet_PT[i]=-999.;
      cSVMVA_BTagJet_PT[i]=-999.;
      sEByIP3d_BTagJet_PT[i]=-999.;
      sEByPt_BTagJet_PT[i]=-999.;
      sM_BTagJet_PT[i]=-999.;
      sMByIP3d_BTagJet_PT[i]=-999.;
      sMByPt_BTagJet_PT[i]=-999.;

      tCHighEff_BTagJet_ETA[i]=-999.;
      tCHighPur_BTagJet_ETA[i]=-999.;
      jPHighEff_BTagJet_ETA[i]=-999.;
      jBP_BTagJet_ETA[i]=-999.;
      sSVHighEff_BTagJet_ETA[i]=-999.;
      sSVHighPur_BTagJet_ETA[i]=-999.;
      cSV_BTagJet_ETA[i]=-999.;
      cSVMVA_BTagJet_ETA[i]=-999.;
      sEByIP3d_BTagJet_ETA[i]=-999.;
      sEByPt_BTagJet_ETA[i]=-999.;
      sM_BTagJet_ETA[i]=-999.;
      sMByIP3d_BTagJet_ETA[i]=-999.;
      sMByPt_BTagJet_ETA[i]=-999.;

      tCHighEff_BTagJet_PHI[i]=-999.;
      tCHighPur_BTagJet_PHI[i]=-999.;
      jPHighEff_BTagJet_PHI[i]=-999.;
      jBP_BTagJet_PHI[i]=-999.;
      sSVHighEff_BTagJet_PHI[i]=-999.;
      sSVHighPur_BTagJet_PHI[i]=-999.;
      cSV_BTagJet_PHI[i]=-999.;
      cSVMVA_BTagJet_PHI[i]=-999.;
      sEByIP3d_BTagJet_PHI[i]=-999.;
      sEByPt_BTagJet_PHI[i]=-999.;
      sM_BTagJet_PHI[i]=-999.;
      sMByIP3d_BTagJet_PHI[i]=-999.;
      sMByPt_BTagJet_PHI[i]=-999.;

      tCHighEff_BTagJet_DISCR[i]=-999.;
      tCHighPur_BTagJet_DISCR[i]=-999.;
      jPHighEff_BTagJet_DISCR[i]=-999.;
      jBP_BTagJet_DISCR[i]=-999.;
      sSVHighEff_BTagJet_DISCR[i]=-999.;
      sSVHighPur_BTagJet_DISCR[i]=-999.;
      cSV_BTagJet_DISCR[i]=-999.;
      cSVMVA_BTagJet_DISCR[i]=-999.;
      sEByIP3d_BTagJet_DISCR[i]=-999.;
      sEByPt_BTagJet_DISCR[i]=-999.;
      sM_BTagJet_DISCR[i]=-999.;
      sMByIP3d_BTagJet_DISCR[i]=-999.;
      sMByPt_BTagJet_DISCR[i]=-999.;

    }
    
  }
  
  void fillPU(const edm::Event& iEvent){
      edm::Handle<vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByLabel(PileupSrc_.label(), PupInfo);

      int i=0;
      for( vector<PileupSummaryInfo>::const_iterator cand = PupInfo->begin();cand != PupInfo->end(); ++ cand ) { 
	std::cout << " Pileup Information: bunchXing, nvtx: " << cand->getBunchCrossing() << " " << cand->getPU_NumInteractions() << std::endl;
	//  const int num_PU_vertices_;          // the number of pileup interactions that have been added to the event
	//  std::vector<float> zpositions_;      // the true primary vertex position along the z axis for each added interaction
	//  std::vector<float> sumpT_lowpT_;     // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > low_cut
	//  std::vector<float> sumpT_highpT_;    // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > high_cut
	//  std::vector<int> ntrks_lowpT_;       // the number of tracks originating from each interaction, where track pT > low_cut
	//  std::vector<int> ntrks_highpT_;      // the number of tracks originating from each interaction, where track pT > high_cut
	num_PU_vertices=cand->getPU_NumInteractions();
	PU_BunchCrossing=cand->getBunchCrossing();
	/* for (int j=0;j<num_PU_vertices;j++){ */
/* 	  PU_zpos[j]=(cand->getPU_zpositions()).at(j); */
/* 	  PU_sumpT_lowpT[j]=(cand->getPU_sumpT_lowpT()).at(j); */
/* 	  PU_sumpT_highpT[j]=(cand->getPU_sumpT_highpT()).at(j); */
/* 	  PU_ntrks_lowpT[j]=(cand->getPU_ntrks_lowpT()).at(j); */
/* 	  PU_ntrks_highpT[j]=(cand->getPU_ntrks_highpT()).at(j); */
/* 	} */
      }	
  }


  void fillHLTFired(const edm::Event& iEvent){
    edm::Handle<vector<std::string> > HLTfired_;
    iEvent.getByLabel(HLTInfoFired,HLTfired_);

    vector<string> HLTimported;
    string tmpstring="";

    for (vector<string>::const_iterator cand=HLTfired_->begin(); cand!=HLTfired_->end(); ++cand){
      unsigned int i=cand-HLTfired_->begin();
      HLTimported.push_back(cand->c_str());
      string newstr=HLTimported.at(i) + ":" + tmpstring;
      tmpstring=newstr;
    }

    //cout << "HLTFiredString========>>>>>>>>>>> " << tmpstring.c_str() << endl;
    sprintf(HLTPathsFired,tmpstring.c_str());
    cout << "HLTFiredString========>>>>>>>>>>> " << HLTPathsFired << endl;
  }

  void fillHLT(const edm::Event& iEvent){
    for (unsigned int i=0;i<flagHLTnames.size();i++){
      edm::Handle<bool> HLTflagHandle;
      iEvent.getByLabel(HLTAnalysisinst.c_str(),flagHLTnames.at(i).label().c_str(),HLTflagHandle);
      const bool *HLTflagValue=HLTflagHandle.product();
      cout << "HLTflagValue " << *HLTflagValue<< endl;
      boolHLTpattern[i]=*HLTflagValue;
    }
      
  }

  void fillSkimEarlyData(const edm::Event& iEvent){

    for (unsigned int i=0;i<flagSkimEarlyDatanames.size();i++){
      edm::Handle<bool> SkimEarlyDataflagHandle;
      iEvent.getByLabel(SkimEarlyDataAnalysisinst.c_str(),flagSkimEarlyDatanames.at(i).label().c_str(),SkimEarlyDataflagHandle);
      const bool *SkimEarlyDataflagValue=SkimEarlyDataflagHandle.product();
      cout << "SkimEarlyDataflagValue " << *SkimEarlyDataflagValue<< endl;

      if(i==0){
	boolSkim_highEnergyMuons=*SkimEarlyDataflagValue;
      }
      else if(i==1){
	boolSkim_highEnergyElectrons=*SkimEarlyDataflagValue;
      }
      else if(i==2){
	boolSkim_recoWMNfromPf=*SkimEarlyDataflagValue;
      }
      else if(i==3){
	boolSkim_recoWMNfromTc=*SkimEarlyDataflagValue;
      }
      else if(i==4){
	boolSkim_recoWENfromPf=*SkimEarlyDataflagValue;
      }
      else if(i==5){
	boolSkim_recoWENfromTc=*SkimEarlyDataflagValue;
      }
      else if(i==6){
        boolSkim_diMuonsJPsi=*SkimEarlyDataflagValue;
      }
      else if(i==7){
        boolSkim_diMuonsZ=*SkimEarlyDataflagValue;
      }
      else if(i==8){
	boolSkim_diElectronsZ=*SkimEarlyDataflagValue;
      }
      else if(i==9){
	boolSkim_triLeptonsMuMuMu=*SkimEarlyDataflagValue;
      }
      else if(i==10){
	boolSkim_triLeptonsMuMuEl=*SkimEarlyDataflagValue;
      }
      else if(i==11){
	boolSkim_triLeptonsMuElEl=*SkimEarlyDataflagValue;
      }
      else if(i==12){
	boolSkim_triLeptonsElElEl=*SkimEarlyDataflagValue;
      }
      else if(i==13){
	boolSkim_quadLeptons4Mu=*SkimEarlyDataflagValue;
      }
      else if(i==14){
	boolSkim_quadLeptons2Mu2El=*SkimEarlyDataflagValue;
      }
      else if(i==15){
	boolSkim_quadLeptons4El=*SkimEarlyDataflagValue;
      }
    }
  }


  void triggermatching(const edm::Event& iEvent){
    
     cout << "Start Trigger matching for muon" << endl;

    edm::Handle<trigger::TriggerEvent> handleTriggerEvent;
    iEvent.getByLabel(triggerEvent, handleTriggerEvent );
    const trigger::TriggerObjectCollection & toc(handleTriggerEvent->getObjects());
    size_t nMuHLT =0;

    std::vector<reco::Particle>  HLTMuMatched;
    for ( size_t ia = 0; ia < handleTriggerEvent->sizeFilters(); ++ ia) {
      std::string fullname = handleTriggerEvent->filterTag(ia).encode();
      //std::cout<< "fullname::== " << fullname<< std::endl;
      std::string name;
      size_t p = fullname.find_first_of(':');
      if ( p != std::string::npos) {
	name = fullname.substr(0, p);
      }
      else {
	name = fullname;
      }
      //std::cout<< "name::== " << name<< std::endl;
      if ( &toc !=0 ) {
	const trigger::Keys & k = handleTriggerEvent->filterKeys(ia);
	for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
	  // looking at all the single muon l3 trigger present, for example hltSingleMu15L3Filtered15.....
	  // cout << "name=" << name << endl;
	  if (name == triggerFilter.c_str() ) {
	    HLTMuMatched.push_back(toc[*ki].particle());
	    cout << "Matching " << triggerFilter.c_str()  << endl;
	    nMuHLT++;
          }
	  else {
	    for (unsigned int l=0;l<triggerFilter_asym.size();l++){
	      if (name == triggerFilter_asym.at(l).c_str() ) {
		HLTMuMatched.push_back(toc[*ki].particle());
		cout << "Matching " << triggerFilter_asym.at(l).c_str()<< endl;
		nMuHLT++;
	      }
	    }
          }
	}
      }
    }


    // Based on PAT processing
    edm::Handle<edm::View<reco::Muon> > MuCandidates;
    iEvent.getByLabel(muonTag_, MuCandidates);

    typedef std::vector< TriggerObjectStandAlone > TriggerObjectStandAloneCollection;
    edm::Handle<edm::Association<TriggerObjectStandAloneCollection> > matches;
    iEvent.getByLabel(triggerMatchObject, matches);

    //edm::Handle<TriggerObjectStandAloneMatch> match;
    //iEvent.getByLabel("muonTriggerMatchHLT", match);

    // edm::Handle<TriggerObjectStandAloneCollection> pattrigger_;
    // iEvent.getByLabel("patTrigger", pattrigger_);
    
 
    float maxDeltaR_=0.2;
    float maxDPtRel_=1.0;
    int nMuHLTMatch=0,nMuHLTMatchPat=0,nMuHLTMatch_asym_Pat=0;

    for (edm::View<reco::Muon>::const_iterator iCand = MuCandidates->begin(); iCand != MuCandidates->end(); ++iCand){

      unsigned int i=iCand-MuCandidates->begin();
      cout << "Muon with pt= " << iCand->pt() << ": check trigger matching" << endl;

      if (IsMuMatchedToHLTMu ( *iCand,  HLTMuMatched ,maxDeltaR_, maxDPtRel_ )==true){
	nMuHLTMatch++;
	cout << "Muon HLT Matched"  << endl;
	RECOMU_PT_MuHLTMatch[i]=iCand->pt();
      }

      edm::Ref<edm::View<reco::Muon> > muonRef(MuCandidates,i);
      edm::Ref<std::vector<pat::TriggerObjectStandAlone> > pattriggerref = (*matches)[muonRef];
      unsigned idx=1;
      //bool pathLastFilterAccepted=true;
      if (pattriggerref.isNonnull() && pattriggerref.isAvailable()){
	//for (unsigned int j=0;j<pattriggerref->pathNames(pathLastFilterAccepted).size();j++){
	//  cout << (pattriggerref->pathNames(pathLastFilterAccepted)).at(j).c_str() << endl;
	//}
	//for (unsigned int j=0;j<pattriggerref->filterLabels().size();j++){
	//  cout << (pattriggerref->filterLabels()).at(j).c_str() << endl;
	//}
	//cout << "Muon HLT Matching PAT: ref= " << pattriggerref->hasPathName(std::string("HLT_DoubleMu3_v*"),idx) << " " 
	//   << pattriggerref->hasFilterLabel(std::string("hltDiMuonL3PreFiltered5")) << " " 
	//   << pattriggerref->hasCollection(std::string("hltL3MuonCandidates"))
	//   << endl;	
	//if (pattriggerref->hasPathName(std::string("HLT_DoubleMu7_v*"),idx)==true ||
	//  (pattriggerref->hasPathName(std::string("HLT_DoubleMu5_v*"),idx)==true  && pattriggerref->hasFilterLabel(std::string("hltDiMuonL3PreFiltered5"))==true)
	//  ) {	
	if (pattriggerref->hasCollection(std::string(triggerHLTcollection.c_str()))==true && 
	    ( pattriggerref->hasPathName(std::string("HLT_DoubleMu7_v*"),idx)==true
	      ||
	      (pattriggerref->hasPathName(std::string("HLT_DoubleMu5_v*"),idx)==true  && pattriggerref->hasFilterLabel(std::string(triggerFilter.c_str()))==true))
	    //	      (pattriggerref->hasPathName(std::string("HLT_DoubleMu5_v*"),idx)==true  && pattriggerref->hasFilterLabel(std::string("hltDiMuonL3PreFiltered5"))==true))
	    )
	  {	    
	    cout << "Muon HLT Matched PAT" << endl;
	    nMuHLTMatchPat++;
	    RECOMU_PT_MuHLTMatchPAT[i]=iCand->pt();
	  }
      }
    }
      
    

    cout << "N. Muons HLT Matched= " << nMuHLTMatch << " and PAT= " << nMuHLTMatchPat << " FiredString:" << HLTPathsFired << endl;

    RECO_nMuHLTMatch    = nMuHLTMatch;
    RECO_nMuHLTMatchPAT = nMuHLTMatchPat;

    // Double Mu asymmetric trigger
    edm::Handle<edm::Association<TriggerObjectStandAloneCollection> > matches_asym;
    iEvent.getByLabel(triggerMatchObject_asym, matches_asym);


    for (edm::View<reco::Muon>::const_iterator iCand = MuCandidates->begin(); iCand != MuCandidates->end(); ++iCand){
      
      unsigned int i=iCand-MuCandidates->begin();
      cout << "Muon with pt= " << iCand->pt() << ": check trigger matching" << endl;
      
      edm::Ref<edm::View<reco::Muon> > muonRef(MuCandidates,i);
      edm::Ref<std::vector<pat::TriggerObjectStandAlone> > pattriggerref = (*matches_asym)[muonRef];
      unsigned idx=1;
      
      //bool pathLastFilterAccepted=true;
      if (pattriggerref.isNonnull() && pattriggerref.isAvailable()){
        if (pattriggerref->hasCollection(std::string(triggerHLTcollection.c_str()))==true &&
            ( (pattriggerref->hasPathName(std::string("HLT_Mu13_Mu8_v*"),idx)==true || pattriggerref->hasPathName(std::string("HLT_Mu17_Mu8_v*"),idx)==true ))){
          cout << "Muon HLT Matched asymmetric PAT" << endl;
          nMuHLTMatch_asym_Pat++;
          RECOMU_PT_MuHLTMatch_asym_PAT[i]=iCand->pt();
        }
        else {
          for(unsigned int l=0;l<triggerFilter_asym.size();l++){
            if (pattriggerref->hasCollection(std::string(triggerHLTcollection.c_str()))==true &&
                ( (pattriggerref->hasPathName(std::string("HLT_Mu13_Mu8_v*"),idx)==true || pattriggerref->hasPathName(std::string("HLT_Mu17_Mu8_v*"),idx)==true )
                  ||
                  ((pattriggerref->hasPathName(std::string("HLT_Mu13_Mu8_v*"),idx)==true || pattriggerref->hasPathName(std::string("HLT_Mu17_Mu8_v*"),idx)==true ) && pattriggerref->hasFilterLabel(std::string(triggerFilter_asym.at(l).c_str()))==true))
                )
              {
                cout << "Muon HLT Matched asymmetric PAT" << endl;
                nMuHLTMatch_asym_Pat++;
                RECOMU_PT_MuHLTMatch_asym_PAT[i]=iCand->pt();
              }
          }
        }
      }
    }
    
    
    
    cout << "N. Muons HLT asym Matched= " << nMuHLTMatch << " and asym PAT= " << nMuHLTMatch_asym_Pat << endl;
    
    RECO_nMuHLTMatch    = nMuHLTMatch;
    RECO_nMuHLTMatch_asym_PAT = nMuHLTMatch_asym_Pat;
  




    cout << "Start Trigger matching for electron" << endl;

    int nEleHLTMatchPat=0;
    edm::Handle<edm::View<reco::GsfElectron> > EleCandidates;
    iEvent.getByLabel(electronEgmTag_, EleCandidates);

    typedef std::vector< TriggerObjectStandAlone > TriggerObjectStandAloneCollection;
    edm::Handle<edm::Association<TriggerObjectStandAloneCollection> > matchesele;
    iEvent.getByLabel(triggerMatchObjectEle, matchesele);


    for (edm::View<reco::GsfElectron>::const_iterator iCand = EleCandidates->begin(); iCand != EleCandidates->end(); ++iCand){

      unsigned int i=iCand-EleCandidates->begin();
      cout << "Electron with pt= " << iCand->pt() << ": check trigger matching" << endl;

      edm::Ref<edm::View<reco::GsfElectron> > electronRef(EleCandidates,i);
      edm::Ref<std::vector<pat::TriggerObjectStandAlone> > pattriggerref = (*matchesele)[electronRef];
      unsigned idx=1;
      bool pathLastFilterAccepted=true;
      if (pattriggerref.isNonnull() && pattriggerref.isAvailable()){

	if (
	    ( pattriggerref->hasPathName(std::string("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*"),idx)==true
	      ||
	      pattriggerref->hasPathName(std::string("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*"),idx)==true
	      ||
	      pattriggerref->hasPathName(std::string("HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v*"),idx)==true
	      ))
	  {	    
	    cout << "Electron HLT Matching PAT" << endl;
	    nEleHLTMatchPat++;
	    RECOELE_PT_EleHLTMatchPAT[i]=iCand->pt();
	  }
      }

    }

    cout << "N. Electrons HLT Matched= " << nEleHLTMatchPat << endl;

    RECO_nEleHLTMatchPAT = nEleHLTMatchPat;


  }



  bool IsMuMatchedToHLTMu ( const reco::Muon &mu, std::vector<reco::Particle> HLTMu , double DR, double DPtRel ) {
    size_t dim =  HLTMu.size();
    size_t nPass=0;
    if (dim==0) return false;
    for (size_t k =0; k< dim; k++ ) {
      if (  (deltaR(HLTMu[k], mu) < DR)   && (fabs(HLTMu[k].pt() - mu.pt())/ HLTMu[k].pt()<DPtRel)){     
	nPass++ ;
      }
    }
    return (nPass>0);
  }
  
  void fillmc2e2mu(const edm::Event& iEvent){
    int k=0, l=0;
    for (unsigned int i=0; i<MCcollName.size(); i++) { 
      edm::Handle<reco::CandidateCollection> Candidates;
      iEvent.getByLabel(MCcollName.at(i), Candidates);
      for( reco::CandidateCollection::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetMCValues2e2mu(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const reco::Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetMCValues2e2mu(theParticle,l);
	  }
	}
      }
    }
    
    fillMCCP2e2mu(iEvent);
  }

  
  void fillmc4mu(const edm::Event& iEvent){
    int k=0, l=0;
    for (unsigned int i=0; i<MCcollName.size(); i++) { 
      edm::Handle<reco::CandidateCollection> Candidates;
      iEvent.getByLabel(MCcollName.at(i), Candidates);
      for( reco::CandidateCollection::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetMCValues4mu(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const reco::Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetMCValues4mu(theParticle,l);
	  }
	}
      }
    }
    
    fillMCCP4mu(iEvent);
  }


  void fillmc4e(const edm::Event& iEvent){
    int k=0, l=0;
    for (unsigned int i=0; i<MCcollName.size(); i++) { 
      edm::Handle<reco::CandidateCollection> Candidates;
      iEvent.getByLabel(MCcollName.at(i), Candidates);
      for( reco::CandidateCollection::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetMCValues4e(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const reco::Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetMCValues4e(theParticle,l);
	  }
	}
      }
    }
    
    fillMCCP4e(iEvent);
  }

  
  struct SortCandByDecreasingPt {
    bool operator()( const reco::Candidate &c1, const reco::Candidate &c2) const {
      return c1.pt() > c2.pt();
    }
  };
  
  void fillRECO2e2mu(const edm::Event& iEvent){
    
    // 2e2mu
    leptonscands2e2mu_->clear();
    
    // Fill Higgs, Z and their daughters if any 
    int k=0, l=0;
    for (unsigned int i=0; i<3; i++) {  
      edm::Handle<edm::View<Candidate> > Candidates;
      iEvent.getByLabel(RECOcollNameBest2e2mu.at(i), Candidates);
      for( edm::View<Candidate>::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetRECOValues2e2mu(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetRECOValues2e2mu(theParticle,l);
	    leptonscands2e2mu_->push_back( cand->daughter(j)->clone());
	  }
	}
      }
    }
    
    // Ordering leptons in pT  and filling the RECO_ORDERED_PT vector
    leptonscands2e2mu_->sort(SortCandByDecreasingPt());
    int il=3;
    for( CandidateCollection::const_iterator p = leptonscands2e2mu_->begin();p != leptonscands2e2mu_->end(); ++ p ) {   
      // cout << "Momentum of leptons pt= " << p->p4().pt() << endl;
      RECOBEST_2e2mu_ORDERED_PT[il]=p->p4().pt();
      il++;
    }
  }

  void fillRECO4mu(const edm::Event& iEvent){
    // 4mu
    leptonscands4mu_->clear();
    
    // Fill Higgs, Z and their daughters if any 
    int k=0, l=0;
    for (unsigned int i=0; i<3; i++) {  
      edm::Handle<edm::View<Candidate> > Candidates;
      iEvent.getByLabel(RECOcollNameBest4mu.at(i), Candidates);
      for( edm::View<Candidate>::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetRECOValues4mu(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetRECOValues4mu(theParticle,l);
	    leptonscands4mu_->push_back( cand->daughter(j)->clone());
	  }
	}
      }
    }
    
    // Ordering leptons in pT  and filling the RECO_ORDERED_PT vector
    leptonscands4mu_->sort(SortCandByDecreasingPt());
    int il=3;
    for( CandidateCollection::const_iterator p = leptonscands4mu_->begin();p != leptonscands4mu_->end(); ++ p ) {   
      // cout << "Momentum of leptons pt= " << p->p4().pt() << endl;
      RECOBEST_4mu_ORDERED_PT[il]=p->p4().pt();
      il++;
    }

  }

  void fillRECO4e(const edm::Event& iEvent){
    
    // 4e
    leptonscands4e_->clear();
    
    // Fill Higgs, Z and their daughters if any 
    int k=0, l=0;
    for (unsigned int i=0; i<3; i++) {  
      edm::Handle<edm::View<Candidate> > Candidates;
      iEvent.getByLabel(RECOcollNameBest4e.at(i), Candidates);
      for( edm::View<Candidate>::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetRECOValues4e(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetRECOValues4e(theParticle,l);
	    leptonscands4e_->push_back( cand->daughter(j)->clone());
	  }
	}
      }
    }
    
    // Ordering leptons in pT  and filling the RECO_ORDERED_PT vector
    leptonscands4e_->sort(SortCandByDecreasingPt());
    int il=3;
    for( CandidateCollection::const_iterator p = leptonscands4e_->begin();p != leptonscands4e_->end(); ++ p ) {   
      // cout << "Momentum of leptons pt= " << p->p4().pt() << endl;
      RECOBEST_4e_ORDERED_PT[il]=p->p4().pt();
      il++;
    }

  }
  
  void fillAdditionalRECO(const edm::Event& iEvent){
   
    // di-leptons OS
    leptonscands_Z0->clear();
    leptonscands_Z1->clear();
    for (unsigned int i=0; i<RECOcollNameZ.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZ;
      iEvent.getByLabel(RECOcollNameZ.at(i), CandidatesZ); 
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZ->begin();cand != CandidatesZ->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_ZMM_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_ZMM_PT[0][kk]=cand->p4().pt();
	if (i==0) RECO_ZMM_ETA[0][kk]=cand->p4().eta();
	if (i==0) RECO_ZMM_PHI[0][kk]=cand->p4().phi();
	if (i==1) RECO_ZEE_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_ZEE_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_ZEE_ETA[0][kk]=cand->p4().eta();
	if (i==1) RECO_ZEE_PHI[0][kk]=cand->p4().phi();
	cout << "di-lepton candidate of type=" << RECOcollNameZ.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEE_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEE_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEE_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Z1->push_back( cand->daughter(j)->clone());
	  }
	}
	kk++;
      }
    }

    // di-leptons SS and cross-leptons
    leptonscands_Zss0->clear();
    leptonscands_Zss1->clear();
    leptonscands_Zcross->clear();

    for (unsigned int i=0; i<RECOcollNameZss.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesZss;
      iEvent.getByLabel(RECOcollNameZss.at(i), CandidatesZss); 
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesZss->begin();cand != CandidatesZss->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0)  RECO_ZMMss_MASS[kk]=cand->p4().mass();
	if (i==0)  RECO_ZMMss_PT[0][kk]=cand->p4().pt();
	if (i==0)  RECO_ZMMss_ETA[0][kk]=cand->p4().eta();
	if (i==0)  RECO_ZMMss_PHI[0][kk]=cand->p4().phi();

	if (i==1)  RECO_ZEEss_MASS[kk]=cand->p4().mass();
	if (i==1)  RECO_ZEEss_PT[0][kk]=cand->p4().pt();
	if (i==1)  RECO_ZEEss_ETA[0][kk]=cand->p4().eta();
	if (i==1)  RECO_ZEEss_PHI[0][kk]=cand->p4().phi();

	if (i==2)  RECO_ZEM_MASS[kk]=cand->p4().mass();
	if (i==2)  RECO_ZEM_PT[0][kk]=cand->p4().pt();
	if (i==2)  RECO_ZEM_ETA[0][kk]=cand->p4().eta();
	if (i==2)  RECO_ZEM_PHI[0][kk]=cand->p4().phi();

	cout << "di-lepton candidate of type=" << RECOcollNameZss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl; 
	  if (i==0) {
	    RECO_ZMMss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZMMss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZMMss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_ZEEss_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEEss_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEEss_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zss1->push_back( cand->daughter(j)->clone());
	  }
	  if (i==2) {
	    RECO_ZEM_PT[j+1][kk]=cand->daughter(j)->p4().pt();
	    RECO_ZEM_ETA[j+1][kk]=cand->daughter(j)->p4().eta();
	    RECO_ZEM_PHI[j+1][kk]=cand->daughter(j)->p4().phi();
	    leptonscands_Zcross->push_back( cand->daughter(j)->clone());	  
	  }
	}
	kk++;
      }
    }


    // di-leptons ALL
    leptonscands_DiLep->clear();

    edm::Handle<edm::View<Candidate> > CandidatesDiLep;
    iEvent.getByLabel(RECOcollNameDiLep, CandidatesDiLep); 
    int kkk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesDiLep->begin();cand != CandidatesDiLep->end(); ++ cand ) { 
      if (kkk>49) continue;
      RECO_DiLep_MASS[kkk]=cand->p4().mass();
      RECO_DiLep_PT[0][kkk]=cand->p4().pt();
      RECO_DiLep_ETA[0][kkk]=cand->p4().eta();
      RECO_DiLep_PHI[0][kkk]=cand->p4().phi();
      
      cout << "di-lepton candidate of type=" << RECOcollNameDiLep.label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	cout << "Daugthter with pt and charge=" << cand->daughter(j)->p4().pt() << " " << cand->daughter(j)->charge() << endl;
	RECO_DiLep_PT[j+1][kkk]=cand->daughter(j)->p4().pt();
	RECO_DiLep_ETA[j+1][kkk]=cand->daughter(j)->p4().eta();
	RECO_DiLep_PHI[j+1][kkk]=cand->daughter(j)->p4().phi();
	leptonscands_DiLep->push_back( cand->daughter(j)->clone());
      }
      kkk++;
    }



   /*  // MuMuMuMu */
/*     int i=1;  */
/*     leptonscands_MMMM->clear(); */
/*     edm::Handle<edm::View<Candidate> > CandidatesMMMM; */
/*     iEvent.getByLabel(RECOcollNameMMMM.at(0), CandidatesMMMM); */
/*     int kk=0; */
/*     for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) {  */
/*       if (kk>49) continue; */
/*       RECO_MMMM_MASS[i-1][kk]=cand->p4().mass(); */
/*       RECO_MMMM_PT[i-1][kk]=cand->p4().pt(); */
/*       RECO_MMMM_ETA[i-1][kk]=cand->p4().eta(); */
/*       RECO_MMMM_PHI[i-1][kk]=cand->p4().phi(); */
/*       int l=0; */
/*       //cout << "index" << i-1 << endl; */
/*       for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) { */
/* 	RECO_MMMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass(); */
/* 	RECO_MMMM_PT[i+j][kk]=cand->daughter(j)->p4().pt(); */
/* 	RECO_MMMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta(); */
/* 	RECO_MMMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi(); */
/* 	//cout << "index" << i+j <<endl; */
/* 	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) { */
/* 	  RECO_MMMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass(); */
/* 	  RECO_MMMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt(); */
/* 	  RECO_MMMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta(); */
/* 	  RECO_MMMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi(); */
/* 	  leptonscands_MMMM->push_back( cand->daughter(j)->daughter(k)->clone()); */
/* 	  //cout << "index" << i+j+k+l+2 <<endl; */
/* 	}  */
/* 	l++; */
/*       } */
/*       kk++; */
/*     } */

    
/*     // EEEE */
/*     i=1; */
/*     leptonscands_EEEE->clear(); */
/*     edm::Handle<edm::View<Candidate> > CandidatesEEEE; */
/*     iEvent.getByLabel(RECOcollNameEEEE.at(0), CandidatesEEEE); */
/*     kk=0; */
/*     for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) { */
/*       if (kk>49) continue; */
/*       RECO_EEEE_MASS[i-1][kk]=cand->p4().mass(); */
/*       RECO_EEEE_PT[i-1][kk]=cand->p4().pt(); */
/*       RECO_EEEE_ETA[i-1][kk]=cand->p4().eta(); */
/*       RECO_EEEE_PHI[i-1][kk]=cand->p4().phi(); */
/*       int l=0; */
/*       //cout << "index" << i-1 << endl; */
/*       for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) { */
/* 	RECO_EEEE_MASS[i+j][kk]=cand->daughter(j)->p4().mass(); */
/* 	RECO_EEEE_PT[i+j][kk]=cand->daughter(j)->p4().pt(); */
/* 	RECO_EEEE_ETA[i+j][kk]=cand->daughter(j)->p4().eta(); */
/* 	RECO_EEEE_PHI[i+j][kk]=cand->daughter(j)->p4().phi(); */
/* 	//cout << "index" << i+j <<endl; */
/* 	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) { */
/* 	  RECO_EEEE_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass(); */
/* 	  RECO_EEEE_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt(); */
/* 	  RECO_EEEE_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta(); */
/* 	  RECO_EEEE_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi(); */
/* 	  leptonscands_EEEE->push_back( cand->daughter(j)->daughter(k)->clone()); */
/* 	  //cout << "index" << i+j+k+l+2 <<endl; */
/* 	} */
/* 	l++; */
/*       } */
/*       kk++; */
/*     } */


/*     // EEMM */
/*     i=1; */
/*     leptonscands_EEMM->clear(); */
/*     edm::Handle<edm::View<Candidate> > CandidatesEEMM; */
/*     iEvent.getByLabel(RECOcollNameEEMM.at(0), CandidatesEEMM); */
/*     kk=0; */
/*     for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) { */
/*       if (kk>49) continue; */
/*       RECO_EEMM_MASS[i-1][kk]=cand->p4().mass(); */
/*       RECO_EEMM_PT[i-1][kk]=cand->p4().pt(); */
/*       RECO_EEMM_ETA[i-1][kk]=cand->p4().eta(); */
/*       RECO_EEMM_PHI[i-1][kk]=cand->p4().phi(); */
/*       int l=0; */
/*       //cout << "index" << i-1 << endl; */
/*       for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) { */
/* 	RECO_EEMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass(); */
/* 	RECO_EEMM_PT[i+j][kk]=cand->daughter(j)->p4().pt(); */
/* 	RECO_EEMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta(); */
/* 	RECO_EEMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi(); */
/* 	//cout << "index" << i+j <<endl; */
/* 	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) { */
/* 	  RECO_EEMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass(); */
/* 	  RECO_EEMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt(); */
/* 	  RECO_EEMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta(); */
/* 	  RECO_EEMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi(); */

/* 	  leptonscands_EEMM->push_back( cand->daughter(j)->daughter(k)->clone()); */
/* 	  //cout << "index" << i+j+k+l+2 <<endl; */
/* 	} */
/* 	l++; */
/*       } */
/*       kk++; */
/*     } */
    
   
    // tri-leptons
    leptonscands_LLL0->clear();
    leptonscands_LLL1->clear();
    leptonscands_LLL2->clear();
    leptonscands_LLL3->clear();

    for (unsigned int i=0; i<RECOcollNameLLL.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesLLL;
      iEvent.getByLabel(RECOcollNameLLL.at(i), CandidatesLLL); 
      int k=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLL->begin();cand != CandidatesLLL->end(); ++ cand ) { 
	if (k>49) continue;
	if (i==0) RECO_LLL0_MASS[k]=cand->p4().mass();
	if (i==0) RECO_LLL0_PT[0][k]=cand->p4().pt();
	if (i==1) RECO_LLL1_MASS[k]=cand->p4().mass();
	if (i==1) RECO_LLL1_PT[0][k]=cand->p4().pt();
	if (i==2) RECO_LLL2_MASS[k]=cand->p4().mass();
	if (i==2) RECO_LLL2_PT[0][k]=cand->p4().pt();
	if (i==3) RECO_LLL3_MASS[k]=cand->p4().mass();
	if (i==4) RECO_LLL3_PT[0][k]=cand->p4().pt();
	cout << "Tri-lepton candidate of type=" << RECOcollNameLLL.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  if (i==0) {
	    RECO_LLL0_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL0->push_back( cand->daughter(j)->clone());
	  }
	  if (i==1) {
	    RECO_LLL1_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL1->push_back( cand->daughter(j)->clone());
	  }
	  if (i==2) {
	    RECO_LLL2_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL2->push_back( cand->daughter(j)->clone()); 
	  }
	  if (i==3) {
	    RECO_LLL3_PT[j+1][k]=cand->daughter(j)->p4().pt();
	    leptonscands_LLL3->push_back( cand->daughter(j)->clone());
	  }
	}
	k++;
      }
    }
    
    // 4-leptons SS
    leptonscands_LLLLss0->clear();
    leptonscands_LLLLss1->clear();
    leptonscands_LLLLss2->clear();
    
    for (unsigned int i=0; i<RECOcollNameLLLLss.size(); i++) {  
      edm::Handle<edm::View<Candidate> > CandidatesLLLLss;
      iEvent.getByLabel(RECOcollNameLLLLss.at(i), CandidatesLLLLss);
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLLLss->begin();cand != CandidatesLLLLss->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_LLLL0ss_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_LLLL0ss_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_LLLL1ss_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_LLLL1ss_PT[0][kk]=cand->p4().pt();
	if (i==2) RECO_LLLL2ss_MASS[kk]=cand->p4().mass();
	if (i==2) RECO_LLLL2ss_PT[0][kk]=cand->p4().pt();
	cout << "4-lepton ss candidate of type=" << RECOcollNameLLLLss.at(i).label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	    if (i==0) {
	      if (j==0) RECO_LLLL0ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL0ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss0->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==1) {
	      if (j==0) RECO_LLLL1ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL1ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss1->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==2) {
	      if (j==0) RECO_LLLL2ss_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      if (j==1) RECO_LLLL2ss_PT[j+k+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLLss2->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	  }
	}
	kk++;
      }      
    }

    // 4-leptons 3l+l
    leptonscands_LLLl0->clear();
    leptonscands_LLLl1->clear(); 
    
    for (unsigned int i=0; i<RECOcollNameLLLl.size(); i++){ 
      edm::Handle<edm::View<Candidate> > CandidatesLLLl;  
      iEvent.getByLabel(RECOcollNameLLLl.at(i), CandidatesLLLl);  
      int kk=0;
      for( edm::View<Candidate>::const_iterator cand = CandidatesLLLl->begin();cand != CandidatesLLLl->end(); ++ cand ) { 
	if (kk>49) continue;
	if (i==0) RECO_LLLl0_MASS[kk]=cand->p4().mass();
	if (i==0) RECO_LLLl0_PT[0][kk]=cand->p4().pt();
	if (i==1) RECO_LLLl1_MASS[kk]=cand->p4().mass();
	if (i==1) RECO_LLLl1_PT[0][kk]=cand->p4().pt();
	cout << "4-lepton from 3l+l candidate of type=" << RECOcollNameLLLl.at(i).label()
	     << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
	for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	  //cout << "j= " << j << endl;
	  if ( i==0 && j==1) {
	    RECO_LLLl0_PT[j+3][kk]=cand->daughter(j)->p4().pt();
	    leptonscands_LLLl0->push_back( cand->daughter(j)->clone());
	  }
	  if ( i==1 && j==1) {
	    RECO_LLLl1_PT[j+3][kk]=cand->daughter(j)->p4().pt(); 
	    leptonscands_LLLl1->push_back( cand->daughter(j)->clone());
	  }
	  for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) { 
	    //cout << "k= " << k << endl;
	    if (i==0) {
	      if (j==0) RECO_LLLl0_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLl0->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    if (i==1) {
	      if (j==0) RECO_LLLl1_PT[j+k+1][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	      leptonscands_LLLl1->push_back( cand->daughter(j)->daughter(k)->clone());
	    }
	    
	  } 
	} 
	kk++;
      }
    }


    // LLLL merged (no flavour, no charge)

    int i=1;
    leptonscands_LLLL->clear();
    edm::Handle<edm::View<Candidate> > CandidatesLLLL;
    iEvent.getByLabel(RECOcollNameLLLL, CandidatesLLLL);
    int kk=0;    
    for( edm::View<Candidate>::const_iterator cand = CandidatesLLLL->begin();cand != CandidatesLLLL->end(); ++ cand ) {
      if (kk>99) continue;
      cout << "4lepton (any flavour and charge) candidate of type=" << RECOcollNameLLLL.label() << " of mass=" << cand->p4().mass() << " and pt=" << cand->p4().pt() <<endl;
      
      RECO_LLLL_MASS[0][kk]=cand->p4().mass();
      RECO_LLLL_PT[0][kk]=cand->p4().pt();
      RECO_LLLL_ETA[0][kk]=cand->p4().eta();
      RECO_LLLL_PHI[0][kk]=cand->p4().phi();
      int l=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_LLLL_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_LLLL_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_LLLL_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_LLLL_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_LLLL_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_LLLL_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_LLLL_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_LLLL_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  
	  leptonscands_LLLL->push_back( cand->daughter(j)->daughter(k)->clone());
	  //cout << "index" << i+j+k+l+2 <<endl;
	}
	l++;
      }
      kk++;
    }
    
  }

  // RF
  void fillRECOrf2e2mu(const edm::Event& iEvent){    
    leptonscands2e2murf_->clear();
    
    int k=0, l=0;
    for (unsigned int i=0; i<3; i++) {  
      edm::Handle<CandidateCollection> Candidates;
      iEvent.getByLabel(RECOcollNameBestRestFrame2e2mu.at(i), Candidates);
      for( CandidateCollection::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetRECORF2e2muValues(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetRECORF2e2muValues(theParticle,l);
	    leptonscands2e2murf_->push_back( cand->daughter(j)->clone());
	  }
	}
      }
    }     
  }

  void fillRECOrf4mu(const edm::Event& iEvent){    
    leptonscands4murf_->clear();
    
    int k=0, l=0;
    for (unsigned int i=0; i<3; i++) {  
      edm::Handle<CandidateCollection> Candidates;
      iEvent.getByLabel(RECOcollNameBestRestFrame4mu.at(i), Candidates);
      for( CandidateCollection::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetRECORF4muValues(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetRECORF4muValues(theParticle,l);
	    leptonscands4murf_->push_back( cand->daughter(j)->clone());
	  }
	}
      }
    }     
  }


  void fillRECOrf4e(const edm::Event& iEvent){    
    leptonscands4erf_->clear();
    
    int k=0, l=0;
    for (unsigned int i=0; i<3; i++) {  
      edm::Handle<CandidateCollection> Candidates;
      iEvent.getByLabel(RECOcollNameBestRestFrame4e.at(i), Candidates);
      for( CandidateCollection::const_iterator cand = Candidates->begin();cand != Candidates->end(); ++ cand ) { 
	const reco::Candidate& theParticle = (*cand);
	SetRECORF4eValues(theParticle,i);
	if (i>0) {
	  k=i-1;
	  for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	    l=i+j+k+2;
	    const Candidate* daughter;
	    daughter = cand->daughter(j);
	    const reco::Candidate& theParticle = (*daughter);
	    SetRECORF4eValues(theParticle,l);
	    leptonscands4erf_->push_back( cand->daughter(j)->clone());
	  }
	}
      }
    }     
  }


  void SetMCValues2e2mu(const reco::Candidate& cand, int nMC){
    MC_2e2mu_E[nMC]     = cand.p4().energy();
    MC_2e2mu_PT[nMC]    = cand.p4().pt();
    MC_2e2mu_ETA[nMC]   = cand.p4().eta();
    MC_2e2mu_THETA[nMC] = cand.p4().theta();
    MC_2e2mu_PHI[nMC]   = cand.p4().phi();
    MC_2e2mu_MASS[nMC]  = cand.p4().mass();
    MC_2e2mu_PDGID[nMC] = cand.pdgId();
  }
  
  void SetMCValues4mu(const reco::Candidate& cand, int nMC){
    MC_4mu_E[nMC]     = cand.p4().energy();
    MC_4mu_PT[nMC]    = cand.p4().pt();
    MC_4mu_ETA[nMC]   = cand.p4().eta();
    MC_4mu_THETA[nMC] = cand.p4().theta();
    MC_4mu_PHI[nMC]   = cand.p4().phi();
    MC_4mu_MASS[nMC]  = cand.p4().mass();
    MC_4mu_PDGID[nMC] = cand.pdgId();
  }


  void SetMCValues4e(const reco::Candidate& cand, int nMC){
    MC_4e_E[nMC]     = cand.p4().energy();
    MC_4e_PT[nMC]    = cand.p4().pt();
    MC_4e_ETA[nMC]   = cand.p4().eta();
    MC_4e_THETA[nMC] = cand.p4().theta();
    MC_4e_PHI[nMC]   = cand.p4().phi();
    MC_4e_MASS[nMC]  = cand.p4().mass();
    MC_4e_PDGID[nMC] = cand.pdgId();
  }


  void SetMCRF2e2muValues(const reco::Candidate& cand, int nMC){  
    MCRF_2e2mu_E[nMC]    = cand.p4().energy();
    MCRF_2e2mu_PT[nMC]   = cand.p4().pt();
    MCRF_2e2mu_ETA[nMC]  = cand.p4().eta();
    MCRF_2e2mu_THETA[nMC]= cand.p4().theta();
    MCRF_2e2mu_PHI[nMC]  = cand.p4().phi();
    MCRF_2e2mu_MASS[nMC] = cand.p4().mass();
    MCRF_2e2mu_PDGID[nMC] = cand.pdgId();
  }

  void SetMCRF4muValues(const reco::Candidate& cand, int nMC){  
    MCRF_4mu_E[nMC]    = cand.p4().energy();
    MCRF_4mu_PT[nMC]   = cand.p4().pt();
    MCRF_4mu_ETA[nMC]  = cand.p4().eta();
    MCRF_4mu_THETA[nMC]= cand.p4().theta();
    MCRF_4mu_PHI[nMC]  = cand.p4().phi();
    MCRF_4mu_MASS[nMC] = cand.p4().mass();
    MCRF_4mu_PDGID[nMC] = cand.pdgId();
  }

  void SetMCRF4eValues(const reco::Candidate& cand, int nMC){  
    MCRF_4e_E[nMC]    = cand.p4().energy();
    MCRF_4e_PT[nMC]   = cand.p4().pt();
    MCRF_4e_ETA[nMC]  = cand.p4().eta();
    MCRF_4e_THETA[nMC]= cand.p4().theta();
    MCRF_4e_PHI[nMC]  = cand.p4().phi();
    MCRF_4e_MASS[nMC] = cand.p4().mass();
    MCRF_4e_PDGID[nMC] = cand.pdgId();
  }

  
  void SetRECOValues2e2mu(const reco::Candidate& cand, int nRECO){
    RECOBEST_2e2mu_E[nRECO]     = cand.p4().energy();
    RECOBEST_2e2mu_PT[nRECO]    = cand.p4().pt();
    RECOBEST_2e2mu_ETA[nRECO]   = cand.p4().eta();
    RECOBEST_2e2mu_THETA[nRECO] = cand.p4().theta();
    RECOBEST_2e2mu_PHI[nRECO]   = cand.p4().phi();
    RECOBEST_2e2mu_MASS[nRECO]  = cand.p4().mass();
    if(cand.isGlobalMuon())     RECOBEST_2e2mu_isGlobalMu[nRECO]=true;
    if(cand.isStandAloneMuon()) RECOBEST_2e2mu_isStandAloneMu[nRECO]=true;
    if(cand.isTrackerMuon())    RECOBEST_2e2mu_isTrackerMu[nRECO]=true;
    if(cand.isCaloMuon())       RECOBEST_2e2mu_isCaloMu[nRECO]=true;
    if(cand.isElectron())       RECOBEST_2e2mu_isElectron[nRECO]=true;
    std::cout << "RECOBEST2e2mu_PT in rootple= " << RECOBEST_2e2mu_PT[nRECO] << " nRECO=" << nRECO << std::endl;
  }

  void SetRECOValues4mu(const reco::Candidate& cand, int nRECO){
    RECOBEST_4mu_E[nRECO]     = cand.p4().energy();
    RECOBEST_4mu_PT[nRECO]    = cand.p4().pt();
    RECOBEST_4mu_ETA[nRECO]   = cand.p4().eta();
    RECOBEST_4mu_THETA[nRECO] = cand.p4().theta();
    RECOBEST_4mu_PHI[nRECO]   = cand.p4().phi();
    RECOBEST_4mu_MASS[nRECO]  = cand.p4().mass();
    if(cand.isGlobalMuon())     RECOBEST_4mu_isGlobalMu[nRECO]=true;
    if(cand.isStandAloneMuon()) RECOBEST_4mu_isStandAloneMu[nRECO]=true;
    if(cand.isTrackerMuon())    RECOBEST_4mu_isTrackerMu[nRECO]=true;
    if(cand.isCaloMuon())       RECOBEST_4mu_isCaloMu[nRECO]=true;
    std::cout << "RECOBEST4mu_PT in rootple= " << RECOBEST_4mu_PT[nRECO] << " nRECO=" << nRECO << std::endl;
  }

  void SetRECOValues4e(const reco::Candidate& cand, int nRECO){
    RECOBEST_4e_E[nRECO]     = cand.p4().energy();
    RECOBEST_4e_PT[nRECO]    = cand.p4().pt();
    RECOBEST_4e_ETA[nRECO]   = cand.p4().eta();
    RECOBEST_4e_THETA[nRECO] = cand.p4().theta();
    RECOBEST_4e_PHI[nRECO]   = cand.p4().phi();
    RECOBEST_4e_MASS[nRECO]  = cand.p4().mass();
    if(cand.isGlobalMuon())     RECOBEST_4e_isGlobalMu[nRECO]=true;
    if(cand.isStandAloneMuon()) RECOBEST_4e_isStandAloneMu[nRECO]=true;
    if(cand.isTrackerMuon())    RECOBEST_4e_isTrackerMu[nRECO]=true;
    if(cand.isCaloMuon())       RECOBEST_4e_isCaloMu[nRECO]=true;
    if(cand.isElectron())       RECOBEST_4e_isElectron[nRECO]=true;
    std::cout << "RECOBEST4e_PT in rootple= " << RECOBEST_4e_PT[nRECO] << " nRECO=" << nRECO << std::endl;
  }
  

  void SetRECORF2e2muValues(const reco::Candidate& cand, int nRECO){
    RECORFBEST_2e2mu_E[nRECO]     = cand.p4().energy();
    RECORFBEST_2e2mu_PT[nRECO]    = cand.p4().pt();
    RECORFBEST_2e2mu_ETA[nRECO]   = cand.p4().eta();
    RECORFBEST_2e2mu_THETA[nRECO] = cand.p4().theta();
    RECORFBEST_2e2mu_PHI[nRECO]   = cand.p4().phi();
    RECORFBEST_2e2mu_MASS[nRECO]  = cand.p4().mass();
    std::cout << "RECORFBEST_2e2mu_PT in rootple= " << RECORFBEST_2e2mu_PT[nRECO] << " nRECO=" << nRECO << std::endl;
  }
  
  void SetRECORF4muValues(const reco::Candidate& cand, int nRECO){
    RECORFBEST_4mu_E[nRECO]     = cand.p4().energy();
    RECORFBEST_4mu_PT[nRECO]    = cand.p4().pt();
    RECORFBEST_4mu_ETA[nRECO]   = cand.p4().eta();
    RECORFBEST_4mu_THETA[nRECO] = cand.p4().theta();
    RECORFBEST_4mu_PHI[nRECO]   = cand.p4().phi();
    RECORFBEST_4mu_MASS[nRECO]  = cand.p4().mass();
    std::cout << "RECORFBEST_4mu_PT in rootple= " << RECORFBEST_4mu_PT[nRECO] << " nRECO=" << nRECO << std::endl;
  }

  void SetRECORF4eValues(const reco::Candidate& cand, int nRECO){
    RECORFBEST_4e_E[nRECO]     = cand.p4().energy();
    RECORFBEST_4e_PT[nRECO]    = cand.p4().pt();
    RECORFBEST_4e_ETA[nRECO]   = cand.p4().eta();
    RECORFBEST_4e_THETA[nRECO] = cand.p4().theta();
    RECORFBEST_4e_PHI[nRECO]   = cand.p4().phi();
    RECORFBEST_4e_MASS[nRECO]  = cand.p4().mass();
    std::cout << "RECORFBEST_4e_PT in rootple= " << RECORFBEST_4e_PT[nRECO] << " nRECO=" << nRECO << std::endl;
  }

  
  bool match(double mass, double pt, int charge, const reco::CandidateCollection *c1Coll){
    
    bool found=false;
    for( reco::CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++ pp ) {
      if ((fabs(pp->p4().mass()-mass)  <0.001 ) &&
	  (fabs(pp->p4().pt()  -pt)    <0.001 ) &&
	  (fabs(pp->charge()   -charge)<0.001 )  ){
	found=true;
	//std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
      }
    }
    return found;
  }
  
  bool matchParticle(double mass, double pt, int charge, const reco::Candidate *c1){
    
    bool found=false;
    if ((fabs(c1->p4().mass()-mass)  <0.001 ) &&
	(fabs(c1->p4().pt()  -pt)    <0.001 ) &&
	(fabs(c1->charge()   -charge)<0.001 )  ){
      found=true;
      //std::cout << "Found lepton in the higgs-like leptons collection" << std::endl;
    }    
    return found;
  }
  
  void fillRECO2l2tau(const edm::Event& iEvent){
    
    // MMMT
    int i=1; 
    leptonscands_MMMT->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMMT;
    iEvent.getByLabel(RECOcollNameMMMT, CandidatesMMMT);
    int kk=0;
    
    //cout<<"MMMTcands=================="<<CandidatesMMMT->size()<<endl;

    for( edm::View<Candidate>::const_iterator cand = CandidatesMMMT->begin();cand != CandidatesMMMT->end(); ++ cand ) { 
      //cout << "Found a MuMuMuTau candidate "<< endl;
      if (kk>199) continue;
      RECO_MMMT_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMMT_PT[i-1][kk]=cand->p4().pt();
      RECO_MMMT_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMMT_PHI[i-1][kk]=cand->p4().phi();
      RECO_MMMT_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index=" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMMT_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMMT_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMMT_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMMT_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_MMMT_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index1=" << i+j <<endl;	
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMMT_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMMT_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMMT_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMMT_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_MMMT_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index2=" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_MMMT_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_MMMT_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_MMMT_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_MMMT_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_MMMT_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_MMMT->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index3=" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }


    // EEET
    i=1; 
    leptonscands_EEET->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEET;
    iEvent.getByLabel(RECOcollNameEEET, CandidatesEEET);

    //cout<<"EEETcands=================="<<CandidatesEEET->size()<<endl;
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEET->begin();cand != CandidatesEEET->end(); ++ cand ) {
      //cout<<"mass of EEET pair====="<<cand->p4().mass()<<endl; 
      //cout << "Found a EEETau candidate "<< endl;
      if (kk>99) continue;
      RECO_EEET_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEET_PT[i-1][kk]=cand->p4().pt();
      RECO_EEET_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEET_PHI[i-1][kk]=cand->p4().phi();
      RECO_EEET_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      //cout<< " mass of RECO_EEET_MASS[0]====="<<RECO_EEET_MASS[i-1][kk]<<"\t mass of cand[0]"<<cand->p4().mass()<<endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEET_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEET_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEET_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEET_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_EEET_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	//cout<< " mass of RECO_EEET_MASS[1]====="<<RECO_EEET_MASS[i+j][kk]<<"\t mass of cand[1]"<<cand->daughter(j)->p4().mass()<<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEET_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEET_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEET_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEET_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_EEET_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  //cout<< " mass of RECO_EEET_MASS[2]====="<<RECO_EEET_MASS[i+j+k+l+2][kk]<<"\t mass of cand[2]"<<cand->daughter(j)->daughter(k)->p4().mass()<<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_EEET_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_EEET_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_EEET_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_EEET_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_EEET_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_EEET->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	    //cout<< " mass of RECO_EEET_MASS[3]====="<<RECO_EEET_MASS[i+j+k+l+m+n+4][kk]<<"\t mass of cand[3]"<<cand->daughter(j)->daughter(k)->daughter(m)->p4().mass()<<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    // MMTT
    i=1; 
    leptonscands_MMTT->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMTT;
    iEvent.getByLabel(RECOcollNameMMTT, CandidatesMMTT);

    //cout<<"MMTTcands=================="<<CandidatesMMTT->size()<<endl;

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMTT->begin();cand != CandidatesMMTT->end(); ++ cand ) { 
      //cout << "Found a MuMuTauTau candidate "<< endl;
      if (kk>199) continue;
      RECO_MMTT_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMTT_PT[i-1][kk]=cand->p4().pt();
      RECO_MMTT_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMTT_PHI[i-1][kk]=cand->p4().phi();
      RECO_MMTT_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMTT_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMTT_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMTT_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMTT_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_MMTT_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMTT_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMTT_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMTT_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMTT_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_MMTT_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_MMTT_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_MMTT_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_MMTT_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_MMTT_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_MMTT_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_MMTT->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    // EETauTau
    i=1; 
    leptonscands_EETT->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEETT;
    iEvent.getByLabel(RECOcollNameEETT, CandidatesEETT);

    //cout<<"EETTcands=================="<<CandidatesEETT->size()<<endl;

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEETT->begin();cand != CandidatesEETT->end(); ++ cand ) { 
      //cout << "Found a EETauTau candidate "<< endl;
      if (kk>199) continue;
      RECO_EETT_MASS[i-1][kk]=cand->p4().mass();
      RECO_EETT_PT[i-1][kk]=cand->p4().pt();
      RECO_EETT_ETA[i-1][kk]=cand->p4().eta();
      RECO_EETT_PHI[i-1][kk]=cand->p4().phi();
      RECO_EETT_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EETT_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EETT_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EETT_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EETT_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_EETT_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EETT_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EETT_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EETT_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EETT_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_EETT_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_EETT_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_EETT_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_EETT_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_EETT_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_EETT_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_EETT->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    // MMET
    i=1; 
    leptonscands_MMET->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMET;
    iEvent.getByLabel(RECOcollNameMMET, CandidatesMMET);

    //cout<<"MMETcands=================="<<CandidatesMMET->size()<<endl;

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMET->begin();cand != CandidatesMMET->end(); ++ cand ) { 
      //cout << "Found a MuMuETau candidate "<< endl;
      if (kk>199) continue;
      RECO_MMET_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMET_PT[i-1][kk]=cand->p4().pt();
      RECO_MMET_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMET_PHI[i-1][kk]=cand->p4().phi();
      RECO_MMET_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMET_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMET_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMET_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMET_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_MMET_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMET_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMET_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMET_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMET_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_MMET_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_MMET_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_MMET_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_MMET_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_MMET_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_MMET_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_MMET->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    // EEMTau
    i=1; 
    leptonscands_EEMT->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEMT;
    iEvent.getByLabel(RECOcollNameEEMT, CandidatesEEMT);

    //cout<<"EEMTcands=================="<<CandidatesEEMT->size()<<endl;

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEMT->begin();cand != CandidatesEEMT->end(); ++ cand ) { 
      //cout << "Found a EEMTau candidate "<< endl;
      if (kk>199) continue;
      RECO_EEMT_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEMT_PT[i-1][kk]=cand->p4().pt();
      RECO_EEMT_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEMT_PHI[i-1][kk]=cand->p4().phi();
      RECO_EEMT_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEMT_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEMT_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEMT_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEMT_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_EEMT_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEMT_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEMT_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEMT_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEMT_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_EEMT_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_EEMT_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_EEMT_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_EEMT_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_EEMT_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_EEMT_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_EEMT->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    // MMME 
    i=1; 
    leptonscands_MMME->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMME;
    iEvent.getByLabel(RECOcollNameMMME, CandidatesMMME);

    //cout<<"MMMEcands=================="<<CandidatesMMME->size()<<endl;

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMME->begin();cand != CandidatesMMME->end(); ++ cand ) { 
      //cout << "Found a MuMuMuE candidate "<< endl;
      if (kk>199) continue;
      RECO_MMME_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMME_PT[i-1][kk]=cand->p4().pt();
      RECO_MMME_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMME_PHI[i-1][kk]=cand->p4().phi();
      RECO_MMME_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMME_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMME_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMME_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMME_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_MMME_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMME_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMME_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMME_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMME_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_MMME_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_MMME_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_MMME_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_MMME_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_MMME_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_MMME_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_MMME->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }
    
    //EEEM
    i=1; 
    leptonscands_EEEM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEEM;
    iEvent.getByLabel(RECOcollNameEEEM, CandidatesEEEM);

    //cout<<"EEEMcands=================="<<CandidatesEEEM->size()<<endl;

    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEEM->begin();cand != CandidatesEEEM->end(); ++ cand ) { 
      //cout << "Found a EEEMu candidate "<< endl;
      if (kk>199) continue;
      RECO_EEEM_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEEM_PT[i-1][kk]=cand->p4().pt();
      RECO_EEEM_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEEM_PHI[i-1][kk]=cand->p4().phi();
      RECO_EEEM_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEEM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEEM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEEM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEEM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_EEEM_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEEM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEEM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEEM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEEM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_EEEM_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_EEEM_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_EEEM_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_EEEM_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_EEEM_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_EEEM_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_EEEM->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }


    //EEMM
    i=1; 
    leptonscands_EEMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEMM;
    iEvent.getByLabel(RECOcollNameEEMM, CandidatesEEMM);
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEMM->begin();cand != CandidatesEEMM->end(); ++ cand ) { 
      //cout << "Found a EEMuMu candidate "<< endl;
      if (kk>199) continue;
      RECO_EEMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEMM_PT[i-1][kk]=cand->p4().pt();
      RECO_EEMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEMM_PHI[i-1][kk]=cand->p4().phi();
      RECO_EEMM_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_EEMM_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_EEMM_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_EEMM_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_EEMM_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_EEMM_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_EEMM_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_EEMM_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_EEMM->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }
    
    //MMEE
    i=1; 
    leptonscands_MMEE->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMEE;
    iEvent.getByLabel(RECOcollNameMMEE, CandidatesMMEE);
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMEE->begin();cand != CandidatesMMEE->end(); ++ cand ) { 
      //cout << "Found a EEMuMu candidate "<< endl;
      if (kk>199) continue;
      RECO_MMEE_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMEE_PT[i-1][kk]=cand->p4().pt();
      RECO_MMEE_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMEE_PHI[i-1][kk]=cand->p4().phi();
      RECO_MMEE_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMEE_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMEE_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMEE_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMEE_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_MMEE_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMEE_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMEE_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMEE_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMEE_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_MMEE_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_MMEE_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_MMEE_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_MMEE_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_MMEE_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_MMEE_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_MMEE->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    //EEEE
    i=1; 
    leptonscands_EEEE->clear();
    edm::Handle<edm::View<Candidate> > CandidatesEEEE;
    iEvent.getByLabel(RECOcollNameEEEE, CandidatesEEEE);
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesEEEE->begin();cand != CandidatesEEEE->end(); ++ cand ) { 
      //cout << "Found a EEEE candidate "<< endl;
      if (kk>199) continue;
      RECO_EEEE_MASS[i-1][kk]=cand->p4().mass();
      RECO_EEEE_PT[i-1][kk]=cand->p4().pt();
      RECO_EEEE_ETA[i-1][kk]=cand->p4().eta();
      RECO_EEEE_PHI[i-1][kk]=cand->p4().phi();
      RECO_EEEE_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_EEEE_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_EEEE_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_EEEE_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_EEEE_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_EEEE_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_EEEE_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_EEEE_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_EEEE_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_EEEE_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_EEEE_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_EEEE_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_EEEE_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_EEEE_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_EEEE_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_EEEE_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_EEEE->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    //MMMM
    i=1; 
    leptonscands_MMMM->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMMMM;
    iEvent.getByLabel(RECOcollNameMMMM, CandidatesMMMM);
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMMMM->begin();cand != CandidatesMMMM->end(); ++ cand ) { 
      //cout << "Found a MMMM candidate "<< endl;
      if (kk>199) continue;
      RECO_MMMM_MASS[i-1][kk]=cand->p4().mass();
      RECO_MMMM_PT[i-1][kk]=cand->p4().pt();
      RECO_MMMM_ETA[i-1][kk]=cand->p4().eta();
      RECO_MMMM_PHI[i-1][kk]=cand->p4().phi();
      RECO_MMMM_CHARGE[i-1][kk]=cand->charge();
      int l=0;
      int n=0;
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MMMM_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MMMM_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MMMM_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MMMM_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_MMMM_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
	for (unsigned k = 0; k < cand->daughter(j)->numberOfDaughters(); ++k ) {
	  RECO_MMMM_MASS[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().mass();
	  RECO_MMMM_PT[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().pt();
	  RECO_MMMM_ETA[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().eta();
	  RECO_MMMM_PHI[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->p4().phi();
	  RECO_MMMM_CHARGE[i+j+k+l+2][kk]=cand->daughter(j)->daughter(k)->charge();
	  //cout << "index" << i+j+k+l+2 <<endl;
	  for (unsigned m = 0; m < cand->daughter(j)->daughter(k)->numberOfDaughters(); ++m ) {
	    //cout << "1,2,3,4= " << cand->daughter(j)->daughter(k)->daughter(m)->p4().pt() << endl; 
	    RECO_MMMM_MASS[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().mass();
	    RECO_MMMM_PT[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().pt();
	    RECO_MMMM_ETA[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().eta();
	    RECO_MMMM_PHI[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->p4().phi();
	    RECO_MMMM_CHARGE[i+j+k+l+m+n+4][kk]=cand->daughter(j)->daughter(k)->daughter(m)->charge();
	    leptonscands_MMMM->push_back( cand->daughter(j)->daughter(k)->daughter(m)->clone());
	    //cout << "index" << i+j+k+l+m+n+4 <<endl;
	  }
	  n++; 
	} 
	l++;
      }
      kk++;
    }

    //MT
    i=1; 
    leptonscands_MT->clear();
    edm::Handle<edm::View<Candidate> > CandidatesMT;
    iEvent.getByLabel(RECOcollNameMT, CandidatesMT);
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesMT->begin();cand != CandidatesMT->end(); ++ cand ) { 
      //cout << "Found a MT candidate "<< endl;
      if (kk>199) continue;
      RECO_MT_MASS[i-1][kk]=cand->p4().mass();
      RECO_MT_PT[i-1][kk]=cand->p4().pt();
      RECO_MT_ETA[i-1][kk]=cand->p4().eta();
      RECO_MT_PHI[i-1][kk]=cand->p4().phi();
      RECO_MT_CHARGE[i-1][kk]=cand->charge();
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_MT_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_MT_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_MT_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_MT_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_MT_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
      }
      kk++;
    }
    
    
    //EE
    edm::Handle<edm::View<Candidate> > CandidatesEE;
    iEvent.getByLabel(RECOcollNameZ.at(1), CandidatesEE);
    for( edm::View<Candidate>::const_iterator cand = CandidatesEE->begin();cand != CandidatesEE->end(); ++ cand ) {
      //cout<<"mass of EE pair====="<<cand->p4().mass()<<endl;
    } 

    //ET
    i=1; 
    leptonscands_ET->clear();
    edm::Handle<edm::View<Candidate> > CandidatesET;
    iEvent.getByLabel(RECOcollNameET, CandidatesET);
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesET->begin();cand != CandidatesET->end(); ++ cand ) {
      //cout<<"mass of ET pair====="<<cand->p4().mass()<<endl;
      //cout << "Found a ET candidate "<< endl;
      if (kk>199) continue;
      RECO_ET_MASS[i-1][kk]=cand->p4().mass();
      RECO_ET_PT[i-1][kk]=cand->p4().pt();
      RECO_ET_ETA[i-1][kk]=cand->p4().eta();
      RECO_ET_PHI[i-1][kk]=cand->p4().phi();
      RECO_ET_CHARGE[i-1][kk]=cand->charge();
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_ET_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_ET_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_ET_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_ET_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_ET_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
      }
      kk++;
    }
    
    //ME
    i=1; 
    leptonscands_ME->clear();
    edm::Handle<edm::View<Candidate> > CandidatesME;
    iEvent.getByLabel(RECOcollNameME, CandidatesME);
    kk=0;
    for( edm::View<Candidate>::const_iterator cand = CandidatesME->begin();cand != CandidatesME->end(); ++ cand ) { 
      //cout << "Found a ME candidate "<< endl;
      if (kk>199) continue;
      RECO_ME_MASS[i-1][kk]=cand->p4().mass();
      RECO_ME_PT[i-1][kk]=cand->p4().pt();
      RECO_ME_ETA[i-1][kk]=cand->p4().eta();
      RECO_ME_PHI[i-1][kk]=cand->p4().phi();
      RECO_ME_CHARGE[i-1][kk]=cand->charge();
      //cout << "index" << i-1 << endl;
      for (unsigned j = 0; j < cand->numberOfDaughters(); ++j ) {
	RECO_ME_MASS[i+j][kk]=cand->daughter(j)->p4().mass();
	RECO_ME_PT[i+j][kk]=cand->daughter(j)->p4().pt();
	RECO_ME_ETA[i+j][kk]=cand->daughter(j)->p4().eta();
	RECO_ME_PHI[i+j][kk]=cand->daughter(j)->p4().phi();
	RECO_ME_CHARGE[i+j][kk]=cand->daughter(j)->charge();
	//cout << "index" << i+j <<endl;
      }
      kk++;
    }  
    
    //Trans Electron & Muon
    edm::Handle<std::vector<reco::PFMET> > pfmet;
    iEvent.getByLabel("pfMet", pfmet);
    double PxMET = (*pfmet)[0].px();
    double PyMET = (*pfmet)[0].py();
    
    edm::Handle <GsfElectronCollection> gsfElectronCollection;
    iEvent.getByLabel("hTozzTo4leptonsElectronSelector", gsfElectronCollection);
    if (!(iEvent.getByLabel("hTozzTo4leptonsElectronSelector", gsfElectronCollection))) {
      std::cout<<"No electron Coll "<<std::endl;
    }
    unsigned int GsfElectronCollectionSize = gsfElectronCollection->size();
    
    for (unsigned int m=0; m<GsfElectronCollectionSize; ++m) {
      const GsfElectron & elec = gsfElectronCollection->at(m);
      reco::Candidate::LorentzVector VisParticleELE = elec.p4();
      double TransEle= compMt(VisParticleELE,PxMET,PyMET);
      RECOELE_TMASS[m]=TransEle;
    }
    
    edm::Handle<MuonCollection> muonCollection;
    iEvent.getByLabel("hTozzTo4leptonsMuonSelector", muonCollection);  
    if (!(iEvent.getByLabel("hTozzTo4leptonsMuonSelector", muonCollection))) {
      std::cout<<"No muon Coll "<<std::endl;
    }
    unsigned int muonCollectionSize = muonCollection->size();
    for (unsigned int m=0; m<muonCollectionSize; ++m) {
      const Muon & mu = muonCollection->at(m);
      reco::Candidate::LorentzVector VisParticleMU = mu.p4();
      double TransMu= compMt(VisParticleMU,PxMET,PyMET);
      RECOMU_TMASS[m]=TransMu;
      cout <<"RECOMU_TMASS[===="<<RECOMU_TMASS[m]<<endl;
    }
  }
  
  

  void fillHPSTaus(const edm::Event& iEvent){
    
    edm::Handle<PFTauCollection > thehpsTauHandle;
    iEvent.getByLabel(edm::InputTag("hpsPFTauProducer",""),thehpsTauHandle);
    const PFTauCollection& myhpsTauCollection = *(thehpsTauHandle.product()); 
    
    edm::Handle<PFTauDiscriminator> hpsdiscLER;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByLooseElectronRejection",""), hpsdiscLER);

    edm::Handle<PFTauDiscriminator> hpsdiscMER;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByMediumElectronRejection",""), hpsdiscMER);

    edm::Handle<PFTauDiscriminator> hpsdiscTER;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByTightElectronRejection",""), hpsdiscTER);
    
    edm::Handle<PFTauDiscriminator> hpsdiscLMR;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByLooseMuonRejection",""), hpsdiscLMR);

    edm::Handle<PFTauDiscriminator> hpsdiscTMR;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByTightMuonRejection",""), hpsdiscTMR);
    
    edm::Handle<PFTauDiscriminator> hpsdiscDecayF;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByDecayModeFinding",""), hpsdiscDecayF);
    
    edm::Handle<PFTauDiscriminator> hpsdiscLooseIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByLooseIsolation",""), hpsdiscLooseIso);

    edm::Handle<PFTauDiscriminator> hpsdiscVLooseIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByVLooseIsolation",""), hpsdiscVLooseIso);
    
    edm::Handle<PFTauDiscriminator> hpsdiscMediumIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByMediumIsolation",""), hpsdiscMediumIso);
    
    edm::Handle<PFTauDiscriminator> hpsdiscTightIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByTightIsolation",""), hpsdiscTightIso);

    edm::Handle<PFTauDiscriminator> hpsdiscLooseChargedIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByLooseChargedIsolation",""), hpsdiscLooseChargedIso);
    
    edm::Handle<PFTauDiscriminator> hpsdiscVLooseChargedIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByVLooseChargedIsolation",""), hpsdiscVLooseChargedIso);
    
    edm::Handle<PFTauDiscriminator> hpsdiscMediumChargedIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByMediumChargedIsolation",""), hpsdiscMediumChargedIso);
    
    edm::Handle<PFTauDiscriminator> hpsdiscTightChargedIso;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByTightChargedIsolation",""), hpsdiscTightChargedIso);
    
    edm::Handle<PFTauDiscriminator> hpsdiscLooseIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr",""), hpsdiscLooseIsoDB);
    
    edm::Handle<PFTauDiscriminator> hpsdiscVLooseIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr",""), hpsdiscVLooseIsoDB);
    
    edm::Handle<PFTauDiscriminator> hpsdiscMediumIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr",""), hpsdiscMediumIsoDB);
    
    edm::Handle<PFTauDiscriminator> hpsdiscTightIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByTightIsolationDBSumPtCorr",""), hpsdiscTightIsoDB);

    edm::Handle<PFTauDiscriminator> hpsdiscLooseCombinedIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr",""), hpsdiscLooseCombinedIsoDB);
    
    edm::Handle<PFTauDiscriminator> hpsdiscVLooseCombinedIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr",""), hpsdiscVLooseCombinedIsoDB);
    
    edm::Handle<PFTauDiscriminator> hpsdiscMediumCombinedIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr",""), hpsdiscMediumCombinedIsoDB);
    
    edm::Handle<PFTauDiscriminator> hpsdiscTightCombinedIsoDB;
    iEvent.getByLabel(edm::InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr",""), hpsdiscTightCombinedIsoDB);


    int indexhps=0;

    for (PFTauCollection::size_type ihpsTau = 0; ihpsTau < myhpsTauCollection.size(); ihpsTau++)  {
    PFTauRef thehpsTau(thehpsTauHandle, ihpsTau);
   
      if(indexhps>=100) continue;

      HPSTAU_MASS[indexhps]    = myhpsTauCollection.at(ihpsTau).p4().mass();
      HPSTAU_PT[indexhps]      = myhpsTauCollection.at(ihpsTau).p4().pt();
      HPSTAU_ETA[indexhps]     = myhpsTauCollection.at(ihpsTau).p4().eta();
      HPSTAU_PHI[indexhps]     = myhpsTauCollection.at(ihpsTau).p4().phi();
      HPSTAU_CHARGE[indexhps]  = myhpsTauCollection.at(ihpsTau).charge();

      double hpstau_discByLER = (*hpsdiscLER)[thehpsTau];
      double hpstau_discByMER = (*hpsdiscMER)[thehpsTau];
      double hpstau_discByTER = (*hpsdiscTER)[thehpsTau];
      double hpstau_discByLMR = (*hpsdiscLMR)[thehpsTau];
      double hpstau_discByTMR = (*hpsdiscTMR)[thehpsTau];
      double hpstau_discByDecayF = (*hpsdiscDecayF)[thehpsTau];
      double hpstau_discLooseIso = (*hpsdiscLooseIso)[thehpsTau];
      double hpstau_discVLooseIso = (*hpsdiscVLooseIso)[thehpsTau];
      double hpstau_discMediumIso = (*hpsdiscMediumIso)[thehpsTau];
      double hpstau_discTightIso = (*hpsdiscTightIso)[thehpsTau];

      double hpstau_discLooseChargedIso = (*hpsdiscLooseChargedIso)[thehpsTau];
      double hpstau_discVLooseChargedIso = (*hpsdiscVLooseChargedIso)[thehpsTau];
      double hpstau_discMediumChargedIso = (*hpsdiscMediumChargedIso)[thehpsTau];
      double hpstau_discTightChargedIso = (*hpsdiscTightChargedIso)[thehpsTau];
      
      double hpstau_discLooseIsoDB = (*hpsdiscLooseIsoDB)[thehpsTau];
      double hpstau_discVLooseIsoDB = (*hpsdiscVLooseIsoDB)[thehpsTau];
      double hpstau_discMediumIsoDB = (*hpsdiscMediumIsoDB)[thehpsTau];
      double hpstau_discTightIsoDB = (*hpsdiscTightIsoDB)[thehpsTau];
      
      double hpstau_discLooseCombinedIsoDB = (*hpsdiscLooseCombinedIsoDB)[thehpsTau];
      double hpstau_discVLooseCombinedIsoDB = (*hpsdiscVLooseCombinedIsoDB)[thehpsTau];
      double hpstau_discMediumCombinedIsoDB = (*hpsdiscMediumCombinedIsoDB)[thehpsTau];
      double hpstau_discTightCombinedIsoDB = (*hpsdiscTightCombinedIsoDB)[thehpsTau];
      
      HPSTAU_discByLER[indexhps]       = hpstau_discByLER;
      HPSTAU_discByMER[indexhps]       = hpstau_discByMER;
      HPSTAU_discByTER[indexhps]       = hpstau_discByTER;
      HPSTAU_discByLMR[indexhps]       = hpstau_discByLMR;
      HPSTAU_discByTMR[indexhps]       = hpstau_discByTMR;
      HPSTAU_discByDecayF[indexhps]    = hpstau_discByDecayF;
      HPSTAU_discLooseIso[indexhps]    = hpstau_discLooseIso;
      HPSTAU_discVLooseIso[indexhps]   = hpstau_discVLooseIso;
      HPSTAU_discMediumIso[indexhps]   = hpstau_discMediumIso;
      HPSTAU_discTightIso[indexhps]    = hpstau_discTightIso;

      HPSTAU_discLooseChargedIso[indexhps]    = hpstau_discLooseChargedIso;
      HPSTAU_discVLooseChargedIso[indexhps]   = hpstau_discVLooseChargedIso;
      HPSTAU_discMediumChargedIso[indexhps]   = hpstau_discMediumChargedIso;
      HPSTAU_discTightChargedIso[indexhps]    = hpstau_discTightChargedIso;
      
      HPSTAU_discLooseIsoDB[indexhps]    = hpstau_discLooseIsoDB;
      HPSTAU_discVLooseIsoDB[indexhps]   = hpstau_discVLooseIsoDB;
      HPSTAU_discMediumIsoDB[indexhps]   = hpstau_discMediumIsoDB;
      HPSTAU_discTightIsoDB[indexhps]    = hpstau_discTightIsoDB;
      
      HPSTAU_discLooseCombinedIsoDB[indexhps]    = hpstau_discLooseCombinedIsoDB;
      HPSTAU_discVLooseCombinedIsoDB[indexhps]   = hpstau_discVLooseCombinedIsoDB;
      HPSTAU_discMediumCombinedIsoDB[indexhps]   = hpstau_discMediumCombinedIsoDB;
      HPSTAU_discTightCombinedIsoDB[indexhps]    = hpstau_discTightCombinedIsoDB;
      
      indexhps++;
    } 
  }

  void fillElectronDisc(const edm::Event& iEvent){
    
    edm::Handle <GsfElectronCollection> ElectronCollection;//edm::View<reco::GsfElectron> 
    iEvent.getByLabel("gsfElectrons", ElectronCollection);
    //unsigned int GsfElectronCollectionSize = GsfElectronCollection->size();
    const GsfElectronCollection& myGsfElectronCollection = *(ElectronCollection.product());
    
    edm::Handle<edm::ValueMap<float> >eIDValueMap95;
    iEvent.getByLabel( "simpleEleId95relIso" , eIDValueMap95);
    
    edm::Handle<edm::ValueMap<float> >eIDValueMap90;
    iEvent.getByLabel( "simpleEleId95relIso" , eIDValueMap90);
    
    edm::Handle<edm::ValueMap<float> >eIDValueMap85;
    iEvent.getByLabel( "simpleEleId95relIso" , eIDValueMap85);
    
    edm::Handle<edm::ValueMap<float> >eIDValueMap80;
    iEvent.getByLabel( "simpleEleId95relIso" , eIDValueMap80); 
    
    int indexElectronDisc=0;
    
    //for (unsigned int i=0; i<GsfElectronCollectionSize; ++i) {//PFTauCollection::size_type iSTau = 0; iSTau < mySTauCollection.size(); iSTau++
    for(GsfElectronCollection::size_type iele=0; iele<myGsfElectronCollection.size(); ++iele) {
      //edm::Ref<edm::View<reco::GsfElectron> >  electronRef(GsfElectronCollection,iele);
      GsfElectronRef  electronRef(ElectronCollection,iele);
      
      if(indexElectronDisc>=100) continue;
      
      double simpleEleId95relIso = (*eIDValueMap95)[electronRef];
      double simpleEleId90relIso = (*eIDValueMap90)[electronRef];
      double simpleEleId85relIso = (*eIDValueMap85)[electronRef];
      double simpleEleId80relIso = (*eIDValueMap80)[electronRef];
      
      EleId95relIso[indexElectronDisc]        = simpleEleId95relIso;
      EleId90relIso[indexElectronDisc]        = simpleEleId90relIso;
      EleId85relIso[indexElectronDisc]        = simpleEleId85relIso;
      EleId80relIso[indexElectronDisc]        = simpleEleId80relIso;
      
      indexElectronDisc++;
    }
  }
  
  void fillElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    
    // Supercluster collection
    edm::Handle<reco::SuperClusterCollection> clusters;
    iEvent.getByLabel(clusterCollectionTag_,clusters);
    
    // GsfTracks
    edm::Handle<reco::GsfTrackCollection> gsfTracks;
    iEvent.getByLabel(gsftrackCollection_,gsfTracks);
    
    // Electrons
   /*  edm::Handle<edm::View<reco::GsfElectron> > EleCandidates; */
/*     iEvent.getByLabel(electronTag_, EleCandidates); */
    
/*     // Isolation */
/*     edm::Handle<edm::ValueMap<float> > isoTkelemap; */
/*     iEvent.getByLabel(electronTkMapTag_, isoTkelemap); */
    
/*     edm::Handle<edm::ValueMap<float> > isoEcalelemap; */
/*     iEvent.getByLabel(electronEcalMapTag_, isoEcalelemap); */
    
/*     edm::Handle<edm::ValueMap<float> > isoHcalelemap; */
/*     iEvent.getByLabel(electronHcalMapTag_, isoHcalelemap); */
    
/*     edm::Handle<edm::ValueMap<float> > isoelemap; */
/*     iEvent.getByLabel(electronMapTag_, isoelemap); */

    // EG isolation
    //edm::Handle<reco::GsfElectronRefVector> EleRefs;
    //iEvent.getByLabel(electronEgmTag_, EleRefs);
    edm::Handle<edm::View<reco::GsfElectron> > EleRefs;
    iEvent.getByLabel(electronEgmTag_, EleRefs);

    edm::Handle<edm::ValueMap<double> > egmisoTkelemap;
    iEvent.getByLabel(electronEgmTkMapTag_, egmisoTkelemap);
    
    edm::Handle<edm::ValueMap<double> > egmisoEcalelemap;
    iEvent.getByLabel(electronEgmEcalMapTag_, egmisoEcalelemap);
    
    edm::Handle<edm::ValueMap<double> > egmisoHcalelemap;
    iEvent.getByLabel(electronEgmHcalMapTag_, egmisoHcalelemap);

    // Particle Flow Isolation
    edm::Handle<edm::ValueMap<double> > isoPFAllelemap;
    iEvent.getByLabel(electronPFIsoValueAllTag_, isoPFAllelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFChargedelemap;
    iEvent.getByLabel(electronPFIsoValueChargedTag_, isoPFChargedelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFNeutralelemap;
    iEvent.getByLabel(electronPFIsoValueNeutralTag_, isoPFNeutralelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFGammaelemap;
    iEvent.getByLabel(electronPFIsoValueGammaTag_, isoPFGammaelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFPUelemap;
    iEvent.getByLabel(electronPFIsoValuePUTag_, isoPFPUelemap);
    
    edm::Handle<edm::ValueMap<double> > isoPFPULowelemap;
    iEvent.getByLabel(electronPFIsoValuePULowTag_, isoPFPULowelemap);

    //Electron ID CiC
    std::vector<edm::Handle<edm::ValueMap<float> > > eleIdCutHandles(4) ;
    iEvent.getByLabel  (EleID_VeryLooseTag_ , eleIdCutHandles[0]) ;
    iEvent.getByLabel  (EleID_LooseTag_ , eleIdCutHandles[1]) ;
    iEvent.getByLabel  (EleID_MediumTag_ , eleIdCutHandles[2]) ;
    iEvent.getByLabel  (EleID_TightTag_ , eleIdCutHandles[3]) ;

    //Electron ID HZZ CiC
    std::vector<edm::Handle<edm::ValueMap<float> > > eleIdHZZCutHandles(4) ;
    iEvent.getByLabel  (EleID_HZZVeryLooseTag_ , eleIdHZZCutHandles[0]) ;
    iEvent.getByLabel  (EleID_HZZLooseTag_ , eleIdHZZCutHandles[1]) ;
    iEvent.getByLabel  (EleID_HZZMediumTag_ , eleIdHZZCutHandles[2]) ;
    iEvent.getByLabel  (EleID_HZZTightTag_ , eleIdHZZCutHandles[3]) ;


    // Vertexing
    //3D DA
    edm::Handle<edm::View<reco::GsfElectron> > VertEleCandidates;
    iEvent.getByLabel(electronTag_Vert, VertEleCandidates);
    
    edm::Handle<edm::ValueMap<float> > vertexelemap;
    iEvent.getByLabel(electronMapTag_Vert, vertexelemap);
    
    edm::Handle<edm::ValueMap<float> > vertexelemapvalue;
    iEvent.getByLabel(electronMapTag_VertValue, vertexelemapvalue);
    
    edm::Handle<edm::ValueMap<float> > vertexelemaperror;
    iEvent.getByLabel(electronMapTag_VertError, vertexelemaperror);

    // 3D KF
    edm::Handle<edm::ValueMap<float> > vertexelemapKF;
    iEvent.getByLabel(electronMapTag_VertKF, vertexelemapKF);
    
    edm::Handle<edm::ValueMap<float> > vertexelemapvalueKF;
    iEvent.getByLabel(electronMapTag_VertValueKF, vertexelemapvalueKF);
    
    edm::Handle<edm::ValueMap<float> > vertexelemaperrorKF;
    iEvent.getByLabel(electronMapTag_VertErrorKF, vertexelemaperrorKF);


    // 3D w.r.t GD 2e2mu vertex
    edm::Handle<edm::ValueMap<float> > vertexelemapGD;
    iEvent.getByLabel(electronMapTag_VertGD, vertexelemapGD);
    // 3D w.r.t GD 4mu vertex 
    edm::Handle<edm::ValueMap<float> > vertexelemapGDEEEE;
    iEvent.getByLabel(electronMapTag_VertGDEEEE, vertexelemapGDEEEE);
    // 3D w.r.t std 2e2mu vertex
    edm::Handle<edm::ValueMap<float> > vertexelemapStd;
    iEvent.getByLabel(electronMapTag_VertStd, vertexelemapStd);
    // 3D w.r.t Std 4mu vertex
    edm::Handle<edm::ValueMap<float> > vertexelemapStdEEEE;
    iEvent.getByLabel(electronMapTag_VertStdEEEE, vertexelemapStdEEEE);
    // 3D w.r.t std 2e2mu vertex
    edm::Handle<edm::ValueMap<float> > vertexelemapKin;
    iEvent.getByLabel(electronMapTag_VertKin, vertexelemapKin);
     // 3D w.r.t Kin 4mu vertex
    edm::Handle<edm::ValueMap<float> > vertexelemapKinEEEE;
    iEvent.getByLabel(electronMapTag_VertKinEEEE, vertexelemapKinEEEE);




    // STIP SLIP
    edm::Handle<edm::ValueMap<float> > stipelemap;
    iEvent.getByLabel(electronSTIPMapTag_Vert, stipelemap);
    
    edm::Handle<edm::ValueMap<float> > slipelemap;
    iEvent.getByLabel(electronSLIPMapTag_Vert, slipelemap);
    
    edm::Handle<edm::ValueMap<float> > stipelemapvalue;
    iEvent.getByLabel(electronSTIPMapTag_VertValue, stipelemapvalue);
    
    edm::Handle<edm::ValueMap<float> > slipelemapvalue;
    iEvent.getByLabel(electronSLIPMapTag_VertValue, slipelemapvalue);
    
    edm::Handle<edm::ValueMap<float> > stipelemaperror;
    iEvent.getByLabel(electronSTIPMapTag_VertError, stipelemaperror);
    
    edm::Handle<edm::ValueMap<float> > slipelemaperror;
    iEvent.getByLabel(electronSLIPMapTag_VertError, slipelemaperror);		
    
    // Conversion Finder
    edm::Handle<edm::ValueMap<float> > conversiondistmap;
    iEvent.getByLabel(ConvMapDistTag_, conversiondistmap);

    edm::Handle<edm::ValueMap<float> > conversiondcotmap;
    iEvent.getByLabel(ConvMapDcotTag_, conversiondcotmap);


    // beam spot
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
    iEvent.getByType(recoBeamSpotHandle) ;
    const reco::BeamSpot bs = *recoBeamSpotHandle ;

    // primary vertex
    edm::Handle<VertexCollection> recoPVCollection;
    iEvent.getByLabel(verticesTag_, recoPVCollection);
    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (recoPVCollection->size() != 0);
    if(pvfound) {
      primVertex = recoPVCollection->front();
      pVertex = math::XYZPoint(primVertex.position().x(), primVertex.position().y(), primVertex.position().z());
    }

    //ECAL Reduced RecHit Input Tag
    edm::Handle< EcalRecHitCollection > reducedEBRecHits;
    edm::Handle< EcalRecHitCollection > reducedEERecHits;

    string EcalRecHitsEB="";
    string EcalRecHitsEE="";

    if (useRECOformat) {
       EcalRecHitsEB="EcalRecHitsEB";
       EcalRecHitsEE="EcalRecHitsEE";
    }
    else {
      EcalRecHitsEB="reducedEcalRecHitsEB";
      EcalRecHitsEE="reducedEcalRecHitsEE";
    }
    
    iEvent.getByLabel( edm::InputTag(EcalRecHitsEB.c_str()), reducedEBRecHits );
    iEvent.getByLabel( edm::InputTag(EcalRecHitsEE.c_str()), reducedEERecHits );
    
    //calo topology
    const CaloTopology * topology ;
    const EcalChannelStatus *chStatus ;
    unsigned long long cacheIDTopo_=0;
    edm::ESHandle<CaloTopology> theCaloTopo;
    if (cacheIDTopo_!=iSetup.get<CaloTopologyRecord>().cacheIdentifier()){
      cacheIDTopo_=iSetup.get<CaloTopologyRecord>().cacheIdentifier();
      iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
    }
    topology = theCaloTopo.product() ;
    edm::ESHandle<EcalChannelStatus> pChannelStatus;
    iSetup.get<EcalChannelStatusRcd>().get(pChannelStatus);
    chStatus = pChannelStatus.product();
    
    
    int index=0;
    RECO_NELE=EleRefs->size();
   

    for (edm::View<reco::GsfElectron>::const_iterator cand = EleRefs->begin(); 
	 cand != EleRefs->end(); ++cand) {
      
      if(index>=100) continue;

      edm::Ref<edm::View<reco::GsfElectron> > theele(EleRefs,index);
      edm::Ref<edm::View<reco::GsfElectron> > eletrackref(EleRefs,index);
      edm::Ref<edm::View<reco::GsfElectron> > eletrackrefv(VertEleCandidates,index);
  
     
      // Global variables
      RECOELE_isEcalDriven[index]    = cand->ecalDrivenSeed();
      RECOELE_isTrackerDriven[index] = cand->trackerDrivenSeed();

      std::cout << "\n Electron in the event: "
		<< "  isEcalDriven="    << RECOELE_isEcalDriven[index]  
		<< "  isTrackerDriven=" << RECOELE_isTrackerDriven[index] 
		<< std::endl;
      
      // kinematic
      RECOELE_E[index]       = cand->p4().energy();
      RECOELE_PT[index]      = cand->p4().pt();
      RECOELE_P[index]=sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz());
      RECOELE_ETA[index]     = cand->p4().eta();
      RECOELE_THETA[index]   = cand->p4().theta();
      RECOELE_PHI[index]     = cand->p4().phi();
      RECOELE_MASS[index]    = cand->p4().mass();
      RECOELE_QUALITY[index] = -999.;
      RECOELE_CHARGE[index]  = cand->charge();


      std::cout << "--kinematic:" 		
		<< "  pT="     << RECOELE_PT[index]  
		<< "  E="      << RECOELE_E[index]  
		<< "  p="      << RECOELE_P[index] 
		<< "  eta="    << RECOELE_ETA[index] 
		<< "  theta="  << RECOELE_THETA[index] 
		<< "  phi="    << RECOELE_PHI[index] 
		<< "  mass="   << RECOELE_MASS[index] 
		<< "  charge=" << RECOELE_CHARGE[index] 
		<< std::endl;

      
      // HZZ ele Isolation 
/*       RECOELE_TRACKISO[index]=(*isoTkelemap)[eletrackref]; */
/*       RECOELE_ECALISO[index]=(*isoEcalelemap)[eletrackref]; */
/*       RECOELE_HCALISO[index]=(*isoHcalelemap)[eletrackref]; */
/*       RECOELE_X[index]=(*isoelemap)[eletrackref]; */

/*       std::cout << "--isolation: electron" */
/* 		<< "  X="     << RECOELE_X[index] */
/* 		<< "  Track=" << RECOELE_TRACKISO[index] */
/* 		<< "  Ecal= " << RECOELE_ECALISO[index] */
/* 		<< "  Hcal=" << RECOELE_HCALISO[index] */
/* 		<< std::endl; */
    
      // Egamma isolation
      RECOELE_EGMTRACKISO[index]=(*egmisoTkelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMECALISO[index]=(*egmisoEcalelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMHCALISO[index]=(*egmisoHcalelemap)[eletrackref]/eletrackref->pt();
      // RECOELE_EGMX[index]=( (*egmisoTkelemap)[eletrackref] + (*egmisoEcalelemap)[eletrackref] + (*egmisoHcalelemap)[eletrackref])/eletrackref->pt();

      // temporary solution for the problem of rechits
      RECOELE_EGMECALISO[index]=(eletrackref->dr03EcalRecHitSumEt())/eletrackref->pt();
      RECOELE_EGMHCALISO[index]=(eletrackref->dr03HcalTowerSumEt())/eletrackref->pt();
      RECOELE_EGMX[index]=RECOELE_EGMTRACKISO[index]+RECOELE_EGMECALISO[index]+RECOELE_EGMHCALISO[index];
      
      std::cout << "--EG isolation: electron"
		<< "  X="     << RECOELE_EGMX[index]
		<< "  Track=" << RECOELE_EGMTRACKISO[index]
		<< "  Ecal= " << RECOELE_EGMECALISO[index]
		<< "  Hcal=" << RECOELE_EGMHCALISO[index]
		<< std::endl;

      // Particle Flow Isolation
      RECOELE_ALLDEPOSITS[index]      =  (*isoPFAllelemap)[theele]/theele->p4().pt(); 
      RECOELE_ChargedDEPOSITS[index]  =  (*isoPFChargedelemap)[theele]/theele->p4().pt();
      RECOELE_NeutralDEPOSITS[index]  =  (*isoPFNeutralelemap)[theele]/theele->p4().pt();
      RECOELE_GammaDEPOSITS[index]    =  (*isoPFGammaelemap)[theele]/theele->p4().pt();
      RECOELE_PUDEPOSITS[index]       =  (*isoPFPUelemap)[theele]/theele->p4().pt();
      RECOELE_PULowDEPOSITS[index]    =  (*isoPFPULowelemap)[theele]/theele->p4().pt();
      RECOELE_PFX_DB[index]           =  RECOELE_ChargedDEPOSITS[index] + max(RECOELE_GammaDEPOSITS[index]+RECOELE_NeutralDEPOSITS[index]-0.5*RECOELE_PUDEPOSITS[index] ,0.);
      
      std::cout << "Particle Flow isolation: electron" 
		<< "  \t RECOELE_ALLDEPOSITS="      << RECOELE_ALLDEPOSITS[index]
		<< "  \t RECOELE_ChargedDEPOSITS="  << RECOELE_ChargedDEPOSITS[index]
		<< "  \t RECOELE_NeutralDEPOSITS="  << RECOELE_NeutralDEPOSITS[index]
		<< "  \t RECOELE_GammaDEPOSITS="    << RECOELE_GammaDEPOSITS[index]
		<< "  \t RECOELE_PUDEPOSITS="       << RECOELE_PUDEPOSITS[index]
		<< "  \t RECOELE_PULowDEPOSITS="    << RECOELE_PULowDEPOSITS[index]
		<< "  \t RECOELE_PFX_DB="           << RECOELE_PFX_DB[index]
		<< std::endl; 
      
      
      // Vertexing DA
      RECOELE_SIP[index]=(*vertexelemap)[eletrackrefv];
      RECOELE_IP[index]=(*vertexelemapvalue)[eletrackrefv];
      RECOELE_IPERROR[index]=(*vertexelemaperror)[eletrackrefv];

      // KF
      RECOELE_SIP_KF[index]=(*vertexelemapKF)[eletrackrefv];
      RECOELE_IP_KF[index]=(*vertexelemapvalueKF)[eletrackrefv];
      RECOELE_IPERROR_KF[index]=(*vertexelemaperrorKF)[eletrackrefv];


      /* RECOELE_SIP_GD[index]=(*vertexelemapGD)[eletrackrefv]; */
/*       if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_GDEEEE[index]=(*vertexelemapGDEEEE)[eletrackrefv]; */
/*       RECOELE_SIP_Std[index]=(*vertexelemapStd)[eletrackrefv];  */
/*       if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_StdEEEE[index]=(*vertexelemapStdEEEE)[eletrackrefv];  */
/*       RECOELE_SIP_Kin[index]=(*vertexelemapKin)[eletrackrefv];  */
/*       if (decaychannel=="4e" || decaychannel=="2e2mu" ) RECOELE_SIP_KinEEEE[index]=(*vertexelemapKinEEEE)[eletrackrefv];  */
      

      RECOELE_STIP[index]=(*stipelemap)[eletrackrefv];
      RECOELE_SLIP[index]=(*slipelemap)[eletrackrefv];
      RECOELE_TIP[index]=(*stipelemapvalue)[eletrackrefv] ;
      RECOELE_LIP[index]=(*slipelemapvalue)[eletrackrefv];
      RECOELE_TIPERROR[index]=(*stipelemaperror)[eletrackrefv] ;
      RECOELE_LIPERROR[index]=(*slipelemaperror)[eletrackrefv];
            
      std::cout << "--vertexing: electron" 
		<< "  Sign_3DIP=" << RECOELE_SIP[index]
		<< "  3DIP="      << RECOELE_IP[index]
		<< "  3DIPerr="   << RECOELE_IPERROR[index] 	     
		<< "  Sign_TIP="  << RECOELE_STIP[index]  
		<< "  TIP="       << RECOELE_TIP[index] 
		<< "  TIPerror="  << RECOELE_TIPERROR[index] 
		<< "  Sign_LIP="  << RECOELE_SLIP[index]
		<< "  LIP="       << RECOELE_LIP[index] 
		<< "  LIPerror="  << RECOELE_LIPERROR[index]
		<< std::endl;


      // GsfTrack
      RECOELE_gsftrack_chi2[index]      = cand->gsfTrack()->normalizedChi2();
      RECOELE_gsftrack_dxyB[index]      = cand->gsfTrack()->dxy(bs.position()) ;
      RECOELE_gsftrack_dxy[index]       = cand->gsfTrack()->dxy(pVertex);
      RECOELE_gsftrack_dxyError[index]  = cand->gsfTrack()->dxyError();
      RECOELE_gsftrack_dzB[index]       = cand->gsfTrack()->dz(bs.position());
      RECOELE_gsftrack_dz[index]        = cand->gsfTrack()->dz(pVertex);
      RECOELE_gsftrack_dzError[index]   = cand->gsfTrack()->dzError();

      //Conversion variables
      RECOELE_gsftrack_losthits[index]  = cand->gsfTrack()->numberOfLostHits();
      RECOELE_gsftrack_validhits[index] = cand->gsfTrack()->numberOfValidHits();
      RECOELE_gsftrack_expected_inner_hits[index] = cand->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

      std::cout << "--gfstrack properties: " 
	<< "  chi2="     << RECOELE_gsftrack_chi2[index] 
	<< "  dxy="      << RECOELE_gsftrack_dxy[index] 
	<< "  dxyB="     << RECOELE_gsftrack_dxyB[index] 
	<< "  dxyError=" << RECOELE_gsftrack_dxyError[index] 
	<< "  dz="       << RECOELE_gsftrack_dz[index]
	<< "  dzB="      << RECOELE_gsftrack_dzB[index]
        << "  dzError="  << RECOELE_gsftrack_dzError[index] 
	<< "  losthits=" << RECOELE_gsftrack_losthits[index] 
	<< "  validhits="<< RECOELE_gsftrack_validhits[index] 
	<< std::endl;

      // Track-Cluster matching attributes
      RECOELE_ep[index]             = cand->eSuperClusterOverP(); 
      RECOELE_eSeedp[index]         = cand->eSeedClusterOverP();     
      RECOELE_eSeedpout[index]      = cand->eSeedClusterOverPout(); 
      RECOELE_eElepout[index]       = cand->eEleClusterOverPout();        
      
      RECOELE_deltaEtaIn [index]  = cand->deltaEtaSuperClusterTrackAtVtx(); 
      RECOELE_deltaEtaSeed[index] = cand->deltaEtaSeedClusterTrackAtCalo(); 
      RECOELE_deltaEtaEle[index]  = cand->deltaEtaEleClusterTrackAtCalo();  
      RECOELE_deltaPhiIn[index]   = cand->deltaPhiSuperClusterTrackAtVtx();
      RECOELE_deltaPhiSeed[index] = cand->deltaPhiSeedClusterTrackAtCalo(); 
      RECOELE_deltaPhiEle[index]  = cand->deltaPhiEleClusterTrackAtCalo() ;  

      std::cout << "--track-cluster matching: " 
	<< "  eSC_p="        << RECOELE_ep[index] 
	<< "  eSeed_p="      << RECOELE_eSeedp[index] 
	<< "  eSeed_pout="   << RECOELE_eSeedpout[index] 
	<< "  eEle_pout="    << RECOELE_eElepout[index] 
	<< "  deltaEtaIn="   << RECOELE_deltaEtaIn[index] 
        << "  deltaEtaSeed=" << RECOELE_deltaEtaSeed[index] 
        << "  deltaEtaEle="  << RECOELE_deltaEtaEle[index] 
        << "  deltaPhiIn="   << RECOELE_deltaPhiIn[index] 
        << "  deltaPhiSeed=" << RECOELE_deltaPhiSeed[index] 
	<< "  deltaPhiEle="  << RECOELE_deltaPhiEle[index] 
	<< std::endl;
      
       // Fiducial flags
      if (cand->isEB()) RECOELE_isbarrel[index] = 1 ; 
      else  RECOELE_isbarrel[index] = 0 ;
      if (cand->isEE()) RECOELE_isendcap[index] = 1 ; 
      else  RECOELE_isendcap[index] = 0 ;
      if (cand->isEBEtaGap())  RECOELE_isEBetaGap[index]  = 1 ;  
      if (cand->isEBPhiGap())  RECOELE_isEBphiGap[index]  = 1 ;  
      if (cand->isEEDeeGap())  RECOELE_isEEdeeGap[index]  = 1 ;  
      if (cand->isEERingGap()) RECOELE_isEEringGap[index] = 1 ;
 
      std::cout << "--fiducial flags: " 
	<< "  isEB="        << RECOELE_isbarrel[index] 
	<< "  isEBetagap="  << RECOELE_isEBetaGap[index] 
	<< "  isEBphigap="  << RECOELE_isEBphiGap[index] 
	<< "  isEE="        << RECOELE_isendcap[index] 
        << "  isEEdeegap="  << RECOELE_isEEdeeGap[index] 
	<< "  isEEringgap=" << RECOELE_isEEringGap[index] 
	<< std::endl;
      

      // Shower shape
      RECOELE_sigmaIetaIeta[index] = cand->sigmaIetaIeta() ; 
      RECOELE_sigmaEtaEta[index]   = cand->sigmaEtaEta() ;
      RECOELE_e15[index]           = cand->e1x5() ;
      RECOELE_e25max[index]        = cand->e2x5Max() ;
      RECOELE_e55[index]           = cand->e5x5() ;
      RECOELE_he[index]            = cand->hadronicOverEm() ;

      std::cout << "--shower shape:" 
	<< "  sigmaIetaIeta=" << RECOELE_sigmaIetaIeta[index] 
	<< "  sigmaEtaEta=" << RECOELE_sigmaEtaEta[index]  
	<< "  e15=" << RECOELE_e15[index] 
	<< "  e25max=" << RECOELE_e25max[index] 
	<< "  e55=" << RECOELE_e55[index]  
	<< "  he=" << RECOELE_he[index] 
	<< std::endl;
      
      // also variables on hcal
      // Isolation
      // Particle flow
      RECOELE_mva[index] = cand->mva();
      std::cout << "--PF: mva = " << RECOELE_mva[index] << std::endl;

      // Brem & Classifaction
      RECOELE_fbrem[index]         = cand->fbrem() ;
      RECOELE_nbrems[index]        = cand->numberOfBrems();
      RECOELE_Class[index]         = cand->classification() ;
      if (useRECOformat) RECOELE_fbrem_mean[index]=1. - cand->gsfTrack()->outerMomentum().R()/cand->gsfTrack()->innerMomentum().R();
      RECOELE_fbrem_mode[index]=cand->fbrem();

      std::cout << "--brem & classification: fbrem/nbrems/Class/fbrem_mean/fbrem_mode = "
		<< "  fbrem="      << RECOELE_fbrem[index]
		<< "  nbrems="     << RECOELE_nbrems[index]
        	<< "  class="      << RECOELE_Class[index]
		<< "  fbrem_mean=" << RECOELE_fbrem_mean[index]
		<< "  fbrem_mode=" << RECOELE_fbrem_mode[index]
		<< std::endl;

      // Corrections
      /*        RECOELE_isEcalEnergyCorrected[index] = cand->isEcalEnergyCorrected(); */
      /*        RECOELE_ecalEnergy[index]            = cand->ecalEnergy(); */
      /*        RECOELE_ecalEnergyError[index]       = cand->ecalEnergyError(); */
      /*        RECOELE_isMomentumCorrected[index]   = cand->isMomentumCorrected(); */
      /*        RECOELE_trackMomentumError[index]    = cand->trackMomentumError(); */
      /*        RECOELE_electronMomentumError[index] = cand->electronMomentumError(); */
      

      // SuperCluster
      reco::SuperClusterRef sclRef = cand->superCluster();
      math::XYZPoint sclPos = cand->superClusterPosition();

      if(!cand->ecalDrivenSeed() && cand->trackerDrivenSeed())
	sclRef = cand->pflowSuperCluster();
      double R  = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y() +sclRef->z()*sclRef->z());
      double Rt = TMath::Sqrt(sclRef->x()*sclRef->x() + sclRef->y()*sclRef->y());
      ele_sclRawE[index]   = sclRef->rawEnergy() ;
      RECOELE_scl_E[index]   = sclRef->energy() ;
      RECOELE_scl_Et[index]  = sclRef->energy()*(Rt/R) ;
      RECOELE_scl_Eta[index] = sclRef->eta() ;
      RECOELE_scl_Phi[index] = sclRef->phi() ;
      ele_sclX[index]  = sclPos.X();
      ele_sclY[index] =  sclPos.Y();
      ele_sclZ[index] =  sclPos.Z();
			
      std::cout << "--supercluster properties:"
	<< "  E_sc="    << RECOELE_scl_E[index]
	<< "  Et_sc="   << RECOELE_scl_Et[index]
	<< "  Eeta_sc=" << RECOELE_scl_Eta[index]
	<< "  phi_sc="  << RECOELE_scl_Phi[index]
	<< std::endl;

      // Seed Collection
      if (useRECOformat) {
	edm::RefToBase<TrajectorySeed> seed = cand->gsfTrack()->extra()->seedRef();
	reco::ElectronSeedRef MyS = seed.castTo<reco::ElectronSeedRef>();
	ele_seedSubdet2[index] = int(MyS->subDet2());
	if(fabs(MyS->dPhi2()) < 100.) ele_seedDphi2[index] = double(MyS->dPhi2());
	if(fabs(MyS->dRz2()) < 100.)  ele_seedDrz2[index]  = double(MyS->dRz2());
	
	ele_seedSubdet1[index] = int(MyS->subDet1());
	if(fabs(MyS->dPhi1()) < 100.) ele_seedDphi1[index] = double(MyS->dPhi1());
	if(fabs(MyS->dRz1()) < 100.)  ele_seedDrz1[index]  = double(MyS->dRz1());
      }
      
      //Spikes Variables
/*       if (useRECOformat) { */
/* 	if (cand->isEB()) { */
/* 	  const EcalRecHitCollection * reducedRecHits = 0 ; */
/* 	  reducedRecHits = reducedEBRecHits.product() ; */
	  
/* 	  //seed cluster analysis */
/* 	  const edm::Ptr<reco::CaloCluster> & seedCluster = cand->superCluster()->seed() ; */
/* 	  std::pair<DetId, float> id = EcalClusterTools::getMaximum(seedCluster->hitsAndFractions(),reducedRecHits); */
/* 	  const EcalRecHit & rh = getRecHit(id.first,reducedRecHits); */
/* 	  int flag = rh.recoFlag(); */
/* 	  if (flag == EcalRecHit::kOutOfTime) */
/* 	    ele_outOfTimeSeed[index] = 1; */
/* 	  else */
/* 	    ele_outOfTimeSeed[index] = 0; */
/* 	  int sev = EcalSeverityLevelAlgo::severityLevel(id.first,*reducedRecHits,*chStatus, 5., EcalSeverityLevelAlgo::kSwissCross,0.95) ; */
/* 	  ele_severityLevelSeed[index] = sev ; */
	  
/* 	  int dummyFlag = 0; */
/* 	  int dummySev = 0; */
/* 	  for (reco::CaloCluster_iterator bc = cand->superCluster()->clustersBegin(); */
/* 	       bc!=cand->superCluster()->clustersEnd(); */
/* 	       ++bc) { */
	    
/* 	    if ( seedCluster==(*bc) ) continue; */
/* 	    std::pair<DetId, float> id = EcalClusterTools::getMaximum((*bc)->hitsAndFractions(),reducedRecHits); */
/* 	    const EcalRecHit & rh = getRecHit(id.first,reducedRecHits); */
/* 	    int flag = rh.recoFlag(); */
/* 	    if (flag == EcalRecHit::kOutOfTime) */
/* 	      dummyFlag = 1 ; */
/* 	    int sev = EcalSeverityLevelAlgo::severityLevel(id.first,*reducedRecHits,*chStatus, 5., EcalSeverityLevelAlgo::kSwissCross,0.95); */
/* 	    if (sev > dummySev) */
/* 	      dummySev = sev ; */
/* 	  } */
/* 	  ele_severityLevelClusters[index] = dummySev ; */
/* 	  ele_outOfTimeClusters[index] = dummyFlag ; */
/* 	  ele_e2overe9[index] = E2overE9( id.first,*reducedRecHits,5,5, true, true); */
/* 	} */
/* 	else { */
/* 	  ele_e2overe9[index] = 0; */
/* 	  ele_severityLevelSeed[index] = 0 ; */
/* 	  ele_outOfTimeSeed[index] = 0 ; */
/* 	  ele_severityLevelClusters[index] = 0 ; */
/* 	  ele_outOfTimeClusters[index] = 0 ; */
/* 	} */
/*       } */
      
     
     
      //Old HZZ electron ID
      /* for (unsigned int i=0; i<eleIDTag_.size(); i++) { */
/* 	edm::Handle< edm::ValueMap<float> >  eIDValueMap; */
/* 	if (iEvent.getByLabel( eleIDTag_.at(i) , eIDValueMap)){ */
/* 	  const edm::ValueMap<float> & eIdmapCuts = * eIDValueMap ; */
/* 	  RECOELE_eleID[index][i]=eIdmapCuts[eletrackref]; */
	  
/* 	  std::cout << "eleID flag (loose,medium) = " */
/* 		    << RECOELE_eleID[index][0] << " "  */
/* 		    << RECOELE_eleID[index][1]  */
/* 		    << std::endl;  */
/* 	} */
/*       } */
 

      // CiC electron IDs
      ele_eidVeryLoose[index] = (*(eleIdCutHandles[0]))[eletrackref]; 
      ele_eidLoose[index] = (*(eleIdCutHandles[1]))[eletrackref]; 
      ele_eidMedium[index] = (*(eleIdCutHandles[2]))[eletrackref]; 
      ele_eidTight[index] = (*(eleIdCutHandles[3]))[eletrackref]; 
      
      std::cout << "CiC eleID flag (veryloose,loose,medium,tight) = "
		<< ele_eidVeryLoose[index] << " " 
		<< ele_eidLoose[index]     << " " 
		<< ele_eidMedium[index]    << " " 
		<< ele_eidTight[index]     << " " 
		<< std::endl; 


      // CiC HZZ eleID
      ele_eidHZZVeryLoose[index] = (*(eleIdHZZCutHandles[0]))[eletrackref]; 
      ele_eidHZZLoose[index] = (*(eleIdHZZCutHandles[1]))[eletrackref]; 
      ele_eidHZZMedium[index] = (*(eleIdHZZCutHandles[2]))[eletrackref]; 
      ele_eidHZZTight[index] = (*(eleIdHZZCutHandles[3]))[eletrackref]; 
      
      std::cout << "CiC HZZ eleID flag (veryloose,loose,medium,tight) = "
		<< ele_eidHZZVeryLoose[index] << " " 
		<< ele_eidHZZLoose[index]     << " " 
		<< ele_eidHZZMedium[index]    << " " 
		<< ele_eidHZZTight[index]     << " " 
		<< std::endl; 
      
      // Conversion Finder
      ConvMapDist[index]=(*conversiondistmap)[eletrackref];
      ConvMapDcot[index]=(*conversiondcotmap)[eletrackref];
      std::cout << "Converisone finder = "
		<< ConvMapDist[index] << " " 
		<< ConvMapDcot[index] << " " 
		<< std::endl;

      if (match(cand->mass(),cand->pt(),cand->charge(),leptonscands2e2mu_)) RECOELEBEST_2e2mu_MATCHED[index]=1;
      if (match(cand->mass(),cand->pt(),cand->charge(),leptonscands4e_))    RECOELEBEST_4e_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_EEMM)) RECOELE_EEMM_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_EEEE)) RECOELE_EEEE_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL0)) RECOELE_LLL0_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL1)) RECOELE_LLL1_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL2)) RECOELE_LLL2_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL3)) RECOELE_LLL3_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLLss0)) RECOELE_LLLLss0_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLLss1)) RECOELE_LLLLss1_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLLss2)) RECOELE_LLLLss2_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLl0)) RECOELE_LLLl0_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLl1)) RECOELE_LLLl1_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLL)) RECOELE_LLLL_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_Z1)) RECOELE_ZEE_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_Zss1)) RECOELE_ZssEE_MATCHED[index]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_Zcross)) RECOELE_ZEM_MATCHED[index]=1;
      
      index ++;
    }    
  }
  
  
  void fillMuons(const edm::Event& iEvent){
    // Muons
    edm::Handle<edm::View<reco::Muon> > MuCandidates;
    iEvent.getByLabel(muonTag_, MuCandidates);

    // Isolation
    edm::Handle<edm::ValueMap<double> > isoTkmumap;
    iEvent.getByLabel(muonTkMapTag_, isoTkmumap);
    
    edm::Handle<edm::ValueMap<double> > isoEcalmumap;
    iEvent.getByLabel(muonEcalMapTag_, isoEcalmumap);
    
    edm::Handle<edm::ValueMap<double> > isoHcalmumap;
    iEvent.getByLabel(muonHcalMapTag_, isoHcalmumap);

    //edm::Handle<edm::ValueMap<double> > isomumap;
    //iEvent.getByLabel(muonMapTag_, isomumap);
    
    // Particle Flow Isolation
    edm::Handle<edm::ValueMap<double> > isoPFAllmumap;
    iEvent.getByLabel(muonPFIsoValueAllTag_, isoPFAllmumap);
    const edm::ValueMap<double>& muPFAllIso = *isoPFAllmumap;

    edm::Handle<edm::ValueMap<double> > isoPFChargedmumap;
    iEvent.getByLabel(muonPFIsoValueChargedTag_, isoPFChargedmumap);

    edm::Handle<edm::ValueMap<double> > isoPFNeutralmumap;
    iEvent.getByLabel(muonPFIsoValueNeutralTag_, isoPFNeutralmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFGammamumap;
    iEvent.getByLabel(muonPFIsoValueGammaTag_, isoPFGammamumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFPUmumap;
    iEvent.getByLabel(muonPFIsoValuePUTag_, isoPFPUmumap);
    
    edm::Handle<edm::ValueMap<double> > isoPFPULowmumap;
    iEvent.getByLabel(muonPFIsoValuePULowTag_, isoPFPULowmumap);
     
    // Vertexing
    // 3D
    edm::Handle<edm::View<reco::Muon> > VertMuCandidates;
    iEvent.getByLabel(muonTag_Vert, VertMuCandidates);
    
    // 3D w.r.t primary vertex DA
    edm::Handle<edm::ValueMap<float> > vertexmumap;
    iEvent.getByLabel(muonMapTag_Vert, vertexmumap);
    
    edm::Handle<edm::ValueMap<float> > vertexmumapvalue;
    iEvent.getByLabel(muonMapTag_VertValue, vertexmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > vertexmumaperror;
    iEvent.getByLabel(muonMapTag_VertError, vertexmumaperror);

    // 3D w.r.t primary vertex KF
    edm::Handle<edm::ValueMap<float> > vertexmumapKF;
    iEvent.getByLabel(muonMapTag_VertKF, vertexmumapKF);
    
    edm::Handle<edm::ValueMap<float> > vertexmumapvalueKF;
    iEvent.getByLabel(muonMapTag_VertValueKF, vertexmumapvalueKF);
    
    edm::Handle<edm::ValueMap<float> > vertexmumaperrorKF;
    iEvent.getByLabel(muonMapTag_VertErrorKF, vertexmumaperrorKF);


    // 3D w.r.t GD 2e2mu vertex
    edm::Handle<edm::ValueMap<float> > vertexmumapGD;
    iEvent.getByLabel(muonMapTag_VertGD, vertexmumapGD);
    // 3D w.r.t GD 4mu vertex 
    edm::Handle<edm::ValueMap<float> > vertexmumapGDMMMM;
    iEvent.getByLabel(muonMapTag_VertGDMMMM, vertexmumapGDMMMM);
    // 3D w.r.t std 2e2mu vertex
    edm::Handle<edm::ValueMap<float> > vertexmumapStd;
    iEvent.getByLabel(muonMapTag_VertStd, vertexmumapStd);
    // 3D w.r.t Std 4mu vertex
    edm::Handle<edm::ValueMap<float> > vertexmumapStdMMMM;
    iEvent.getByLabel(muonMapTag_VertStdMMMM, vertexmumapStdMMMM);
    // 3D w.r.t std 2e2mu vertex
    edm::Handle<edm::ValueMap<float> > vertexmumapKin;
    iEvent.getByLabel(muonMapTag_VertKin, vertexmumapKin);
     // 3D w.r.t Kin 4mu vertex
    edm::Handle<edm::ValueMap<float> > vertexmumapKinMMMM;
    iEvent.getByLabel(muonMapTag_VertKinMMMM, vertexmumapKinMMMM);



    // STIP SLIP
    edm::Handle<edm::ValueMap<float> > stipmumap;
    iEvent.getByLabel(muonSTIPMapTag_Vert, stipmumap);
    
    edm::Handle<edm::ValueMap<float> > slipmumap;
    iEvent.getByLabel(muonSLIPMapTag_Vert, slipmumap);
    
    edm::Handle<edm::ValueMap<float> > stipmumapvalue;
    iEvent.getByLabel(muonSTIPMapTag_VertValue, stipmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > slipmumapvalue;
    iEvent.getByLabel(muonSLIPMapTag_VertValue, slipmumapvalue);
    
    edm::Handle<edm::ValueMap<float> > stipmumaperror;
    iEvent.getByLabel(muonSTIPMapTag_VertError, stipmumaperror);
    
    edm::Handle<edm::ValueMap<float> > slipmumaperror;
    iEvent.getByLabel(muonSLIPMapTag_VertError, slipmumaperror);
    
    // beam spot
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle ;
    iEvent.getByType(recoBeamSpotHandle) ;
    const reco::BeamSpot bs = *recoBeamSpotHandle ;
    
    // primary vertex
    edm::Handle<VertexCollection> recoPVCollection;
    iEvent.getByLabel(verticesTag_, recoPVCollection);
    Vertex primVertex;
    math::XYZPoint pVertex(0., 0., 0.);
    bool pvfound = (recoPVCollection->size() != 0);
    if(pvfound) {
      primVertex = recoPVCollection->front();
      pVertex = math::XYZPoint(primVertex.position().x(), primVertex.position().y(), primVertex.position().z());
    }
    
   
    
    
    int indexbis=0;
    RECO_NMU=MuCandidates->size();
    
    for (edm::View<reco::Muon>::const_iterator cand = MuCandidates->begin(); cand != MuCandidates->end(); ++cand) {
            
      if(indexbis>=100) continue; 
      
      unsigned int s=cand-MuCandidates->begin();

      edm::Ref<edm::View<reco::Muon> > themu(MuCandidates, indexbis);
      edm::Ref<edm::View<reco::Muon> > mutrackref(MuCandidates,indexbis); 
      edm::Ref<edm::View<reco::Muon> > mutrackrefv(VertMuCandidates,indexbis);
      
      // Global variables
      RECOMU_isGlobalMu[indexbis]=cand->isGlobalMuon();
      RECOMU_isStandAloneMu[indexbis]=cand->isStandAloneMuon();
      RECOMU_isTrackerMu[indexbis]=cand->isTrackerMuon();
      RECOMU_isCaloMu[indexbis]=cand->isCaloMuon();

      std::cout << "\n Muon in the event: "
		<<   "  isGB=" << RECOMU_isGlobalMu[indexbis]
		<< "   isSTA=" << RECOMU_isStandAloneMu[indexbis]
		<<   "  isTM=" << RECOMU_isTrackerMu[indexbis]
		<< "  isCalo=" << RECOMU_isCaloMu[indexbis]
		<< std::endl;

      // Kinematic of muon
      RECOMU_E[indexbis]=cand->p4().energy();
      RECOMU_PT[indexbis]=cand->p4().pt();
      RECOMU_P[indexbis]=sqrt(cand->p4().px()*cand->p4().px()+cand->p4().py()*cand->p4().py()+cand->p4().pz()*cand->p4().pz());
      RECOMU_ETA[indexbis]=cand->p4().eta();
      RECOMU_THETA[indexbis]=cand->p4().theta();
      RECOMU_PHI[indexbis]=cand->p4().phi();
      RECOMU_MASS[indexbis]=cand->p4().mass();
      RECOMU_QUALITY[indexbis]=-999.;
      RECOMU_CHARGE[indexbis]=cand->charge();

      std::cout << "--kinematic:"
		<< "  pT="     << RECOMU_PT[indexbis]
	        << "  E="      << RECOMU_E[indexbis]
		<< "  p="      << RECOMU_P[indexbis]
		<< "  eta="    << RECOMU_ETA[indexbis]
		<< "  theta="  << RECOMU_THETA[indexbis]
		<< "  phi="    << RECOMU_PHI[indexbis]
		<< "  mass="   << RECOMU_MASS[indexbis]
		<< "  charge=" << RECOMU_CHARGE[indexbis]
		<< std::endl;
	

      // Isolation
      RECOMU_TRACKISO[indexbis]=(*isoTkmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_ECALISO[indexbis]=(*isoEcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_HCALISO[indexbis]=(*isoHcalmumap)[mutrackref]/cand->p4().pt();
      //RECOMU_X[indexbis]=[(*isoTkmumap)[mutrackref]+(*isoEcalmumap)[mutrackref]+(*isoHcalmumap)[mutrackref]]/cand->p4().pt();
      // temporary solution for reducedrechit problem
      RECOMU_ECALISO[indexbis]=(cand->isolationR03().emEt)/cand->p4().pt();
      RECOMU_HCALISO[indexbis]=(cand->isolationR03().hadEt)/cand->p4().pt();
      RECOMU_X[indexbis]=RECOMU_TRACKISO[indexbis]+RECOMU_ECALISO[indexbis]+RECOMU_HCALISO[indexbis];

      std::cout << "--isolation: muon"
		<< "  X="      << RECOMU_X[indexbis]
		<< "  Track="  << RECOMU_TRACKISO[indexbis]
		<< "  Ecal="   << RECOMU_ECALISO[indexbis]
		<< "  Hcal="   << RECOMU_HCALISO[indexbis]
		<< std::endl;


      // Particle Flow Isolation
      RECOMU_ALLDEPOSITS[indexbis]      =  (*isoPFAllmumap)[themu]/themu->p4().pt(); 
      RECOMU_ChargedDEPOSITS[indexbis]  =  (*isoPFChargedmumap)[themu]/themu->p4().pt();
      RECOMU_NeutralDEPOSITS[indexbis]  =  (*isoPFNeutralmumap)[themu]/themu->p4().pt();
      RECOMU_GammaDEPOSITS[indexbis]    =  (*isoPFGammamumap)[themu]/themu->p4().pt();
      RECOMU_PUDEPOSITS[indexbis]       =  (*isoPFPUmumap)[themu]/themu->p4().pt();
      RECOMU_PULowDEPOSITS[indexbis]    =  (*isoPFPULowmumap)[themu]/themu->p4().pt();
      RECOMU_PFX_DB[indexbis]           =  RECOMU_ChargedDEPOSITS[indexbis] + max(RECOMU_GammaDEPOSITS[indexbis]+RECOMU_NeutralDEPOSITS[indexbis]-0.5*RECOMU_PUDEPOSITS[indexbis] ,0.);
      
      std::cout << "Particle Flow isolation: muon" 
		<< "  \t RECOMU_ALLDEPOSITS="      << RECOMU_ALLDEPOSITS[indexbis]
		<< "  \t RECOMU_ChargedDEPOSITS="  << RECOMU_ChargedDEPOSITS[indexbis]
		<< "  \t RECOMU_NeutralDEPOSITS="  << RECOMU_NeutralDEPOSITS[indexbis]
		<< "  \t RECOMU_GammaDEPOSITS="    << RECOMU_GammaDEPOSITS[indexbis]
		<< "  \t RECOMU_PUDEPOSITS="       << RECOMU_PUDEPOSITS[indexbis]
		<< "  \t RECOMU_PULowDEPOSITS="    << RECOMU_PULowDEPOSITS[indexbis]
		<< "  \t RECOMU_PFX_DB="           << RECOMU_PFX_DB[indexbis]
		<< std::endl; 
          
      // Vertexing
      RECOMU_SIP[indexbis]=(*vertexmumap)[mutrackrefv];
      RECOMU_IP[indexbis]=(*vertexmumapvalue)[mutrackrefv];
      RECOMU_IPERROR[indexbis]=(*vertexmumaperror)[mutrackrefv];
      
      RECOMU_SIP_KF[indexbis]=(*vertexmumapKF)[mutrackrefv];
      RECOMU_IP_KF[indexbis]=(*vertexmumapvalueKF)[mutrackrefv];
      RECOMU_IPERROR_KF[indexbis]=(*vertexmumaperrorKF)[mutrackrefv];
      

/*       /\* RECOMU_SIP_GD[indexbis]=(*vertexmumapGD)[mutrackrefv]; *\/ */
/* /\*       if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_GDMMMM[indexbis]=(*vertexmumapGDMMMM)[mutrackrefv]; *\/ */
/* /\*       RECOMU_SIP_Std[indexbis]=(*vertexmumapStd)[mutrackrefv]; *\/ */
/* /\*       if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_StdMMMM[indexbis]=(*vertexmumapStdMMMM)[mutrackrefv]; *\/ */
/* /\*       RECOMU_SIP_Kin[indexbis]=(*vertexmumapKin)[mutrackrefv]; *\/ */
/* /\*       if (decaychannel=="4mu" || decaychannel=="2e2mu" ) RECOMU_SIP_KinMMMM[indexbis]=(*vertexmumapKinMMMM)[mutrackrefv]; *\/ */


/*       RECOMU_STIP[indexbis]=(*stipmumap)[mutrackrefv]; */
/*       RECOMU_SLIP[indexbis]=(*slipmumap)[mutrackrefv]; */
/*       RECOMU_TIP[indexbis]=(*stipmumapvalue)[mutrackrefv]; */
/*       RECOMU_LIP[indexbis]=(*slipmumapvalue)[mutrackrefv]; */
/*       RECOMU_TIPERROR[indexbis]=(*stipmumaperror)[mutrackrefv]; */
/*       RECOMU_LIPERROR[indexbis]=(*slipmumaperror)[mutrackrefv]; */
      

/*       std::cout << "--vertexing: muon" */
/* 		<< "  Sign3DIP=" << RECOMU_SIP[indexbis] */
/* 		<< "  3DIP="     << RECOMU_IP[indexbis] */
/* 		<< "  3DIPerr="  << RECOMU_IPERROR[indexbis] */
/* 		<< "  SignTIP="  << RECOMU_STIP[indexbis] */
/* 		<< "  TIP="      << RECOMU_TIP[indexbis] */
/* 		<< "  TIPerr="   << RECOMU_TIPERROR[indexbis] */
/* 		<< "  SignLIP="  << RECOMU_SLIP[indexbis] */
/* 		<< "  LIP="      << RECOMU_LIP[indexbis] */
/* 		<< "  LIPerror=" << RECOMU_LIPERROR[indexbis] */
/* 		<< std::endl; */
      
      // Other properties
      RECOMU_numberOfMatches[indexbis]=cand->numberOfMatches();
      RECOMU_caloCompatibility[indexbis]=cand->caloCompatibility();
      RECOMU_segmentCompatibility[indexbis]=(muon::segmentCompatibility( (*cand)));
      RECOMU_glbmuPromptTight[indexbis]=(muon::isGoodMuon( (*cand),muon::GlobalMuonPromptTight));


      std::cout	<< "--other properties:"
		<< "  n.matches="   << RECOMU_numberOfMatches[indexbis]
		<< "  caloComp="    << RECOMU_caloCompatibility[indexbis]
		<< "  segmentComp=" << RECOMU_segmentCompatibility[indexbis]
	        << "  glbmuPromptTight=" << RECOMU_glbmuPromptTight[indexbis]
		<< std::endl;

      if (match(cand->mass(),cand->pt(),cand->charge(),leptonscands2e2mu_)) RECOMUBEST_2e2mu_MATCHED[indexbis]=1;
      if (match(cand->mass(),cand->pt(),cand->charge(),leptonscands4mu_)) RECOMUBEST_4mu_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_MMMM)) RECOMU_MMMM_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_EEMM)) RECOMU_EEMM_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL0)) RECOMU_LLL0_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL1)) RECOMU_LLL1_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL2)) RECOMU_LLL2_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLL3)) RECOMU_LLL3_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLLss0)) RECOMU_LLLLss0_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLLss1)) RECOMU_LLLLss1_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLLss2)) RECOMU_LLLLss2_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLl0)) RECOMU_LLLl0_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLl1)) RECOMU_LLLl1_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_LLLL)) RECOMU_LLLL_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_Z0)) RECOMU_ZMM_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_Zss0)) RECOMU_ZssMM_MATCHED[indexbis]=1;
      if (useAdditionalRECO==true && match(cand->mass(),cand->pt(),cand->charge(),leptonscands_Zcross)) RECOMU_ZEM_MATCHED[indexbis]=1;
      

      // Track properties
      if(cand->globalTrack().isAvailable()){
	RECOMU_mutrkDxy[indexbis]=cand->globalTrack()->dxy(pVertex);
	RECOMU_mutrkDxyError[indexbis]=cand->globalTrack()->dxyError();
	RECOMU_mutrkDxyB[indexbis]=cand->globalTrack()->dxy(bs.position()) ;
	RECOMU_mutrkDz[indexbis]=cand->globalTrack()->dz(pVertex);
	RECOMU_mutrkDzError[indexbis]=cand->globalTrack()->dzError();
	RECOMU_mutrkDzB[indexbis]=cand->globalTrack()->dz(bs.position());
	RECOMU_mutrkChi2PerNdof[indexbis]=cand->globalTrack()->normalizedChi2();
	RECOMU_mutrkCharge[indexbis]=cand->globalTrack()->charge();
 	RECOMU_mutrkNHits[indexbis]=cand->globalTrack()->numberOfValidHits();
	RECOMU_mutrkNPixHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidPixelHits();
        RECOMU_mutrkNStripHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidStripHits();
	RECOMU_mutrkNMuonHits[indexbis]=cand->globalTrack()->hitPattern().numberOfValidMuonHits();


	std::cout << "--muon track properties: "
		  << "  dxy="       << RECOMU_mutrkDxy[indexbis]
		  << "  dxyError="  << RECOMU_mutrkDxyError[indexbis]
		  << "  dxyB="      << RECOMU_mutrkDxyB[indexbis]
		  << "  dz="        << RECOMU_mutrkDz[indexbis]
		  << "  dzError="   << RECOMU_mutrkDzError[indexbis]
		  << "  dzB="       << RECOMU_mutrkDzB[indexbis]
		  << "  chi2_nodf=" << RECOMU_mutrkChi2PerNdof[indexbis]
		  << "  charge="    << RECOMU_mutrkCharge[indexbis]
		  << "  nhits="     << RECOMU_mutrkNHits[indexbis]
		  << "  nPixhits="   << RECOMU_mutrkNPixHits[indexbis]
	          << "  nStriphits=" << RECOMU_mutrkNStripHits[indexbis]
	          << "  nMuonhits="  << RECOMU_mutrkNMuonHits[indexbis]
		  << std::endl;


	// Tracker muon properties
	// RECOMU_trkmuArbitration[indexbis]=(muon::isGoodMuon( (*cand),muon::TrackerMuonArbitrated));
	RECOMU_trkmuArbitration[indexbis]=(muon::segmentCompatibility((*cand),reco::Muon::SegmentAndTrackArbitration));
	RECOMU_trkmu2DCompatibilityLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityLoose));
	RECOMU_trkmu2DCompatibilityTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TM2DCompatibilityTight));
	RECOMU_trkmuOneStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationLoose));
	RECOMU_trkmuOneStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationTight));
	RECOMU_trkmuLastStationLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationLoose));
	RECOMU_trkmuLastStationTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationTight));
	RECOMU_trkmuOneStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngLoose));
	RECOMU_trkmuOneStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMOneStationAngTight));
	RECOMU_trkmuLastStationAngLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngLoose));
	RECOMU_trkmuLastStationAngTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationAngTight));
	RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtLoose));
	RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]=(muon::isGoodMuon( (*cand),muon::TMLastStationOptimizedLowPtTight));
	
	std::cout << "--tracker muon properties:"
		  << "  arbitration="              << RECOMU_trkmuArbitration[indexbis]
		  << "  2DCompLoose="              << RECOMU_trkmu2DCompatibilityLoose[indexbis]
		  << "  2DCompTight="              << RECOMU_trkmu2DCompatibilityTight[indexbis]
		  << "  1StationLoose="            << RECOMU_trkmuOneStationLoose[indexbis]
		  << "  1StationTight="            << RECOMU_trkmuOneStationTight[indexbis]
		  << "  LastStationLoose="         << RECOMU_trkmuLastStationLoose[indexbis]
		  << "  LastStationTight="         << RECOMU_trkmuLastStationTight[indexbis]
		  << "  1StationAngLoose="         << RECOMU_trkmuOneStationAngLoose[indexbis]
		  << "  1StationAngTight="         << RECOMU_trkmuOneStationAngTight[indexbis]
		  << "  LastStationAngLoose="      << RECOMU_trkmuLastStationAngLoose[indexbis]
		  << "  LastStationAngTight="      << RECOMU_trkmuLastStationAngTight[indexbis]
		  << "  LastStationOptLowptLoose=" << RECOMU_trkmuLastStationOptimizedLowPtLoose[indexbis]
		  << "  LastStationOptLowptTight=" << RECOMU_trkmuLastStationOptimizedLowPtTight[indexbis]
		  << std::endl;
      }
      
      indexbis++;
    } 
  }
  
  double SetMET(const edm::Event& iEvent, edm::InputTag myTag_ ){
    double met=-999.;
    edm::Handle<reco::METCollection> metHandle;
    iEvent.getByLabel(myTag_,metHandle);
    for ( METCollection::const_iterator iMet=metHandle->begin(); iMet!=metHandle->end(); iMet++) {
      met = iMet->pt();
    }
    return met;
  }
  
  double SetCaloMET(const edm::Event& iEvent, edm::InputTag myTag_ ){
    double met=-999.;
    edm::Handle<reco::CaloMETCollection> metHandle;
    iEvent.getByLabel(myTag_,metHandle);
    for ( CaloMETCollection::const_iterator iMet=metHandle->begin(); iMet!=metHandle->end(); iMet++) {
      met = iMet->pt();
    }
    return met;
  }
  

  void fillMET(const edm::Event& iEvent){
    //TCMET
    tcmet = SetMET(iEvent,trackermetTag_);
    //CALOMET
    calomet = SetCaloMET(iEvent,calometTag_);
    calometnohf = SetCaloMET(iEvent,calometnohfTag_);

    if (useAdditionalMET_==true){
      calometho = SetCaloMET(iEvent,calomethoTag_);
      calometopt = SetCaloMET(iEvent,calometoptTag_);
      calometoptnohf =  SetCaloMET(iEvent,calometoptnohfTag_);
      calometoptnohfho =  SetCaloMET(iEvent,calometoptnohfhoTag_);
      calometoptho = SetCaloMET(iEvent,calomethoTag_);
      calometnohfho = SetCaloMET(iEvent,calometnohfhoTag_);
    }


    //htMET
    htmetic5 = SetMET(iEvent,htmetic5Tag_);
    htmetkt4 = SetMET(iEvent,htmetkt4Tag_);
    htmetkt6 = SetMET(iEvent,htmetkt6Tag_);
    htmetsc5 = SetMET(iEvent,htmetsc5Tag_);
    htmetsc7 = SetMET(iEvent,htmetsc7Tag_);
    //JES Correction MET
    /* if (fillMCTruth){ */
    /*       jescormetic5 = SetCaloMET(iEvent,jescormetic5Tag_); */
    /*       jescormetkt4 = SetCaloMET(iEvent,jescormetkt4Tag_); */
    /*       jescormetkt6 = SetCaloMET(iEvent,jescormetkt6Tag_); */
    /*       jescormetsc5 = SetCaloMET(iEvent,jescormetsc5Tag_); */
    /*       jescormetsc7 = SetCaloMET(iEvent,jescormetsc7Tag_); */
    /*     } */
    //Type I Muon Correction MET
    cormetmuons = SetCaloMET(iEvent,cormetMuTag_);
    
    //PFMET
    edm::Handle<reco::PFMETCollection> pfmetHandle;
    iEvent.getByLabel(pfmetTag_,pfmetHandle);
    for ( reco::PFMETCollection::const_iterator i=pfmetHandle->begin(); i!=pfmetHandle->end(); i++) {
      pfmet       = i->pt();     
      pfmet_x     = i->px();
      pfmet_y     = i->py();
      pfmet_phi   = i->phi();
      pfmet_theta = i->theta();
    }

    std::cout << "MET:"
	      << "  tcMET="            << tcmet 
	      << "  caloMET="          << calomet 
	      << "  caloMETopt="       << calometopt
	      << "  caloMEToptnohf="   << calometoptnohf
	      << "  caloMEToptnohfho=" << calometoptnohfho
	      << "  caloMEToptho="     << calometoptho
	      << "  caloMETnohf="      << calometnohf
	      << "  caloMETnohfho="    << calometnohfho
	      << "  caloMETho="        << calometho
	      << "  TypeIMucorrMET="   << cormetmuons
	      << "  pfMET="            << pfmet
	      << std::endl;
    
  }

  void fillGD2e2mu(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD;
    iEvent.getByLabel(ftsigma_Vert, GeomD);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD->begin(); cand!=GeomD->end(); ++cand){
      ftsigma[jjj]=(*cand);
      jjj++;
    }

   
    edm::Handle<vector<double> > GeomDlag;
    iEvent.getByLabel(ftsigmalag_Vert, GeomDlag);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag->begin(); cand!=GeomDlag->end(); ++cand){
      ftsigmalag[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdX_;
    iEvent.getByLabel(gdX_Vert, gdX_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdX_->begin(); cand!=gdX_->end(); ++cand){
      gdX[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagX_;
    iEvent.getByLabel(gdlagX_Vert, gdlagX_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagX_->begin(); cand!=gdlagX_->end(); ++cand){
      gdlagX[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdY_;
    iEvent.getByLabel(gdY_Vert, gdY_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdY_->begin(); cand!=gdY_->end(); ++cand){
      gdY[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagY_;
    iEvent.getByLabel(gdlagY_Vert, gdlagY_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagY_->begin(); cand!=gdlagY_->end(); ++cand){
      gdlagY[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZ_;
    iEvent.getByLabel(gdZ_Vert, gdZ_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZ_->begin(); cand!=gdZ_->end(); ++cand){
      gdZ[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZ_;
    iEvent.getByLabel(gdlagZ_Vert, gdlagZ_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZ_->begin(); cand!=gdlagZ_->end(); ++cand){
      gdlagZ[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProb_;
    iEvent.getByLabel(gdlagProb_Vert, gdlagProb_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProb_->begin(); cand!=gdlagProb_->end(); ++cand){
      gdlagProb[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdof_;
    iEvent.getByLabel(gdlagNdof_Vert, gdlagNdof_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdof_->begin(); cand!=gdlagNdof_->end(); ++cand){
      gdlagNdof[jjj]=(*cand);
      jjj++;
    }

  }

  void fillGD4mu(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD_MMMM;
    iEvent.getByLabel(ftsigma_VertMMMM, GeomD_MMMM);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD_MMMM->begin(); cand!=GeomD_MMMM->end(); ++cand){
      ftsigmaMMMM[jjj]=(*cand);
      jjj++;
    }
  
    edm::Handle<vector<double> > GeomDlag_MMMM;
    iEvent.getByLabel(ftsigmalag_VertMMMM, GeomDlag_MMMM);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag_MMMM->begin(); cand!=GeomDlag_MMMM->end(); ++cand){
      ftsigmalagMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdXMMMM_;
    iEvent.getByLabel(gdX_VertMMMM, gdXMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdXMMMM_->begin(); cand!=gdXMMMM_->end(); ++cand){
      gdXMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagXMMMM_;
    iEvent.getByLabel(gdlagX_VertMMMM, gdlagXMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagXMMMM_->begin(); cand!=gdlagXMMMM_->end(); ++cand){
      gdlagXMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdYMMMM_;
    iEvent.getByLabel(gdY_VertMMMM, gdYMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdYMMMM_->begin(); cand!=gdYMMMM_->end(); ++cand){
      gdYMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagYMMMM_;
    iEvent.getByLabel(gdlagY_VertMMMM, gdlagYMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagYMMMM_->begin(); cand!=gdlagYMMMM_->end(); ++cand){
      gdlagYMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZMMMM_;
    iEvent.getByLabel(gdZ_VertMMMM, gdZMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZMMMM_->begin(); cand!=gdZMMMM_->end(); ++cand){
      gdZMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZMMMM_;
    iEvent.getByLabel(gdlagZ_VertMMMM, gdlagZMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZMMMM_->begin(); cand!=gdlagZMMMM_->end(); ++cand){
      gdlagZMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProbMMMM_;
    iEvent.getByLabel(gdlagProb_VertMMMM, gdlagProbMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProbMMMM_->begin(); cand!=gdlagProbMMMM_->end(); ++cand){
      gdlagProbMMMM[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdofMMMM_;
    iEvent.getByLabel(gdlagNdof_VertMMMM, gdlagNdofMMMM_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdofMMMM_->begin(); cand!=gdlagNdofMMMM_->end(); ++cand){
      gdlagNdofMMMM[jjj]=(*cand);
      jjj++;
    }

  }

  void fillGD4e(const edm::Event& iEvent){
    edm::Handle<vector<double> > GeomD_EEEE;
    iEvent.getByLabel(ftsigma_VertEEEE, GeomD_EEEE);
    int jjj=0;
    for (vector<double>::const_iterator cand=GeomD_EEEE->begin(); cand!=GeomD_EEEE->end(); ++cand){
      ftsigmaEEEE[jjj]=(*cand);
      jjj++;
    }
  
    edm::Handle<vector<double> > GeomDlag_EEEE;
    iEvent.getByLabel(ftsigmalag_VertEEEE, GeomDlag_EEEE);
    jjj=0;
    for (vector<double>::const_iterator cand=GeomDlag_EEEE->begin(); cand!=GeomDlag_EEEE->end(); ++cand){
      ftsigmalagEEEE[jjj]=(*cand);
      jjj++;
    }


    edm::Handle<vector<double> > gdXEEEE_;
    iEvent.getByLabel(gdX_VertEEEE, gdXEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdXEEEE_->begin(); cand!=gdXEEEE_->end(); ++cand){
      gdXEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagXEEEE_;
    iEvent.getByLabel(gdlagX_VertEEEE, gdlagXEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagXEEEE_->begin(); cand!=gdlagXEEEE_->end(); ++cand){
      gdlagXEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdYEEEE_;
    iEvent.getByLabel(gdY_VertEEEE, gdYEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdYEEEE_->begin(); cand!=gdYEEEE_->end(); ++cand){
      gdYEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagYEEEE_;
    iEvent.getByLabel(gdlagY_VertEEEE, gdlagYEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagYEEEE_->begin(); cand!=gdlagYEEEE_->end(); ++cand){
      gdlagYEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdZEEEE_;
    iEvent.getByLabel(gdZ_VertEEEE, gdZEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdZEEEE_->begin(); cand!=gdZEEEE_->end(); ++cand){
      gdZEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagZEEEE_;
    iEvent.getByLabel(gdlagZ_VertEEEE, gdlagZEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagZEEEE_->begin(); cand!=gdlagZEEEE_->end(); ++cand){
      gdlagZEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagProbEEEE_;
    iEvent.getByLabel(gdlagProb_VertEEEE, gdlagProbEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagProbEEEE_->begin(); cand!=gdlagProbEEEE_->end(); ++cand){
      gdlagProbEEEE[jjj]=(*cand);
      jjj++;
    }

    edm::Handle<vector<double> > gdlagNdofEEEE_;
    iEvent.getByLabel(gdlagNdof_VertEEEE, gdlagNdofEEEE_);
    jjj=0;
    for (vector<double>::const_iterator cand=gdlagNdofEEEE_->begin(); cand!=gdlagNdofEEEE_->end(); ++cand){
      gdlagNdofEEEE[jjj]=(*cand);
      jjj++;
    }

  }

  void fillConstraintVtx2e2mu(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertex, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      StdFitVertexX[jjj]=cand->position().x();
      StdFitVertexY[jjj]=cand->position().y();
      StdFitVertexZ[jjj]=cand->position().z();
      StdFitVertexChi2r[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProb[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit: " 
	   <<  StdFitVertexX[jjj] << " " 
	   <<  StdFitVertexY[jjj] << " " 
	   <<  StdFitVertexZ[jjj] << " " 
	   <<  StdFitVertexChi2r[jjj] << " " 
	   <<  StdFitVertexProb[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertex, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      KinFitVertexX[jjj]=cand->position().x();
      KinFitVertexY[jjj]=cand->position().y();
      KinFitVertexZ[jjj]=cand->position().z();
      KinFitVertexChi2r[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProb[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Kin Fit: " 
	   <<  KinFitVertexX[jjj] << " " 
	   <<  KinFitVertexY[jjj] << " " 
	   <<  KinFitVertexZ[jjj] << " " 
	   <<  KinFitVertexChi2r[jjj] << " " 
	   <<  KinFitVertexProb[jjj] << endl;
      jjj++;
    }
  }

  void fillConstraintVtx4mu(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertexMMMM, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      StdFitVertexXMMMM[jjj]=cand->position().x();
      StdFitVertexYMMMM[jjj]=cand->position().y();
      StdFitVertexZMMMM[jjj]=cand->position().z();
      StdFitVertexChi2rMMMM[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MMMM: " 
	   <<  StdFitVertexXMMMM[jjj] << " " 
	   <<  StdFitVertexYMMMM[jjj] << " " 
	   <<  StdFitVertexZMMMM[jjj] << " " 
	   <<  StdFitVertexChi2rMMMM[jjj] << " " 
	   <<  StdFitVertexProbMMMM[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertexMMMM, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      KinFitVertexXMMMM[jjj]=cand->position().x();
      KinFitVertexYMMMM[jjj]=cand->position().y();
      KinFitVertexZMMMM[jjj]=cand->position().z();
      KinFitVertexChi2rMMMM[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProbMMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Kin Fit MMMM: " 
	   <<  KinFitVertexXMMMM[jjj] << " " 
	   <<  KinFitVertexYMMMM[jjj] << " " 
	   <<  KinFitVertexZMMMM[jjj] << " " 
	   <<  KinFitVertexChi2rMMMM[jjj] << " " 
	   <<  KinFitVertexProbMMMM[jjj] << endl;
      jjj++;
    }
  }

  void fillConstraintVtx4e(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtx_;
    iEvent.getByLabel(StandardFitVertexEEEE, StandardFitVtx_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtx_->begin(); cand!=StandardFitVtx_->end(); ++cand){
      StdFitVertexXEEEE[jjj]=cand->position().x();
      StdFitVertexYEEEE[jjj]=cand->position().y();
      StdFitVertexZEEEE[jjj]=cand->position().z();
      StdFitVertexChi2rEEEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbEEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit EEEE: " 
	   <<  StdFitVertexXEEEE[jjj] << " " 
	   <<  StdFitVertexYEEEE[jjj] << " " 
	   <<  StdFitVertexZEEEE[jjj] << " " 
	   <<  StdFitVertexChi2rEEEE[jjj] << " " 
	   <<  StdFitVertexProbEEEE[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> KinematicFitVtx_;
    iEvent.getByLabel(KinematicFitVertexEEEE, KinematicFitVtx_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=KinematicFitVtx_->begin(); cand!=KinematicFitVtx_->end(); ++cand){
      KinFitVertexXEEEE[jjj]=cand->position().x();
      KinFitVertexYEEEE[jjj]=cand->position().y();
      KinFitVertexZEEEE[jjj]=cand->position().z();
      KinFitVertexChi2rEEEE[jjj]=cand->chi2()/cand->ndof();
      KinFitVertexProbEEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Kin Fit EEEE: " 
	   <<  KinFitVertexXEEEE[jjj] << " " 
	   <<  KinFitVertexYEEEE[jjj] << " " 
	   <<  KinFitVertexZEEEE[jjj] << " " 
	   <<  KinFitVertexChi2rEEEE[jjj] << " " 
	   <<  KinFitVertexProbEEEE[jjj] << endl;
      jjj++;
    }
  }


  void fillConstraintVtxDiLeptons(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtxDiLep_;
    iEvent.getByLabel(StandardFitVertexDiLep, StandardFitVtxDiLep_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxDiLep_->begin(); cand!=StandardFitVtxDiLep_->end(); ++cand){
      StdFitVertexChi2rDiLep[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbDiLep[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit DiLeptons: " 
	   <<  StdFitVertexChi2rDiLep[jjj] << " " 
	   <<  StdFitVertexProbDiLep[jjj] << endl;
      jjj++;
    }
  }


  void fillConstraintVtxTriLeptons(const edm::Event& iEvent){
    edm::Handle<reco::VertexCollection> StandardFitVtxMMM_;
    iEvent.getByLabel(StandardFitVertexMMM, StandardFitVtxMMM_);
    int jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMMM_->begin(); cand!=StandardFitVtxMMM_->end(); ++cand){
      StdFitVertexChi2rMMM[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMMM[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MMM: " 
	   <<  StdFitVertexChi2rMMM[jjj] << " " 
	   <<  StdFitVertexProbMMM[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxMME_;
    iEvent.getByLabel(StandardFitVertexMME, StandardFitVtxMME_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMME_->begin(); cand!=StandardFitVtxMME_->end(); ++cand){
      StdFitVertexChi2rMME[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMME[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MME: " 
	   <<  StdFitVertexChi2rMME[jjj] << " " 
	   <<  StdFitVertexProbMME[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxEEE_;
    iEvent.getByLabel(StandardFitVertexEEE, StandardFitVtxEEE_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxEEE_->begin(); cand!=StandardFitVtxEEE_->end(); ++cand){
      StdFitVertexChi2rEEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbEEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit EEE: " 
	   <<  StdFitVertexChi2rEEE[jjj] << " " 
	   <<  StdFitVertexProbEEE[jjj] << endl;
      jjj++;
    }

    edm::Handle<reco::VertexCollection> StandardFitVtxMEE_;
    iEvent.getByLabel(StandardFitVertexMEE, StandardFitVtxMEE_);
    jjj=0;
    for (vector<reco::Vertex>::const_iterator cand=StandardFitVtxMEE_->begin(); cand!=StandardFitVtxMEE_->end(); ++cand){
      StdFitVertexChi2rMEE[jjj]=cand->chi2()/cand->ndof();
      StdFitVertexProbMEE[jjj]=TMath::Prob(cand->chi2(),cand->ndof());
      cout << "Std Fit MEE: " 
	   <<  StdFitVertexChi2rMEE[jjj] << " " 
	   <<  StdFitVertexProbMEE[jjj] << endl;
      jjj++;
    }

    
  }


  void fillCP2e2mu(const edm::Event& iEvent){
    edm::Handle<std::vector<float> > cosTheta1Handle;
    iEvent.getByLabel(CP2e2mu_cosTheta1Tag_,cosTheta1Handle);
    
    edm::Handle<std::vector<float> > cosTheta2Handle;
    iEvent.getByLabel(CP2e2mu_cosTheta2Tag_,cosTheta2Handle);
    
    edm::Handle<std::vector<float> > cosThetaStarHandle;
    iEvent.getByLabel(CP2e2mu_cosThetaStarTag_,cosThetaStarHandle);
    
    edm::Handle<std::vector<float> > PhiHandle;
    iEvent.getByLabel(CP2e2mu_PhiTag_,PhiHandle);
    
    edm::Handle<std::vector<float> > Phi1Handle;
    iEvent.getByLabel(CP2e2mu_Phi1Tag_,Phi1Handle);
    
    edm::Handle<std::vector<float> > Phi2Handle;
    iEvent.getByLabel(CP2e2mu_Phi2Tag_,Phi2Handle);
    
    edm::Handle<std::vector<float> > phi1RFHandle;
    iEvent.getByLabel(CP2e2mu_phi1RFTag_,phi1RFHandle);
    
    edm::Handle<std::vector<float> > phi2RFHandle;
    iEvent.getByLabel(CP2e2mu_phi2RFTag_,phi2RFHandle);
    
    
    int iii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta1Handle->begin(); candTheta!=cosTheta1Handle->end(); ++candTheta){
      RECORFBEST_2e2mu_cosTheta1_spin[iii]=(*candTheta);
      iii++;
    }
    
    int ii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta2Handle->begin(); candTheta!=cosTheta2Handle->end(); ++candTheta){
      RECORFBEST_2e2mu_cosTheta2_spin[ii]=(*candTheta);
      ii++;
    }
    
    int iij=0;
    for (std::vector<float>::const_iterator candThetaStar=cosThetaStarHandle->begin(); candThetaStar!=cosThetaStarHandle->end(); ++candThetaStar){
      RECORFBEST_2e2mu_cosThetaStar_spin[iij]=(*candThetaStar);
      iij++;
    }
    
    int iiii=0;
    for (std::vector<float>::const_iterator candPhi=PhiHandle->begin(); candPhi!=PhiHandle->end(); ++candPhi){
      RECORFBEST_2e2mu_Phi_spin[iiii]=(*candPhi);
      iiii++;
    }
    
    int iiiii=0;
    for (std::vector<float>::const_iterator candPhi1=Phi1Handle->begin(); candPhi1!=Phi1Handle->end(); ++candPhi1){
      RECORFBEST_2e2mu_Phi1_spin[iiiii]=(*candPhi1);
      iiiii++;
    }
    
    int iiiiii=0;
    for (std::vector<float>::const_iterator candPhi2=Phi2Handle->begin(); candPhi2!=Phi2Handle->end(); ++candPhi2){
      RECORFBEST_2e2mu_Phi2_spin[iiiiii]=(*candPhi2);
      iiiiii++;
    }
    
    int jjiii=0;
    for (std::vector<float>::const_iterator candPhi1=phi1RFHandle->begin(); candPhi1!=phi1RFHandle->end(); ++candPhi1){
      RECORFBEST_2e2mu_phi1RF_spin[jjiii]=(*candPhi1);
      jjiii++;
    }
    
    int jjiiii=0;
    for (std::vector<float>::const_iterator candPhi2=phi2RFHandle->begin(); candPhi2!=phi2RFHandle->end(); ++candPhi2){
      RECORFBEST_2e2mu_phi2RF_spin[jjiiii]=(*candPhi2);
      jjiiii++;
    }
  }



  void fillCP4e(const edm::Event& iEvent){
    edm::Handle<std::vector<float> > cosTheta1Handle;
    iEvent.getByLabel(CP4e_cosTheta1Tag_,cosTheta1Handle);
    
    edm::Handle<std::vector<float> > cosTheta2Handle;
    iEvent.getByLabel(CP4e_cosTheta2Tag_,cosTheta2Handle);
    
    edm::Handle<std::vector<float> > cosThetaStarHandle;
    iEvent.getByLabel(CP4e_cosThetaStarTag_,cosThetaStarHandle);
    
    edm::Handle<std::vector<float> > PhiHandle;
    iEvent.getByLabel(CP4e_PhiTag_,PhiHandle);
    
    edm::Handle<std::vector<float> > Phi1Handle;
    iEvent.getByLabel(CP4e_Phi1Tag_,Phi1Handle);
    
    edm::Handle<std::vector<float> > Phi2Handle;
    iEvent.getByLabel(CP4e_Phi2Tag_,Phi2Handle);
    
    edm::Handle<std::vector<float> > phi1RFHandle;
    iEvent.getByLabel(CP4e_phi1RFTag_,phi1RFHandle);
    
    edm::Handle<std::vector<float> > phi2RFHandle;
    iEvent.getByLabel(CP4e_phi2RFTag_,phi2RFHandle);
    
    
    int iii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta1Handle->begin(); candTheta!=cosTheta1Handle->end(); ++candTheta){
      RECORFBEST_4e_cosTheta1_spin[iii]=(*candTheta);
      iii++;
    }
    
    int ii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta2Handle->begin(); candTheta!=cosTheta2Handle->end(); ++candTheta){
      RECORFBEST_4e_cosTheta2_spin[ii]=(*candTheta);
      ii++;
    }
    
    int iij=0;
    for (std::vector<float>::const_iterator candThetaStar=cosThetaStarHandle->begin(); candThetaStar!=cosThetaStarHandle->end(); ++candThetaStar){
      RECORFBEST_4e_cosThetaStar_spin[iij]=(*candThetaStar);
      iij++;
    }
    
    int iiii=0;
    for (std::vector<float>::const_iterator candPhi=PhiHandle->begin(); candPhi!=PhiHandle->end(); ++candPhi){
      RECORFBEST_4e_Phi_spin[iiii]=(*candPhi);
      iiii++;
    }
    
    int iiiii=0;
    for (std::vector<float>::const_iterator candPhi1=Phi1Handle->begin(); candPhi1!=Phi1Handle->end(); ++candPhi1){
      RECORFBEST_4e_Phi1_spin[iiiii]=(*candPhi1);
      iiiii++;
    }
    
    int iiiiii=0;
    for (std::vector<float>::const_iterator candPhi2=Phi2Handle->begin(); candPhi2!=Phi2Handle->end(); ++candPhi2){
      RECORFBEST_4e_Phi2_spin[iiiiii]=(*candPhi2);
      iiiiii++;
    }
    
    int jjiii=0;
    for (std::vector<float>::const_iterator candPhi1=phi1RFHandle->begin(); candPhi1!=phi1RFHandle->end(); ++candPhi1){
      RECORFBEST_4e_phi1RF_spin[jjiii]=(*candPhi1);
      jjiii++;
    }
    
    int jjiiii=0;
    for (std::vector<float>::const_iterator candPhi2=phi2RFHandle->begin(); candPhi2!=phi2RFHandle->end(); ++candPhi2){
      RECORFBEST_4e_phi2RF_spin[jjiiii]=(*candPhi2);
      jjiiii++;
    }
  }





 void fillCP4mu(const edm::Event& iEvent){
    edm::Handle<std::vector<float> > cosTheta1Handle;
    iEvent.getByLabel(CP4mu_cosTheta1Tag_,cosTheta1Handle);
    
    edm::Handle<std::vector<float> > cosTheta2Handle;
    iEvent.getByLabel(CP4mu_cosTheta2Tag_,cosTheta2Handle);
    
    edm::Handle<std::vector<float> > cosThetaStarHandle;
    iEvent.getByLabel(CP4mu_cosThetaStarTag_,cosThetaStarHandle);
    
    edm::Handle<std::vector<float> > PhiHandle;
    iEvent.getByLabel(CP4mu_PhiTag_,PhiHandle);
    
    edm::Handle<std::vector<float> > Phi1Handle;
    iEvent.getByLabel(CP4mu_Phi1Tag_,Phi1Handle);
    
    edm::Handle<std::vector<float> > Phi2Handle;
    iEvent.getByLabel(CP4mu_Phi2Tag_,Phi2Handle);
    
    edm::Handle<std::vector<float> > phi1RFHandle;
    iEvent.getByLabel(CP4mu_phi1RFTag_,phi1RFHandle);
    
    edm::Handle<std::vector<float> > phi2RFHandle;
    iEvent.getByLabel(CP4mu_phi2RFTag_,phi2RFHandle);
    
    
    int iii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta1Handle->begin(); candTheta!=cosTheta1Handle->end(); ++candTheta){
      RECORFBEST_4mu_cosTheta1_spin[iii]=(*candTheta);
      iii++;
    }
    
    int ii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta2Handle->begin(); candTheta!=cosTheta2Handle->end(); ++candTheta){
      RECORFBEST_4mu_cosTheta2_spin[ii]=(*candTheta);
      ii++;
    }
    
    int iij=0;
    for (std::vector<float>::const_iterator candThetaStar=cosThetaStarHandle->begin(); candThetaStar!=cosThetaStarHandle->end(); ++candThetaStar){
      RECORFBEST_4mu_cosThetaStar_spin[iij]=(*candThetaStar);
      iij++;
    }
    
    int iiii=0;
    for (std::vector<float>::const_iterator candPhi=PhiHandle->begin(); candPhi!=PhiHandle->end(); ++candPhi){
      RECORFBEST_4mu_Phi_spin[iiii]=(*candPhi);
      iiii++;
    }
    
    int iiiii=0;
    for (std::vector<float>::const_iterator candPhi1=Phi1Handle->begin(); candPhi1!=Phi1Handle->end(); ++candPhi1){
      RECORFBEST_4mu_Phi1_spin[iiiii]=(*candPhi1);
      iiiii++;
    }
    
    int iiiiii=0;
    for (std::vector<float>::const_iterator candPhi2=Phi2Handle->begin(); candPhi2!=Phi2Handle->end(); ++candPhi2){
      RECORFBEST_4mu_Phi2_spin[iiiiii]=(*candPhi2);
      iiiiii++;
    }
    
    int jjiii=0;
    for (std::vector<float>::const_iterator candPhi1=phi1RFHandle->begin(); candPhi1!=phi1RFHandle->end(); ++candPhi1){
      RECORFBEST_4mu_phi1RF_spin[jjiii]=(*candPhi1);
      jjiii++;
    }
    
    int jjiiii=0;
    for (std::vector<float>::const_iterator candPhi2=phi2RFHandle->begin(); candPhi2!=phi2RFHandle->end(); ++candPhi2){
      RECORFBEST_4mu_phi2RF_spin[jjiiii]=(*candPhi2);
      jjiiii++;
    }
  }




  
  void fillMCCP2e2mu(const edm::Event& iEvent){
    edm::Handle<std::vector<float> > cosTheta1Handle;
    iEvent.getByLabel(MCCP_2e2mu_cosTheta1Tag_,cosTheta1Handle);
    
    edm::Handle<std::vector<float> > cosTheta2Handle;
    iEvent.getByLabel(MCCP_2e2mu_cosTheta2Tag_,cosTheta2Handle);
    
    edm::Handle<std::vector<float> > cosThetaStarHandle;
    iEvent.getByLabel(MCCP_2e2mu_cosThetaStarTag_,cosThetaStarHandle);
    
    edm::Handle<std::vector<float> > PhiHandle;
    iEvent.getByLabel(MCCP_2e2mu_PhiTag_,PhiHandle);
    
    edm::Handle<std::vector<float> > Phi1Handle;
    iEvent.getByLabel(MCCP_2e2mu_Phi1Tag_,Phi1Handle);
    
    edm::Handle<std::vector<float> > Phi2Handle;
    iEvent.getByLabel(MCCP_2e2mu_Phi2Tag_,Phi2Handle);
    
    edm::Handle<std::vector<float> > phi1RFHandle;
    iEvent.getByLabel(MCCP_2e2mu_phi1RFTag_,phi1RFHandle);
    
    edm::Handle<std::vector<float> > phi2RFHandle;
    iEvent.getByLabel(MCCP_2e2mu_phi2RFTag_,phi2RFHandle);
    
    
    int iii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta1Handle->begin(); candTheta!=cosTheta1Handle->end(); ++candTheta){
      MCRF_2e2mu_cosTheta1_spin[iii]=(*candTheta);
      iii++;
    }
    
    int ii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta2Handle->begin(); candTheta!=cosTheta2Handle->end(); ++candTheta){
      MCRF_2e2mu_cosTheta2_spin[ii]=(*candTheta);
      ii++;
    }
    
    int iij=0;
    for (std::vector<float>::const_iterator candThetaStar=cosThetaStarHandle->begin(); candThetaStar!=cosThetaStarHandle->end(); ++candThetaStar){
      MCRF_2e2mu_cosThetaStar_spin[iij]=(*candThetaStar);
      iij++;
    }
    
    int iiii=0;
    for (std::vector<float>::const_iterator candPhi=PhiHandle->begin(); candPhi!=PhiHandle->end(); ++candPhi){
      MCRF_2e2mu_Phi_spin[iiii]=(*candPhi);
      iiii++;
    }
    
    int iiiii=0;
    for (std::vector<float>::const_iterator candPhi1=Phi1Handle->begin(); candPhi1!=Phi1Handle->end(); ++candPhi1){
      MCRF_2e2mu_Phi1_spin[iiiii]=(*candPhi1);
      iiiii++;
    }
    
    int iiiiii=0;
    for (std::vector<float>::const_iterator candPhi2=Phi2Handle->begin(); candPhi2!=Phi2Handle->end(); ++candPhi2){
      MCRF_2e2mu_Phi2_spin[iiiiii]=(*candPhi2);
      iiiiii++;
    }
    
    int jjiii=0;
    for (std::vector<float>::const_iterator candPhi1=phi1RFHandle->begin(); candPhi1!=phi1RFHandle->end(); ++candPhi1){
      MCRF_2e2mu_phi1RF_spin[jjiii]=(*candPhi1);
      jjiii++;
    }
    
    int jjiiii=0;
    for (std::vector<float>::const_iterator candPhi2=phi2RFHandle->begin(); candPhi2!=phi2RFHandle->end(); ++candPhi2){
      MCRF_2e2mu_phi2RF_spin[jjiiii]=(*candPhi2);
      jjiiii++;
    }
  }
  
  
  void fillMCCP4mu(const edm::Event& iEvent){
    edm::Handle<std::vector<float> > cosTheta1Handle;
    iEvent.getByLabel(MCCP_4mu_cosTheta1Tag_,cosTheta1Handle);
    
    edm::Handle<std::vector<float> > cosTheta2Handle;
    iEvent.getByLabel(MCCP_4mu_cosTheta2Tag_,cosTheta2Handle);
    
    edm::Handle<std::vector<float> > cosThetaStarHandle;
    iEvent.getByLabel(MCCP_4mu_cosThetaStarTag_,cosThetaStarHandle);
    
    edm::Handle<std::vector<float> > PhiHandle;
    iEvent.getByLabel(MCCP_4mu_PhiTag_,PhiHandle);
    
    edm::Handle<std::vector<float> > Phi1Handle;
    iEvent.getByLabel(MCCP_4mu_Phi1Tag_,Phi1Handle);
    
    edm::Handle<std::vector<float> > Phi2Handle;
    iEvent.getByLabel(MCCP_4mu_Phi2Tag_,Phi2Handle);
    
    edm::Handle<std::vector<float> > phi1RFHandle;
    iEvent.getByLabel(MCCP_4mu_phi1RFTag_,phi1RFHandle);
    
    edm::Handle<std::vector<float> > phi2RFHandle;
    iEvent.getByLabel(MCCP_4mu_phi2RFTag_,phi2RFHandle);
    
    
    int iii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta1Handle->begin(); candTheta!=cosTheta1Handle->end(); ++candTheta){
      MCRF_4mu_cosTheta1_spin[iii]=(*candTheta);
      iii++;
    }
    
    int ii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta2Handle->begin(); candTheta!=cosTheta2Handle->end(); ++candTheta){
      MCRF_4mu_cosTheta2_spin[ii]=(*candTheta);
      ii++;
    }
    
    int iij=0;
    for (std::vector<float>::const_iterator candThetaStar=cosThetaStarHandle->begin(); candThetaStar!=cosThetaStarHandle->end(); ++candThetaStar){
      MCRF_4mu_cosThetaStar_spin[iij]=(*candThetaStar);
      iij++;
    }
    
    int iiii=0;
    for (std::vector<float>::const_iterator candPhi=PhiHandle->begin(); candPhi!=PhiHandle->end(); ++candPhi){
      MCRF_4mu_Phi_spin[iiii]=(*candPhi);
      iiii++;
    }
    
    int iiiii=0;
    for (std::vector<float>::const_iterator candPhi1=Phi1Handle->begin(); candPhi1!=Phi1Handle->end(); ++candPhi1){
      MCRF_4mu_Phi1_spin[iiiii]=(*candPhi1);
      iiiii++;
    }
    
    int iiiiii=0;
    for (std::vector<float>::const_iterator candPhi2=Phi2Handle->begin(); candPhi2!=Phi2Handle->end(); ++candPhi2){
      MCRF_4mu_Phi2_spin[iiiiii]=(*candPhi2);
      iiiiii++;
    }
    
    int jjiii=0;
    for (std::vector<float>::const_iterator candPhi1=phi1RFHandle->begin(); candPhi1!=phi1RFHandle->end(); ++candPhi1){
      MCRF_4mu_phi1RF_spin[jjiii]=(*candPhi1);
      jjiii++;
    }
    
    int jjiiii=0;
    for (std::vector<float>::const_iterator candPhi2=phi2RFHandle->begin(); candPhi2!=phi2RFHandle->end(); ++candPhi2){
      MCRF_4mu_phi2RF_spin[jjiiii]=(*candPhi2);
      jjiiii++;
    }
  }


  



  void fillMCCP4e(const edm::Event& iEvent){
    edm::Handle<std::vector<float> > cosTheta1Handle;
    iEvent.getByLabel(MCCP_4e_cosTheta1Tag_,cosTheta1Handle);
    
    edm::Handle<std::vector<float> > cosTheta2Handle;
    iEvent.getByLabel(MCCP_4e_cosTheta2Tag_,cosTheta2Handle);
    
    edm::Handle<std::vector<float> > cosThetaStarHandle;
    iEvent.getByLabel(MCCP_4e_cosThetaStarTag_,cosThetaStarHandle);
    
    edm::Handle<std::vector<float> > PhiHandle;
    iEvent.getByLabel(MCCP_4e_PhiTag_,PhiHandle);
    
    edm::Handle<std::vector<float> > Phi1Handle;
    iEvent.getByLabel(MCCP_4e_Phi1Tag_,Phi1Handle);
    
    edm::Handle<std::vector<float> > Phi2Handle;
    iEvent.getByLabel(MCCP_4e_Phi2Tag_,Phi2Handle);
    
    edm::Handle<std::vector<float> > phi1RFHandle;
    iEvent.getByLabel(MCCP_4e_phi1RFTag_,phi1RFHandle);
    
    edm::Handle<std::vector<float> > phi2RFHandle;
    iEvent.getByLabel(MCCP_4e_phi2RFTag_,phi2RFHandle);
    
    
    int iii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta1Handle->begin(); candTheta!=cosTheta1Handle->end(); ++candTheta){
      MCRF_4e_cosTheta1_spin[iii]=(*candTheta);
      iii++;
    }
    
    int ii=0;
    for (std::vector<float>::const_iterator candTheta=cosTheta2Handle->begin(); candTheta!=cosTheta2Handle->end(); ++candTheta){
      MCRF_4e_cosTheta2_spin[ii]=(*candTheta);
      ii++;
    }
    
    int iij=0;
    for (std::vector<float>::const_iterator candThetaStar=cosThetaStarHandle->begin(); candThetaStar!=cosThetaStarHandle->end(); ++candThetaStar){
      MCRF_4e_cosThetaStar_spin[iij]=(*candThetaStar);
      iij++;
    }
    
    int iiii=0;
    for (std::vector<float>::const_iterator candPhi=PhiHandle->begin(); candPhi!=PhiHandle->end(); ++candPhi){
      MCRF_4e_Phi_spin[iiii]=(*candPhi);
      iiii++;
    }
    
    int iiiii=0;
    for (std::vector<float>::const_iterator candPhi1=Phi1Handle->begin(); candPhi1!=Phi1Handle->end(); ++candPhi1){
      MCRF_4e_Phi1_spin[iiiii]=(*candPhi1);
      iiiii++;
    }
    
    int iiiiii=0;
    for (std::vector<float>::const_iterator candPhi2=Phi2Handle->begin(); candPhi2!=Phi2Handle->end(); ++candPhi2){
      MCRF_4e_Phi2_spin[iiiiii]=(*candPhi2);
      iiiiii++;
    }
    
    int jjiii=0;
    for (std::vector<float>::const_iterator candPhi1=phi1RFHandle->begin(); candPhi1!=phi1RFHandle->end(); ++candPhi1){
      MCRF_4e_phi1RF_spin[jjiii]=(*candPhi1);
      jjiii++;
    }
    
    int jjiiii=0;
    for (std::vector<float>::const_iterator candPhi2=phi2RFHandle->begin(); candPhi2!=phi2RFHandle->end(); ++candPhi2){
      MCRF_4e_phi2RF_spin[jjiiii]=(*candPhi2);
      jjiiii++;
    }
  }


  void filljets(const edm::Event& iEvent){
    edm::Handle<reco::PFJetCollection> pfjets;
    iEvent.getByLabel(jetsTag_, pfjets);
    RECO_PFJET_N = pfjets->size();
    cout << "Number of CaloJets in the event= " << RECO_PFJET_N << endl;
    
    int index_jets = 0;
    
    for ( PFJetCollection::const_iterator i=pfjets->begin(); i!=pfjets->end(); i++) {  
      if (index_jets>=99) continue;
      
      RECO_PFJET_CHARGE[index_jets] = i->charge();
      RECO_PFJET_ET[index_jets]     = i->et();
      RECO_PFJET_PT[index_jets]     = i->pt();
      RECO_PFJET_ETA[index_jets]    = i->eta();
      RECO_PFJET_PHI[index_jets]    = i->phi();
      
      index_jets++;
    } // for loop on PFJets jets
  
    if(index_jets>=99) { RECO_PFJET_N = 100; cout << "Number of pfjets>100, RECO_PFJETS_N set to 100" << endl;}

    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(rhojetsTag_, rhoHandle);
    RHO=*rhoHandle;
    
  }

  void fillBeamSpot(const edm::Event& iEvent){
    // Beamspot 
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByType(recoBeamSpotHandle);
    const BeamSpot bs = *recoBeamSpotHandle;
    
    BeamSpot_X=bs.position().x();
    BeamSpot_Y=bs.position().y();
    BeamSpot_Z=bs.position().z();

    cout << "BeamSpot:"
	 << "  bs_X=" << BeamSpot_X
	 << "  bs_Y=" << BeamSpot_Y
	 << "  bs_Z=" << BeamSpot_Z
	 << endl;
  }


  void fillBTagging(const edm::Event& iEvent){

    // trackCountingHighEffBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle;
    iEvent.getByLabel(tCHighEff_bTag_,bTagHandle);
    int l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle->begin(); btagIter!=bTagHandle->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighEffBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	tCHighEff_BTagJet_PT[l]=btagIter->first->pt();
	tCHighEff_BTagJet_ETA[l]=btagIter->first->eta();
	tCHighEff_BTagJet_PHI[l]=btagIter->first->phi();
	tCHighEff_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }

    // trackCountingHighPurBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_b;
    iEvent.getByLabel(tCHighPur_bTag_,bTagHandle_b);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_b->begin(); btagIter!=bTagHandle_b->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator trackCountingHighPurBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	tCHighPur_BTagJet_PT[l]=btagIter->first->pt();
	tCHighPur_BTagJet_ETA[l]=btagIter->first->eta();
	tCHighPur_BTagJet_PHI[l]=btagIter->first->phi();
	tCHighPur_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }


    // jetProbabilityHighEffBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_c;
    iEvent.getByLabel(jPHighEff_bTag_,bTagHandle_c);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_c->begin(); btagIter!=bTagHandle_c->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator jetProbabilityHighEffBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;	
	jPHighEff_BTagJet_PT[l]=btagIter->first->pt();
	jPHighEff_BTagJet_ETA[l]=btagIter->first->eta();
	jPHighEff_BTagJet_PHI[l]=btagIter->first->phi();
	jPHighEff_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }

    // jetBProbabilityBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_d;
    iEvent.getByLabel(jBP_bTag_,bTagHandle_d);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_d->begin(); btagIter!=bTagHandle_d->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator jetBProbabilityBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	jBP_BTagJet_PT[l]=btagIter->first->pt();
	jBP_BTagJet_ETA[l]=btagIter->first->eta();
	jBP_BTagJet_PHI[l]=btagIter->first->phi();
	jBP_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }


    // simpleSecondaryVertexHighEffBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_e;
    iEvent.getByLabel(sSVHighEff_bTag_,bTagHandle_e);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_e->begin(); btagIter!=bTagHandle_e->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator simpleSecondaryVertexHighEffBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	sSVHighEff_BTagJet_PT[l]=btagIter->first->pt();
	sSVHighEff_BTagJet_ETA[l]=btagIter->first->eta();
	sSVHighEff_BTagJet_PHI[l]=btagIter->first->phi();
	sSVHighEff_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }

    // simpleSecondaryVertexHighPurBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_f;
    iEvent.getByLabel(sSVHighPur_bTag_,bTagHandle_f);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_f->begin(); btagIter!=bTagHandle_f->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator simpleSecondaryVertexHighPurBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	sSVHighPur_BTagJet_PT[l]=btagIter->first->pt();
	sSVHighPur_BTagJet_ETA[l]=btagIter->first->eta();
	sSVHighPur_BTagJet_PHI[l]=btagIter->first->phi();
	sSVHighPur_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }


    
    // combinedSecondaryVertexBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_g;
    iEvent.getByLabel(cSV_bTag_,bTagHandle_g);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_g->begin(); btagIter!=bTagHandle_g->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator combinedSecondaryVertexBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	cSV_BTagJet_PT[l]=btagIter->first->pt();
	cSV_BTagJet_ETA[l]=btagIter->first->eta();
	cSV_BTagJet_PHI[l]=btagIter->first->phi();
	cSV_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }

    // combinedSecondaryVertexMVABJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_h;
    iEvent.getByLabel(cSVMVA_bTag_,bTagHandle_h);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_h->begin(); btagIter!=bTagHandle_h->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator combinedSecondaryVertexMVABJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	cSVMVA_BTagJet_PT[l]=btagIter->first->pt();
	cSVMVA_BTagJet_ETA[l]=btagIter->first->eta();
	cSVMVA_BTagJet_PHI[l]=btagIter->first->phi();
	cSVMVA_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }


   // softElectronByIP3dBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_i;
    iEvent.getByLabel(sEByIP3d_bTag_,bTagHandle_i);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_i->begin(); btagIter!=bTagHandle_i->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator softElectronByIP3dBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;	
	sEByIP3d_BTagJet_PT[l]=btagIter->first->pt();
	sEByIP3d_BTagJet_ETA[l]=btagIter->first->eta();
	sEByIP3d_BTagJet_PHI[l]=btagIter->first->phi();
	sEByIP3d_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }

    // softElectronByPtBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_l;
    iEvent.getByLabel(sEByPt_bTag_,bTagHandle_l);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_l->begin(); btagIter!=bTagHandle_l->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator softElectronByPtBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;
	sEByPt_BTagJet_PT[l]=btagIter->first->pt();
	sEByPt_BTagJet_ETA[l]=btagIter->first->eta();
	sEByPt_BTagJet_PHI[l]=btagIter->first->phi();
	sEByPt_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }


    // softMuonByIP3dBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_m;
    iEvent.getByLabel(sMByIP3d_bTag_,bTagHandle_m);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_m->begin(); btagIter!=bTagHandle_m->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator softMuonByIP3dBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;
	sMByIP3d_BTagJet_PT[l]=btagIter->first->pt();
	sMByIP3d_BTagJet_ETA[l]=btagIter->first->eta();
	sMByIP3d_BTagJet_PHI[l]=btagIter->first->phi();
	sMByIP3d_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }

    // softMuonByPtBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_n;
    iEvent.getByLabel(sMByPt_bTag_,bTagHandle_n);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_n->begin(); btagIter!=bTagHandle_n->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator softMuonByPtBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	sMByPt_BTagJet_PT[l]=btagIter->first->pt();
	sMByPt_BTagJet_ETA[l]=btagIter->first->eta();
	sMByPt_BTagJet_PHI[l]=btagIter->first->phi();
	sMByPt_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }


    // softMuonBJetTags
    edm::Handle<reco::JetTagCollection> bTagHandle_o;
    iEvent.getByLabel(sM_bTag_,bTagHandle_o);
    l=0;
    for (reco::JetTagCollection::const_iterator btagIter=bTagHandle_o->begin(); btagIter!=bTagHandle_o->end();++btagIter) {
      if(l>=49) continue;
      if (btagIter->second>-99999. && btagIter->second<99999.){
	cout<<" Jet "<< l
	    <<" has b tag discriminator softMuonBJetTags = "<<btagIter->second
	    << " and jet Pt = "<<btagIter->first->pt()<<endl;      
	sM_BTagJet_PT[l]=btagIter->first->pt();
	sM_BTagJet_ETA[l]=btagIter->first->eta();
	sM_BTagJet_PHI[l]=btagIter->first->phi();
	sM_BTagJet_DISCR[l]=btagIter->second;
      }
      l++;
    }  



  }

  
 private:
  
  // ROOT definition
  TFile *theFile_ ;
  TTree *theTree_ ;
  
  bool useRECOformat;
  
  // Input tags
  std::string decaychannel;
  std::string flaginst;
  std::vector<std::string> flagtags;
  std::string rootFileName;

  // PU
  bool fillPUinfo;
  int num_PU_vertices;
  int PU_BunchCrossing; // bunch crossing for the PU event added
  float PU_zpos[50];      // the true primary vertex position along the z axis for each added interaction
  float PU_sumpT_lowpT[50];     // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > low_cut
  float PU_sumpT_highpT[50];    // the sum of the transverse momentum of the tracks originating from each interaction, where track pT > high_cut
  int  PU_ntrks_lowpT[50];       // the number of tracks originating from each interaction, where track pT > low_cut
  int  PU_ntrks_highpT[50];      // the number of tracks originating from each interaction, where track pT > high_cut
  edm::InputTag PileupSrc_;

  // HLT
  bool fillHLTinfo;
  edm::InputTag HLTInfoFired;
  std::string HLTAnalysisinst;
  std::vector<edm::InputTag> flagHLTnames; 

  edm::InputTag triggerEvent,triggerMatchObject,triggerMatchObject_asym,triggerMatchObjectEle;
  std::string triggerFilter,triggerHLTcollection;
  std::vector<std::string> triggerFilter_asym;
  
  // SkimEarlyData
  std::string SkimEarlyDataAnalysisinst;
  std::vector<edm::InputTag> flagSkimEarlyDatanames; 
 
  // MC truth
  bool fillMCTruth;
  std::vector<edm::InputTag> MCcollName,MCcollNameRestFrame;

  // RECO
  bool useAdditionalRECO;
  std::vector<edm::InputTag> RECOcollNameBest2e2mu,RECOcollNameBest4e,RECOcollNameBest4mu,
    RECOcollNameBestRestFrame2e2mu,RECOcollNameBestRestFrame4mu,RECOcollNameBestRestFrame4e;
  std::vector<edm::InputTag> RECOcollNameZ,RECOcollNameZss,
    //RECOcollNameMMMM,RECOcollNameEEEE,RECOcollNameEEMM,
    RECOcollNameLLLLss,RECOcollNameLLL,RECOcollNameLLLl,RECOcollNameLLLLssos;
  edm::InputTag RECOcollNameLLLL,RECOcollNameDiLep;

  edm::InputTag RECOcollNameMMMM,RECOcollNameEEEE,RECOcollNameEEMM,RECOcollNameMMEE,RECOcollNameEEEM,RECOcollNameMMME,RECOcollNameMMMT,RECOcollNameEEET,RECOcollNameMMTT,RECOcollNameEETT,RECOcollNameMMET,RECOcollNameEEMT,RECOcollNameET,RECOcollNameEE,RECOcollNameME,RECOcollNameMT;
  

  // electron and muon tags
  bool useBestCandidate;
  edm::InputTag BestCandidatesLeptonsTag_;
  edm::InputTag electronTag_,electronEgmTag_,muonTag_,clusterCollectionTag_,gsftrackCollection_;
  edm::InputTag electronMapTag_,muonMapTag_;
  edm::InputTag electronTkMapTag_,electronEgmTkMapTag_,muonTkMapTag_;
  edm::InputTag electronEcalMapTag_,electronEgmEcalMapTag_,muonEcalMapTag_;
  edm::InputTag electronHcalMapTag_,electronEgmHcalMapTag_,muonHcalMapTag_;

  // Particle Flow Isolation
  edm::InputTag muonPFIsoValueAllTag_,muonPFIsoValueChargedTag_,muonPFIsoValueNeutralTag_,muonPFIsoValueGammaTag_,muonPFIsoValuePUTag_,muonPFIsoValuePULowTag_;
  edm::InputTag electronPFIsoValueAllTag_,electronPFIsoValueChargedTag_,electronPFIsoValueNeutralTag_,electronPFIsoValueGammaTag_,electronPFIsoValuePUTag_,electronPFIsoValuePULowTag_;
  

  // vertexing 3D 
  edm::InputTag electronTag_Vert,muonTag_Vert;
  edm::InputTag electronMapTag_Vert,electronMapTag_VertKF,muonMapTag_Vert,muonMapTag_VertKF;
  edm::InputTag electronMapTag_VertValue,electronMapTag_VertValueKF,muonMapTag_VertValue,muonMapTag_VertValueKF;
  edm::InputTag electronMapTag_VertError,electronMapTag_VertErrorKF,muonMapTag_VertError,muonMapTag_VertErrorKF;

  edm::InputTag electronMapTag_VertGD,muonMapTag_VertGD;
  edm::InputTag electronMapTag_VertGDEEEE,muonMapTag_VertGDMMMM;
  edm::InputTag electronMapTag_VertStd,muonMapTag_VertStd;
  edm::InputTag electronMapTag_VertStdEEEE,muonMapTag_VertStdMMMM;
  edm::InputTag electronMapTag_VertKin,muonMapTag_VertKin;
  edm::InputTag electronMapTag_VertKinEEEE,muonMapTag_VertKinMMMM;

  // vertexing STIP / SLIP
  edm::InputTag electronSTIPMapTag_Vert,muonSTIPMapTag_Vert;
  edm::InputTag electronSLIPMapTag_Vert,muonSLIPMapTag_Vert;
  edm::InputTag electronSTIPMapTag_VertValue,muonSTIPMapTag_VertValue;
  edm::InputTag electronSLIPMapTag_VertValue,muonSLIPMapTag_VertValue;
  edm::InputTag electronSTIPMapTag_VertError,muonSTIPMapTag_VertError;
  edm::InputTag electronSLIPMapTag_VertError,muonSLIPMapTag_VertError;
  // vertexing GD
  edm::InputTag ftsigma_Vert,ftsigmalag_Vert,
                ftsigma_VertMMMM,ftsigmalag_VertMMMM,
                ftsigma_VertEEEE,ftsigmalag_VertEEEE;
  edm::InputTag 
    gdX_Vert,gdY_Vert,gdZ_Vert,
    gdX_VertMMMM,gdY_VertMMMM,gdZ_VertMMMM,
    gdX_VertEEEE,gdY_VertEEEE,gdZ_VertEEEE;
   edm::InputTag 
    gdlagX_Vert,gdlagY_Vert,gdlagZ_Vert,gdlagProb_Vert,gdlagNdof_Vert,
     gdlagX_VertMMMM,gdlagY_VertMMMM,gdlagZ_VertMMMM,gdlagProb_VertMMMM,gdlagNdof_VertMMMM,
     gdlagX_VertEEEE,gdlagY_VertEEEE,gdlagZ_VertEEEE,gdlagProb_VertEEEE,gdlagNdof_VertEEEE;
  //

   // ConstraintVtx
   edm::InputTag 
     StandardFitVertex,StandardFitVertexMMMM,StandardFitVertexEEEE,
     KinematicFitVertex,KinematicFitVertexMMMM,KinematicFitVertexEEEE;

   edm::InputTag 
     StandardFitVertexEEE,StandardFitVertexMMM,StandardFitVertexMEE,StandardFitVertexMME,StandardFitVertexDiLep;

  //electronID
  std::vector<edm::InputTag> eleIDTag_;
  edm::InputTag EleID_VeryLooseTag_ ;
  edm::InputTag EleID_LooseTag_  ;
  edm::InputTag EleID_MediumTag_  ;
  edm::InputTag EleID_TightTag_ ;

  edm::InputTag EleID_HZZVeryLooseTag_ ;
  edm::InputTag EleID_HZZLooseTag_  ;
  edm::InputTag EleID_HZZMediumTag_  ;
  edm::InputTag EleID_HZZTightTag_ ;
  
  edm::InputTag allelectronsColl, allmuonsColl;
  edm::InputTag isoVarTagElectronsCal,isoVarTag,isoVarTagElectronsTracker,isoVarTagElectronsX;
  edm::InputTag isoVarTagMuonsHCalIso,isoVarTagMuonsECalIso,
    isoVarTagMuonsCalIso,isoVarTagMuonsTracker,isoVarTagMuonsX;
  
  edm::InputTag theECALIsoDepositLabel;    //EM calorimeter Isolation deposit label
  edm::InputTag theHCALIsoDepositLabel;    //Hadron calorimeter Isolation deposit label
  edm::InputTag theHOCALIsoDepositLabel;   //Outer calorimeter Isolation deposit label
  edm::InputTag theTrackerIsoDepositLabel; //Tracker Isolation deposit label 
  
  // Photon, Tracks, Jets, Vertices
  edm::InputTag photonsTag_,tracksTag_,jetsTag_,rhojetsTag_,verticesTag_;
  // MET
  edm::InputTag genmetTag_,trackermetTag_,pfmetTag_;
  edm::InputTag calometTag_,calometoptTag_,calometoptnohfTag_,calometoptnohfhoTag_;
  edm::InputTag calometopthoTag_,calometnohfTag_,calometnohfhoTag_,calomethoTag_;
  bool useAdditionalMET_;
  edm::InputTag htmetic5Tag_,htmetkt4Tag_,htmetkt6Tag_,htmetsc5Tag_,htmetsc7Tag_;
  edm::InputTag jescormetic5Tag_,jescormetkt4Tag_,jescormetkt6Tag_,jescormetsc5Tag_,jescormetsc7Tag_;
  edm::InputTag cormetMuTag_;
  
  // CP variables tag
  edm::InputTag MCCP_2e2mu_cosTheta1Tag_,MCCP_2e2mu_cosTheta2Tag_,MCCP_2e2mu_cosThetaStarTag_,
    MCCP_2e2mu_PhiTag_,MCCP_2e2mu_Phi1Tag_,MCCP_2e2mu_Phi2Tag_,MCCP_2e2mu_phi1RFTag_,MCCP_2e2mu_phi2RFTag_;
  edm::InputTag MCCP_4mu_cosTheta1Tag_,MCCP_4mu_cosTheta2Tag_,MCCP_4mu_cosThetaStarTag_,
    MCCP_4mu_PhiTag_,MCCP_4mu_Phi1Tag_,MCCP_4mu_Phi2Tag_,MCCP_4mu_phi1RFTag_,MCCP_4mu_phi2RFTag_;
  edm::InputTag MCCP_4e_cosTheta1Tag_,MCCP_4e_cosTheta2Tag_,MCCP_4e_cosThetaStarTag_,
    MCCP_4e_PhiTag_,MCCP_4e_Phi1Tag_,MCCP_4e_Phi2Tag_,MCCP_4e_phi1RFTag_,MCCP_4e_phi2RFTag_;

  edm::InputTag CP2e2mu_cosTheta1Tag_,CP2e2mu_cosTheta2Tag_,CP2e2mu_cosThetaStarTag_,CP2e2mu_PhiTag_,CP2e2mu_Phi1Tag_,CP2e2mu_Phi2Tag_,CP2e2mu_phi1RFTag_,CP2e2mu_phi2RFTag_;
  edm::InputTag CP4mu_cosTheta1Tag_,CP4mu_cosTheta2Tag_,CP4mu_cosThetaStarTag_,CP4mu_PhiTag_,CP4mu_Phi1Tag_,CP4mu_Phi2Tag_,CP4mu_phi1RFTag_,CP4mu_phi2RFTag_;
  edm::InputTag CP4e_cosTheta1Tag_,CP4e_cosTheta2Tag_,CP4e_cosThetaStarTag_,CP4e_PhiTag_,CP4e_Phi1Tag_,CP4e_Phi2Tag_,CP4e_phi1RFTag_,CP4e_phi2RFTag_;
  
  // Conversion
  edm::InputTag ConvMapDistTag_,ConvMapDcotTag_;

  // counters
  int irun, ievt,ils,nevt,flag[10];

  // HLT
  bool boolMuonIso, boolMuonNonIso, bool2MuonNonIso; 
  bool boolElectron, boolElectronRelaxed, bool2Electron, bool2ElectronRelaxed, boolCandHLT2ElectronStartup;
  bool boolXElectronMuon, boolHLTMinBiasBSC,boolHLTaccept;
  bool boolHLTpattern[50];
  int RECO_nMuHLTMatch,RECO_nMuHLTMatchPAT,RECO_nMuHLTMatch_asym_PAT,RECO_nEleHLTMatchPAT;
  float RECOMU_PT_MuHLTMatch[100],RECOMU_PT_MuHLTMatchPAT[100],RECOMU_PT_MuHLTMatch_asym_PAT[100],RECOELE_PT_EleHLTMatchPAT[100];
  char HLTPathsFired[20000];

  // Skim Early Data
  bool useSkimEarlyData;
  bool boolSkim_highEnergyMuons,boolSkim_highEnergyElectrons;
  bool boolSkim_recoWMNfromPf,boolSkim_recoWMNfromTc,boolSkim_recoWENfromPf,boolSkim_recoWENfromTc;
  bool boolSkim_diMuonsJPsi,boolSkim_diMuonsZ,boolSkim_diElectronsZ;
  bool boolSkim_triLeptonsMuMuMu,boolSkim_triLeptonsMuMuEl,boolSkim_triLeptonsMuElEl,boolSkim_triLeptonsElElEl;
  bool boolSkim_quadLeptons4Mu,boolSkim_quadLeptons2Mu2El,boolSkim_quadLeptons4El;
  
  // tmp candidate collections
  reco::CandidateCollection *leptonscands2e2mu_;
  reco::CandidateCollection *leptonscands2e2murf_;
  reco::CandidateCollection *leptonscands4mu_;
  reco::CandidateCollection *leptonscands4murf_;
  reco::CandidateCollection *leptonscands4e_;
  reco::CandidateCollection *leptonscands4erf_;
  
  reco::CandidateCollection *leptonscands_Z0;
  reco::CandidateCollection *leptonscands_Z1;
  reco::CandidateCollection *leptonscands_Zss0;
  reco::CandidateCollection *leptonscands_Zss1;
  reco::CandidateCollection *leptonscands_Zcross;
  reco::CandidateCollection *leptonscands_DiLep;
  //reco::CandidateCollection *leptonscands_MMMM;
  //reco::CandidateCollection *leptonscands_EEEE;
  //reco::CandidateCollection *leptonscands_EEMM;
  reco::CandidateCollection *leptonscands_LLL0;
  reco::CandidateCollection *leptonscands_LLL1;
  reco::CandidateCollection *leptonscands_LLL2;
  reco::CandidateCollection *leptonscands_LLL3;
  reco::CandidateCollection *leptonscands_LLLLss0;
  reco::CandidateCollection *leptonscands_LLLLss1;
  reco::CandidateCollection *leptonscands_LLLLss2;
  reco::CandidateCollection *leptonscands_LLLl0;
  reco::CandidateCollection *leptonscands_LLLl1;
  reco::CandidateCollection *leptonscands_LLLL;
  reco::CandidateCollection *leptonscands_MMMM;
  reco::CandidateCollection *leptonscands_EEEE;
  reco::CandidateCollection *leptonscands_EEMM;
  reco::CandidateCollection *leptonscands_MMEE;
  reco::CandidateCollection *leptonscands_MMMT;
  reco::CandidateCollection *leptonscands_EEET;
  reco::CandidateCollection *leptonscands_MMTT;
  reco::CandidateCollection *leptonscands_EETT;
  reco::CandidateCollection *leptonscands_MMET;
  reco::CandidateCollection *leptonscands_EEMT;
  reco::CandidateCollection *leptonscands_MMME;
  reco::CandidateCollection *leptonscands_EEEM;
  reco::CandidateCollection *leptonscands_MT;
  reco::CandidateCollection *leptonscands_ET;
  reco::CandidateCollection *leptonscands_ME;


  // MC info
  double auto_cross_section,external_cross_section,filter_eff;
  double processID,ALPGENid;
  double weight;
  
  // MC truth
  float MC_2e2mu_E[7],MC_2e2mu_PT[7],MC_2e2mu_ETA[7],MC_2e2mu_THETA[7],MC_2e2mu_PHI[7],MC_2e2mu_MASS[7],MC_2e2mu_PDGID[7];
  float MC_4e_E[7],MC_4e_PT[7],MC_4e_ETA[7],MC_4e_THETA[7],MC_4e_PHI[7],MC_4e_MASS[7],MC_4e_PDGID[7];
  float MC_4mu_E[7],MC_4mu_PT[7],MC_4mu_ETA[7],MC_4mu_THETA[7],MC_4mu_PHI[7],MC_4mu_MASS[7],MC_4mu_PDGID[7];

  float MCRF_2e2mu_E[7],MCRF_2e2mu_PT[7],MCRF_2e2mu_ETA[7],MCRF_2e2mu_THETA[7],MCRF_2e2mu_PHI[7],MCRF_2e2mu_MASS[7],MCRF_2e2mu_PDGID[7];
  float MCRF_2e2mu_cosTheta1_spin[4], MCRF_2e2mu_cosTheta2_spin[4], MCRF_2e2mu_cosThetaStar_spin[4], MCRF_2e2mu_Phi_spin[4], 
    MCRF_2e2mu_Phi1_spin[4], MCRF_2e2mu_Phi2_spin[4], MCRF_2e2mu_phi1RF_spin[4], MCRF_2e2mu_phi2RF_spin[4];
  
  float MCRF_4e_E[7],MCRF_4e_PT[7],MCRF_4e_ETA[7],MCRF_4e_THETA[7],MCRF_4e_PHI[7],MCRF_4e_MASS[7],MCRF_4e_PDGID[7];
  float MCRF_4e_cosTheta1_spin[4], MCRF_4e_cosTheta2_spin[4], MCRF_4e_cosThetaStar_spin[4], MCRF_4e_Phi_spin[4], 
    MCRF_4e_Phi1_spin[4], MCRF_4e_Phi2_spin[4], MCRF_4e_phi1RF_spin[4], MCRF_4e_phi2RF_spin[4];
  
  float MCRF_4mu_E[7],MCRF_4mu_PT[7],MCRF_4mu_ETA[7],MCRF_4mu_THETA[7],MCRF_4mu_PHI[7],MCRF_4mu_MASS[7],MCRF_4mu_PDGID[7];
  float MCRF_4mu_cosTheta1_spin[4], MCRF_4mu_cosTheta2_spin[4], MCRF_4mu_cosThetaStar_spin[4], MCRF_4mu_Phi_spin[4], 
    MCRF_4mu_Phi1_spin[4], MCRF_4mu_Phi2_spin[4], MCRF_4mu_phi1RF_spin[4], MCRF_4mu_phi2RF_spin[4];
  
  
  // RECO collection
  float RECOBEST_2e2mu_E[7],RECOBEST_2e2mu_PT[7],RECOBEST_2e2mu_ETA[7],RECOBEST_2e2mu_THETA[7],RECOBEST_2e2mu_PHI[7],RECOBEST_2e2mu_MASS[7],RECOBEST_2e2mu_ORDERED_PT[7];
  bool RECOBEST_2e2mu_isGlobalMu[7],RECOBEST_2e2mu_isStandAloneMu[7],RECOBEST_2e2mu_isTrackerMu[7],RECOBEST_2e2mu_isCaloMu[7],RECOBEST_2e2mu_isElectron[7];
 
  float RECOBEST_4mu_E[7],RECOBEST_4mu_PT[7],RECOBEST_4mu_ETA[7],RECOBEST_4mu_THETA[7],RECOBEST_4mu_PHI[7],RECOBEST_4mu_MASS[7],RECOBEST_4mu_ORDERED_PT[7];
  bool RECOBEST_4mu_isGlobalMu[7],RECOBEST_4mu_isStandAloneMu[7],RECOBEST_4mu_isTrackerMu[7],RECOBEST_4mu_isCaloMu[7];
 
  float RECOBEST_4e_E[7],RECOBEST_4e_PT[7],RECOBEST_4e_ETA[7],RECOBEST_4e_THETA[7],RECOBEST_4e_PHI[7],RECOBEST_4e_MASS[7],RECOBEST_4e_ORDERED_PT[7];
  bool RECOBEST_4e_isGlobalMu[7],RECOBEST_4e_isStandAloneMu[7],RECOBEST_4e_isTrackerMu[7],RECOBEST_4e_isCaloMu[7],RECOBEST_4e_isElectron[7];
 
  
  // RECORF
  bool fillRF;
  float RECORFBEST_2e2mu_E[7],RECORFBEST_2e2mu_PT[7],RECORFBEST_2e2mu_ETA[7],RECORFBEST_2e2mu_THETA[7],RECORFBEST_2e2mu_PHI[7],RECORFBEST_2e2mu_MASS[7];
    
  float RECORFBEST_2e2mu_cosTheta1_spin[4], RECORFBEST_2e2mu_cosTheta2_spin[4], RECORFBEST_2e2mu_cosThetaStar_spin[4], RECORFBEST_2e2mu_Phi_spin[4], 
    RECORFBEST_2e2mu_Phi1_spin[4], RECORFBEST_2e2mu_Phi2_spin[4], RECORFBEST_2e2mu_phi1RF_spin[4], RECORFBEST_2e2mu_phi2RF_spin[4];
  

  float RECORFBEST_4e_E[7],RECORFBEST_4e_PT[7],RECORFBEST_4e_ETA[7],RECORFBEST_4e_THETA[7],RECORFBEST_4e_PHI[7],RECORFBEST_4e_MASS[7];
    
  float RECORFBEST_4e_cosTheta1_spin[4], RECORFBEST_4e_cosTheta2_spin[4], RECORFBEST_4e_cosThetaStar_spin[4], RECORFBEST_4e_Phi_spin[4], 
    RECORFBEST_4e_Phi1_spin[4], RECORFBEST_4e_Phi2_spin[4], RECORFBEST_4e_phi1RF_spin[4], RECORFBEST_4e_phi2RF_spin[4];
  
  float RECORFBEST_4mu_E[7],RECORFBEST_4mu_PT[7],RECORFBEST_4mu_ETA[7],RECORFBEST_4mu_THETA[7],RECORFBEST_4mu_PHI[7],RECORFBEST_4mu_MASS[7];
    
  float RECORFBEST_4mu_cosTheta1_spin[4], RECORFBEST_4mu_cosTheta2_spin[4], RECORFBEST_4mu_cosThetaStar_spin[4], RECORFBEST_4mu_Phi_spin[4], 
    RECORFBEST_4mu_Phi1_spin[4], RECORFBEST_4mu_Phi2_spin[4], RECORFBEST_4mu_phi1RF_spin[4], RECORFBEST_4mu_phi2RF_spin[4];
  

  int leptonflavor;

  // RECO additional
  float 
    RECO_ZMM_MASS[50],RECO_ZMM_PT[3][50],RECO_ZMM_ETA[3][50],RECO_ZMM_PHI[3][50],
    RECO_ZEE_MASS[50],RECO_ZEE_PT[3][50],RECO_ZEE_ETA[3][50],RECO_ZEE_PHI[3][50],
    RECO_ZMMss_MASS[50],RECO_ZMMss_PT[3][50],RECO_ZMMss_ETA[3][50],RECO_ZMMss_PHI[3][50],
    RECO_ZEEss_MASS[50],RECO_ZEEss_PT[3][50],RECO_ZEEss_ETA[3][50],RECO_ZEEss_PHI[3][50],
    RECO_ZEM_MASS[50],RECO_ZEM_PT[3][50],RECO_ZEM_ETA[3][50],RECO_ZEM_PHI[3][50],
    RECO_DiLep_MASS[50],RECO_DiLep_PT[3][50],RECO_DiLep_ETA[3][50],RECO_DiLep_PHI[3][50];
  

  float 
    //RECO_EEMM_MASS[7][50],RECO_MMMM_MASS[7][50],RECO_EEEE_MASS[7][50],
    //RECO_EEMM_PT[7][50],  RECO_MMMM_PT[7][50],  RECO_EEEE_PT[7][50],
    //RECO_EEMM_ETA[7][50],RECO_MMMM_ETA[7][50],RECO_EEEE_ETA[7][50],
    //RECO_EEMM_PHI[7][50],  RECO_MMMM_PHI[7][50],  RECO_EEEE_PHI[7][50],
    RECO_LLLL_MASS[7][100],RECO_LLLL_PT[7][100],RECO_LLLL_ETA[7][100],RECO_LLLL_PHI[7][100];

  float	RECO_LLL0_MASS[50],RECO_LLL1_MASS[50],RECO_LLL2_MASS[50],RECO_LLL3_MASS[50];
  float	RECO_LLL0_PT[4][50],RECO_LLL1_PT[4][50],RECO_LLL2_PT[4][50],RECO_LLL3_PT[4][50];

  float	RECO_LLLl0_MASS[50],RECO_LLLl1_MASS[50];
  float	RECO_LLLl0_PT[5][50],RECO_LLLl1_PT[5][50];

  float 
    RECO_LLLL0ss_MASS[50],RECO_LLLL0ss_PT[5][50],
    RECO_LLLL1ss_MASS[50],RECO_LLLL1ss_PT[5][50],
    RECO_LLLL2ss_MASS[50],RECO_LLLL2ss_PT[5][50];

  float 
    RECO_MMMT_MASS[9][100],RECO_MMMT_PT[9][100],RECO_MMMT_ETA[9][100],RECO_MMMT_PHI[9][100],RECO_MMMT_CHARGE[9][100],
    
    RECO_EEET_MASS[9][100],RECO_EEET_PT[9][100],RECO_EEET_ETA[9][100],RECO_EEET_PHI[9][100],RECO_EEET_CHARGE[9][100],
    
    RECO_MMTT_MASS[9][200],RECO_MMTT_PT[9][200],RECO_MMTT_ETA[9][200],RECO_MMTT_PHI[9][200],RECO_MMTT_CHARGE[9][200],
    
    RECO_EETT_MASS[9][200],RECO_EETT_PT[9][200],RECO_EETT_ETA[9][200],RECO_EETT_PHI[9][200],RECO_EETT_CHARGE[9][200],

    RECO_MMET_MASS[9][100],RECO_MMET_PT[9][100],RECO_MMET_ETA[9][100],RECO_MMET_PHI[9][100],RECO_MMET_CHARGE[9][100],

    RECO_EEMT_MASS[9][100],RECO_EEMT_PT[9][100],RECO_EEMT_ETA[9][100],RECO_EEMT_PHI[9][100],RECO_EEMT_CHARGE[9][100],

    RECO_MMME_MASS[9][100],RECO_MMME_PT[9][100],RECO_MMME_ETA[9][100],RECO_MMME_PHI[9][100],RECO_MMME_CHARGE[9][100],

    RECO_EEEM_MASS[9][100],RECO_EEEM_PT[9][100],RECO_EEEM_ETA[9][100],RECO_EEEM_PHI[9][100],RECO_EEEM_CHARGE[9][100],
    
    RECO_MMMM_MASS[9][100],RECO_MMMM_PT[9][100],RECO_MMMM_ETA[9][100],RECO_MMMM_PHI[9][100],RECO_MMMM_CHARGE[9][100],

    RECO_EEEE_MASS[9][100],RECO_EEEE_PT[9][100],RECO_EEEE_ETA[9][100],RECO_EEEE_PHI[9][100],RECO_EEEE_CHARGE[9][100],

    RECO_EEMM_MASS[9][100],RECO_EEMM_PT[9][100],RECO_EEMM_ETA[9][100],RECO_EEMM_PHI[9][100],RECO_EEMM_CHARGE[9][100],

    RECO_MMEE_MASS[9][100],RECO_MMEE_PT[9][100],RECO_MMEE_ETA[9][100],RECO_MMEE_PHI[9][100],RECO_MMEE_CHARGE[9][100],

    RECO_MT_MASS[5][100],RECO_MT_PT[5][100],RECO_MT_ETA[5][100],RECO_MT_PHI[5][100],RECO_MT_CHARGE[5][100],

    RECO_ET_MASS[5][100],RECO_ET_PT[5][100],RECO_ET_ETA[5][100],RECO_ET_PHI[5][100],RECO_ET_CHARGE[5][100],

    RECO_ME_MASS[5][100],RECO_ME_PT[5][100],RECO_ME_ETA[5][100],RECO_ME_PHI[5][100],RECO_ME_CHARGE[5][100];


  //hpstaus
  double HPSTAU_discByLER[100],HPSTAU_discByMER[100],HPSTAU_discByTER[100],HPSTAU_discByLMR[100],HPSTAU_discByTMR[100],HPSTAU_discByDecayF[100],HPSTAU_discLooseIso[100],HPSTAU_discVLooseIso[100],HPSTAU_discMediumIso[100],HPSTAU_discTightIso[100],HPSTAU_discLooseChargedIso[100],HPSTAU_discVLooseChargedIso[100],HPSTAU_discMediumChargedIso[100],HPSTAU_discTightChargedIso[100],HPSTAU_discLooseIsoDB[100],HPSTAU_discVLooseIsoDB[100],HPSTAU_discMediumIsoDB[100],HPSTAU_discTightIsoDB[100],HPSTAU_discLooseCombinedIsoDB[100],HPSTAU_discVLooseCombinedIsoDB[100],HPSTAU_discMediumCombinedIsoDB[100],HPSTAU_discTightCombinedIsoDB[100];
  float  HPSTAU_MASS[100],HPSTAU_PT[100],HPSTAU_ETA[100],HPSTAU_PHI[100],HPSTAU_CHARGE[100];
  
  //ElectronDisc
  double EleId95relIso[100],EleId90relIso[100],EleId80relIso[100],EleId85relIso[100];
  
  
  // RECO electrons
  edm::ESHandle<CaloGeometry> theCaloGeom_;  
  float RECOELE_E[100],RECOELE_PT[100],RECOELE_P[100],RECOELE_ETA[100],RECOELE_THETA[100],RECOELE_PHI[100],RECOELE_MASS[100],RECOELE_TMASS[100],RECOELE_QUALITY[100];
  float RECOELE_CHARGE[100];
  bool RECOELE_isEcalDriven[100], RECOELE_isTrackerDriven[100];
  float 
    RECOELE_gsftrack_chi2[100], 
    RECOELE_gsftrack_dxyB[100], RECOELE_gsftrack_dxy[100], RECOELE_gsftrack_dxyError[100],
    RECOELE_gsftrack_dzB[100], RECOELE_gsftrack_dz[100],RECOELE_gsftrack_dzError[100];
  int RECOELE_gsftrack_losthits[100],RECOELE_gsftrack_validhits[100],RECOELE_gsftrack_expected_inner_hits[100]; 
  float
    RECOELE_scl_E[100],RECOELE_scl_Et[100],RECOELE_scl_Eta[100],RECOELE_scl_Phi[100]; 
  //float RECOELE_eSuperClusterOverP[100], RECOELE_eSeedClusterOverPout[100], RECOELE_deltaEtaSuperClusterTrackAtVtx[100], RECOELE_deltaPhiSuperClusterTrackAtVtx[100];
  float 
    RECOELE_ep[100], RECOELE_eSeedp[100], RECOELE_eSeedpout[100], RECOELE_eElepout[100],
    RECOELE_deltaEtaIn[100],RECOELE_deltaEtaSeed[100],RECOELE_deltaEtaEle[100],RECOELE_deltaPhiIn[100],
    RECOELE_deltaPhiSeed[100],RECOELE_deltaPhiEle[100];
  int RECOELE_isbarrel[100], RECOELE_isendcap[100], RECOELE_isEBetaGap[100], RECOELE_isEBphiGap[100], RECOELE_isEEdeeGap[100], RECOELE_isEEringGap[100]; 
  float RECOELE_sigmaIetaIeta[100], RECOELE_sigmaEtaEta[100], RECOELE_e15[100], RECOELE_e25max[100], RECOELE_e55[100], RECOELE_he[100];
  float RECOELE_mva[100], RECOELE_fbrem[100],RECOELE_fbrem_mean[100],RECOELE_fbrem_mode[100];
  int RECOELE_nbrems[100], RECOELE_Class[100];
  
  float RECOELE_TRACKISO[100],RECOELE_ECALISO[100],RECOELE_HCALISO[100], RECOELE_X[100],RECOELE_EGMTRACKISO[100],RECOELE_EGMECALISO[100],RECOELE_EGMHCALISO[100], RECOELE_EGMX[100],
    RECOELE_IP[100],RECOELE_SIP[100],RECOELE_IPERROR[100],
    RECOELE_IP_KF[100],RECOELE_SIP_KF[100],RECOELE_IPERROR_KF[100],
    RECOELE_STIP[100],RECOELE_TIP[100],RECOELE_TIPERROR[100],
    RECOELE_SLIP[100],RECOELE_LIP[100],RECOELE_LIPERROR[100],
    RECOELE_SIP_GD[100], RECOELE_SIP_GDEEEE[100],
    RECOELE_SIP_Std[100], RECOELE_SIP_StdEEEE[100], 
    RECOELE_SIP_Kin[100], RECOELE_SIP_KinEEEE[100];
  int RECOELEBEST_2e2mu_MATCHED[100],RECOELEBEST_4e_MATCHED[100];
  int RECOELE_EEEE_MATCHED[100],RECOELE_EEMM_MATCHED[100],RECOELE_ZEE_MATCHED[100],RECOELE_ZssEE_MATCHED[100],RECOELE_ZEM_MATCHED[100],
    RECOELE_LLL0_MATCHED[100],RECOELE_LLL1_MATCHED[100],RECOELE_LLL2_MATCHED[100],RECOELE_LLL3_MATCHED[100],
    RECOELE_LLLLss0_MATCHED[100],RECOELE_LLLLss1_MATCHED[100],RECOELE_LLLLss2_MATCHED[100],
    RECOELE_LLLl0_MATCHED[100],RECOELE_LLLl1_MATCHED[100],RECOELE_LLLL_MATCHED[100];
  bool RECOELE_eleID[100][2];
  
  double ele_sclRawE[100] ;
  double ele_sclX[100], ele_sclY[100], ele_sclZ[100];
  int ele_seedSubdet1[100];
  double ele_seedDphi1[100], ele_seedDrz1[100];
  int ele_seedSubdet2[100];
  double ele_seedDphi2[100], ele_seedDrz2[100];
  int ele_severityLevelSeed[100], ele_severityLevelClusters[100], ele_outOfTimeSeed[100], ele_outOfTimeClusters[100];
  double ele_e2overe9[100] ;
  double ele_eidVeryLoose[100], ele_eidLoose[100], ele_eidMedium[100], ele_eidTight[100] ;
  double ele_eidHZZVeryLoose[100], ele_eidHZZLoose[100], ele_eidHZZMedium[100], ele_eidHZZTight[100] ;
  

  // RECO muons
  bool RECOMU_isGlobalMu[100],RECOMU_isStandAloneMu[100],RECOMU_isTrackerMu[100],RECOMU_isCaloMu[100];
  float RECOMU_E[100],RECOMU_PT[100],RECOMU_P[100],RECOMU_ETA[100],RECOMU_THETA[100],RECOMU_PHI[100],RECOMU_MASS[100],RECOMU_TMASS[100],RECOMU_QUALITY[100],RECOMU_CHARGE[100];
  float 
    RECOMU_TRACKISO[100],RECOMU_ECALISO[100],RECOMU_HCALISO[100], RECOMU_X[100],
    RECOMU_IP[100],RECOMU_SIP[100],RECOMU_IPERROR[100],
    RECOMU_IP_KF[100],RECOMU_SIP_KF[100],RECOMU_IPERROR_KF[100],
    RECOMU_STIP[100],RECOMU_TIP[100],RECOMU_TIPERROR[100],
    RECOMU_SLIP[100],RECOMU_LIP[100],RECOMU_LIPERROR[100],
    RECOMU_SIP_GD[100], RECOMU_SIP_GDMMMM[100],
    RECOMU_SIP_Std[100], RECOMU_SIP_StdMMMM[100], 
    RECOMU_SIP_Kin[100], RECOMU_SIP_KinMMMM[100],;
  
  // Particle Flow Isolation
  float RECOMU_ALLDEPOSITS[100],RECOMU_ChargedDEPOSITS[100],RECOMU_NeutralDEPOSITS[100],RECOMU_GammaDEPOSITS[100],RECOMU_PUDEPOSITS[100],RECOMU_PULowDEPOSITS[100],RECOMU_PFX_DB[100];
  
  float RECOELE_ALLDEPOSITS[100],RECOELE_ChargedDEPOSITS[100],RECOELE_NeutralDEPOSITS[100],RECOELE_GammaDEPOSITS[100],RECOELE_PUDEPOSITS[100],RECOELE_PULowDEPOSITS[100],RECOELE_PFX_DB[100];
  

  float
    RECOMU_caloCompatibility[100],RECOMU_segmentCompatibility[100];
  bool RECOMU_glbmuPromptTight[100];
  int RECOMUBEST_2e2mu_MATCHED[100],RECOMUBEST_4mu_MATCHED[100]; 
  int RECOMU_MMMM_MATCHED[100],RECOMU_EEMM_MATCHED[100],
      RECOMU_ZMM_MATCHED[100],RECOMU_ZssMM_MATCHED[100],RECOMU_ZEM_MATCHED[100],
      RECOMU_LLL0_MATCHED[100],RECOMU_LLL1_MATCHED[100],RECOMU_LLL2_MATCHED[100],RECOMU_LLL3_MATCHED[100],
    RECOMU_LLLLss0_MATCHED[100],RECOMU_LLLLss1_MATCHED[100],RECOMU_LLLLss2_MATCHED[100],
    RECOMU_LLLl0_MATCHED[100],RECOMU_LLLl1_MATCHED[100],RECOMU_LLLL_MATCHED[100];
  int RECOMU_numberOfMatches[100];
  
  float 
    RECOMU_mutrkDxy[100],RECOMU_mutrkDxyError[100],RECOMU_mutrkDxyB[100],
    RECOMU_mutrkDz[100],RECOMU_mutrkDzError[100],RECOMU_mutrkDzB[100],
    RECOMU_mutrkChi2PerNdof[100],RECOMU_mutrkCharge[100],RECOMU_mutrkNHits[100],RECOMU_mutrkNPixHits[100],RECOMU_mutrkNStripHits[100],RECOMU_mutrkNMuonHits[100];
  bool RECOMU_trkmuArbitration[100],RECOMU_trkmu2DCompatibilityLoose[100],RECOMU_trkmu2DCompatibilityTight[100];
  bool RECOMU_trkmuOneStationLoose[100],RECOMU_trkmuOneStationTight[100];
  bool RECOMU_trkmuLastStationLoose[100],RECOMU_trkmuLastStationTight[100];
  bool RECOMU_trkmuLastStationAngLoose[100],RECOMU_trkmuLastStationAngTight[100];
  bool RECOMU_trkmuOneStationAngLoose[100],RECOMU_trkmuOneStationAngTight[100];
  bool RECOMU_trkmuLastStationOptimizedLowPtLoose[100],RECOMU_trkmuLastStationOptimizedLowPtTight[100];
  
  // Isolation 
  float isolation[4],isoele[2],isomu[2];
  
  // Vertexing
  float vertexing[4];
  double ftsigma[4],ftsigmalag[4],ftsigmaMMMM[4],ftsigmalagMMMM[4],ftsigmaEEEE[4],ftsigmalagEEEE[4];
  double gdX[4],gdY[4],gdZ[4],gdXMMMM[4],gdYMMMM[4],gdZMMMM[4],gdXEEEE[4],gdYEEEE[4],gdZEEEE[4];
  double gdlagX[4],gdlagY[4],gdlagZ[4],gdlagProb[4],gdlagNdof[4],
    gdlagXMMMM[4],gdlagYMMMM[4],gdlagZMMMM[4],gdlagProbMMMM[4],gdlagNdofMMMM[4],
    gdlagXEEEE[4],gdlagYEEEE[4],gdlagZEEEE[4],gdlagProbEEEE[4],gdlagNdofEEEE[4];

  //ConstraintFit 4l
  double  StdFitVertexX[4], StdFitVertexY[4], StdFitVertexZ[4], StdFitVertexChi2r[4], StdFitVertexProb[4];
  double  KinFitVertexX[4], KinFitVertexY[4], KinFitVertexZ[4], KinFitVertexChi2r[4], KinFitVertexProb[4];
  double  StdFitVertexXMMMM[4], StdFitVertexYMMMM[4], StdFitVertexZMMMM[4], StdFitVertexChi2rMMMM[4], StdFitVertexProbMMMM[4];
  double  KinFitVertexXMMMM[4], KinFitVertexYMMMM[4], KinFitVertexZMMMM[4], KinFitVertexChi2rMMMM[4], KinFitVertexProbMMMM[4];
  double  StdFitVertexXEEEE[4], StdFitVertexYEEEE[4], StdFitVertexZEEEE[4], StdFitVertexChi2rEEEE[4], StdFitVertexProbEEEE[4];
  double  KinFitVertexXEEEE[4], KinFitVertexYEEEE[4], KinFitVertexZEEEE[4], KinFitVertexChi2rEEEE[4], KinFitVertexProbEEEE[4];

   //ConstraintFit 2l
  double  StdFitVertexChi2rDiLep[40], StdFitVertexProbDiLep[40];

  //ConstraintFit 3l
  double  StdFitVertexChi2rMMM[20], StdFitVertexProbMMM[20];
  double  StdFitVertexChi2rMME[20], StdFitVertexProbMME[20];
  double  StdFitVertexChi2rEEE[20], StdFitVertexProbEEE[20];
  double  StdFitVertexChi2rMEE[20], StdFitVertexProbMEE[20];
  

  // RECO counters
  int RECO_NMU, RECO_NELE, RECO_NTRACK, RECO_NPHOT, RECO_NJET, RECO_NVTX;
  float RECO_TRACK_PT[200], RECO_TRACK_ETA[200], RECO_TRACK_PHI[200],
    RECO_TRACK_CHI2[200],RECO_TRACK_CHI2RED[200],RECO_TRACK_CHI2PROB[200], 
    RECO_TRACK_DXY[200],RECO_TRACK_DXYERR[200], 
    RECO_TRACK_DZ[200],RECO_TRACK_DZERR[200];
  int RECO_TRACK_NHITS[200];
  
  // Primary Vertices
  float RECO_VERTEX_x[10], RECO_VERTEX_y[10], RECO_VERTEX_z[10],RECO_VERTEX_ndof[10],RECO_VERTEX_chi2[10],RECO_VERTEXPROB[10],RECO_VERTEX_TRACK_PT[10][30];
  bool RECO_VERTEX_isValid[10];
  int RECO_VERTEX_ntracks[10];
  
  // RECO JETS
  int RECO_PFJET_N, RECO_PFJET_CHARGE[100];
  float RECO_PFJET_ET[100], RECO_PFJET_PT[100], RECO_PFJET_ETA[100], RECO_PFJET_PHI[100];
  double RHO;

  // RECO MET
  float genmet;
  float calomet,calometopt,calometoptnohf,calometoptnohfho,calometoptho,calometnohf,calometnohfho,calometho;       
  float pfmet,pfmet_x,pfmet_y,pfmet_phi,pfmet_theta,htmetic5,htmetkt4,htmetkt6,htmetsc5,htmetsc7;        
  float tcmet,jescormetic5,jescormetkt4,jescormetkt6,jescormetsc5,jescormetsc7;    
  float cormetmuons;  
  
  // Beam Spot
  double BeamSpot_X,BeamSpot_Y,BeamSpot_Z;	
  
  // bTagging
  edm::InputTag    
    tCHighEff_bTag_, tCHighPur_bTag_,
    jPHighEff_bTag_,jBP_bTag_,
    sSVHighEff_bTag_,sSVHighPur_bTag_,
    cSV_bTag_,cSVMVA_bTag_, 
    sEByIP3d_bTag_,sEByPt_bTag_,     
    sM_bTag_,sMByIP3d_bTag_,sMByPt_bTag_;

  float tCHighEff_BTagJet_PT[100],
    tCHighPur_BTagJet_PT[100],
    jPHighEff_BTagJet_PT[100],
    jBP_BTagJet_PT[100],
    sSVHighEff_BTagJet_PT[100],
    sSVHighPur_BTagJet_PT[100],
    cSV_BTagJet_PT[100],
    cSVMVA_BTagJet_PT[100],
    sEByIP3d_BTagJet_PT[100],
    sEByPt_BTagJet_PT[100],
    sM_BTagJet_PT[100],
    sMByIP3d_BTagJet_PT[100],
    sMByPt_BTagJet_PT[100];

  float tCHighEff_BTagJet_ETA[100],
    tCHighPur_BTagJet_ETA[100],
    jPHighEff_BTagJet_ETA[100],
    jBP_BTagJet_ETA[100],
    sSVHighEff_BTagJet_ETA[100],
    sSVHighPur_BTagJet_ETA[100],
    cSV_BTagJet_ETA[100],
    cSVMVA_BTagJet_ETA[100],
    sEByIP3d_BTagJet_ETA[100],
    sEByPt_BTagJet_ETA[100],
    sM_BTagJet_ETA[100],
    sMByIP3d_BTagJet_ETA[100],
    sMByPt_BTagJet_ETA[100];

  float tCHighEff_BTagJet_PHI[100],
    tCHighPur_BTagJet_PHI[100],
    jPHighEff_BTagJet_PHI[100],
    jBP_BTagJet_PHI[100],
    sSVHighEff_BTagJet_PHI[100],
    sSVHighPur_BTagJet_PHI[100],
    cSV_BTagJet_PHI[100],
    cSVMVA_BTagJet_PHI[100],
    sEByIP3d_BTagJet_PHI[100],
    sEByPt_BTagJet_PHI[100],
    sM_BTagJet_PHI[100],
    sMByIP3d_BTagJet_PHI[100],
    sMByPt_BTagJet_PHI[100];

  float tCHighEff_BTagJet_DISCR[100],
    tCHighPur_BTagJet_DISCR[100],
    jPHighEff_BTagJet_DISCR[100],
    jBP_BTagJet_DISCR[100],
    sSVHighEff_BTagJet_DISCR[100],
    sSVHighPur_BTagJet_DISCR[100],
    cSV_BTagJet_DISCR[100],
    cSVMVA_BTagJet_DISCR[100],
    sEByIP3d_BTagJet_DISCR[100],
    sEByPt_BTagJet_DISCR[100],
    sM_BTagJet_DISCR[100],
    sMByIP3d_BTagJet_DISCR[100],
    sMByPt_BTagJet_DISCR[100];


  float ConvMapDist[100],ConvMapDcot[100];

  double compMt(const reco::Candidate::LorentzVector& visParticle, double metPx, double metPy) { 
    double px = visParticle.px() + metPx;
    double py = visParticle.py() + metPy;
    double et = visParticle.Et() + TMath::Sqrt(metPx*metPx + metPy*metPy);
    
    double mt2 = et*et - (px*px + py*py);
    cout<<"transverse mass = "<<mt2<<endl;
    if ( mt2 < 0 ) { 
      cout<<"compMt " << " mt2 = " << mt2 << " must not be negative !!"<<endl;
      return 0.;
    }
    return TMath::Sqrt(mt2);
  }

};

#endif


