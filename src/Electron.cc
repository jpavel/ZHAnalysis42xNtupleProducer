#include "Analysis/NtupleProducer/interface/NtupleProducer.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

void NtupleProducer::DoElectronAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    (m->PreSelectedElectrons).clear();

    edm::Handle< pat::TriggerEvent > triggerEvent;
    iEvent.getByLabel(triggerEvent_, triggerEvent);
    // PAT trigger helper for trigger matching information
    const pat::helper::TriggerMatchHelper matchHelper;

     edm::Handle<edm::ValueMap<float> >  mvaTrigV0_;
    iEvent.getByLabel("mvaTrigV0", mvaTrigV0_);

    edm::Handle<edm::ValueMap<float> >  mvaNonTrigV0_;
    iEvent.getByLabel("mvaNonTrigV0", mvaNonTrigV0_);


    edm::Handle<edm::ValueMap<double> > isoAllEleMap;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueChargedAll04PFIdPFIso"), isoAllEleMap);

    edm::Handle<edm::ValueMap<double> > isoChargedEleMap;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueCharged04PFIdPFIso"), isoChargedEleMap);

    edm::Handle<edm::ValueMap<double> > isoNeutralEleMap;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueNeutral04PFIdPFIso"), isoNeutralEleMap);

    edm::Handle<edm::ValueMap<double> > isoGammaEleMap;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueGamma04PFIdPFIso"), isoGammaEleMap);

    edm::Handle<edm::ValueMap<double> > isoPUEleMap;
    iEvent.getByLabel(edm::InputTag("elPFIsoValuePU04PFIdPFIso"), isoPUEleMap);



    edm::Handle<edm::ValueMap<double> > isoAllEleMap_NoPFId;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueChargedAll04NoPFIdPFIso"), isoAllEleMap_NoPFId);

    edm::Handle<edm::ValueMap<double> > isoChargedEleMap_NoPFId;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueCharged04NoPFIdPFIso"), isoChargedEleMap_NoPFId);

    edm::Handle<edm::ValueMap<double> > isoNeutralEleMap_NoPFId;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueNeutral04NoPFIdPFIso"), isoNeutralEleMap_NoPFId);

    edm::Handle<edm::ValueMap<double> > isoGammaEleMap_NoPFId;
    iEvent.getByLabel(edm::InputTag("elPFIsoValueGamma04NoPFIdPFIso"), isoGammaEleMap_NoPFId);

    edm::Handle<edm::ValueMap<double> > isoPUEleMap_NoPFId;
    iEvent.getByLabel(edm::InputTag("elPFIsoValuePU04NoPFIdPFIso"), isoPUEleMap_NoPFId);

    // for conversion veto selection
   edm::Handle<reco::ConversionCollection> hConversions;
   iEvent.getByLabel("allConversions", hConversions);



//    edm::Handle<edm::ValueMap<double> > isoPULowEleMap;
//    iEvent.getByLabel(edm::InputTag("electronPFIsoValuePULow"), isoPULowEleMap);

//        edm::Handle<edm::ValueMap<double> > isoAllEleMap;
//    iEvent.getByLabel(edm::InputTag("electronPFIsoValueAll"), isoAllEleMap);
//
//    edm::Handle<edm::ValueMap<double> > isoChargedEleMap;
//    iEvent.getByLabel(edm::InputTag("electronPFIsoValueCharged"), isoChargedEleMap);
//
//    edm::Handle<edm::ValueMap<double> > isoNeutralEleMap;
//    iEvent.getByLabel(edm::InputTag("electronPFIsoValueNeutral"), isoNeutralEleMap);
//
//    edm::Handle<edm::ValueMap<double> > isoGammaEleMap;
//    iEvent.getByLabel(edm::InputTag("electronPFIsoValueGamma"), isoGammaEleMap);
//
//    edm::Handle<edm::ValueMap<double> > isoPUEleMap;
//    iEvent.getByLabel(edm::InputTag("electronPFIsoValuePU"), isoPUEleMap);
//
//    edm::Handle<edm::ValueMap<double> > isoPULowEleMap;
//    iEvent.getByLabel(edm::InputTag("electronPFIsoValuePULow"), isoPULowEleMap);

    Handle<pat::ElectronCollection> ElectronsHandle;
    iEvent.getByLabel(PreSelectedElectronCollection_, ElectronsHandle);
    const pat::ElectronCollection gsfElectrons = *(ElectronsHandle.product());

    pat::ElectronCollection::const_iterator iElectron;


    int index = 0;


    // Get the B-field
    ESHandle<MagneticField> B;
    iSetup.get<IdealMagneticFieldRecord > ().get(B);
    const MagneticField* magField = B.product();

    //BeamSpot
    Handle<reco::BeamSpot> theBeamSpotHandle;
    iEvent.getByLabel(std::string("offlineBeamSpot"), theBeamSpotHandle);
    const reco::BeamSpot* beamSpot = theBeamSpotHandle.product();




    // Get the geometry
    ESHandle<GlobalTrackingGeometry> geomHandle;
    iSetup.get<GlobalTrackingGeometryRecord > ().get(geomHandle);

    // get the track builder
    ESHandle<TransientTrackBuilder> trackBuilder;
    iSetup.get<TransientTrackRecord > ().get("TransientTrackBuilder", trackBuilder);

    // Beamspot
    Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByType(recoBeamSpotHandle);
    const BeamSpot bs = *recoBeamSpotHandle;

    GlobalPoint BSVertex(bs.position().x(), bs.position().y(), bs.position().z());
    Basic3DVector<double> BSVertexErr(bs.x0Error(), bs.y0Error(), bs.z0Error());

    reco::Vertex::Point BSp(bs.position().x(), bs.position().y(), bs.position().z());
    reco::Vertex::Error BSe;

    BSe(0, 0) = bs.x0Error() * bs.x0Error();
    BSe(1, 1) = bs.y0Error() * bs.y0Error();
    BSe(2, 2) = bs.z0Error() * bs.z0Error();
    reco::Vertex BSprimVertex = reco::Vertex(BSp, BSe, 1, 1, 1);


    // Get primary vertex collection
    Handle<reco::VertexCollection> recoPVCollection;
    iEvent.getByLabel(vertexCollectionForLeptonIP_, recoPVCollection);

    //  useBeamSpot_  = pset.getParameter<bool>("useBeamSpot");
    bool useBeamSpot_ = true;

    //
    reco::Vertex primVertex;
    bool pvfound = (recoPVCollection->size() != 0);

    if (pvfound) {
        PrimaryVertexSorter pvs;
        vector<reco::Vertex> sortedList = pvs.sortedList(*(recoPVCollection.product()));
        primVertex = (sortedList.front());
    } else {
        //creating a dummy PV
        reco::Vertex::Point p(0, 0, 0);
        reco::Vertex::Error e;
        e(0, 0) = 0.0015 * 0.0015;
        e(1, 1) = 0.0015 * 0.0015;
        e(2, 2) = 15. * 15.;
        primVertex = reco::Vertex(p, e, 1, 1, 1);
    }
    //
    GlobalPoint pVertex(primVertex.position().x(), primVertex.position().y(), primVertex.position().z());
    Basic3DVector<double> pVertexErr(primVertex.xError(), primVertex.yError(), primVertex.zError());


    //--track refs
    GsfTrackRef eletrack;
    //--transient tracks
    reco::TransientTrack eletranstrack;
    //--signif. vars
    float eleSignificance3D;
    float eleValue3D, eleError3D;



    for (iElectron = gsfElectrons.begin(); iElectron != gsfElectrons.end(); iElectron++, index++) {

        myobject elo;
        elo.pt = iElectron->pt();
        elo.eta = iElectron->eta();
        elo.phi = iElectron->phi();
        elo.px = iElectron->px();
        elo.py = iElectron->py();
        elo.pz = iElectron->pz();
        elo.z = iElectron->vz();


        elo.E = iElectron->p();
        elo.et = iElectron->et();
        elo.mass = iElectron->mass();
        elo.mt = iElectron->mt();
        elo.Energy = iElectron->energy();
        elo.charge = iElectron->charge();

        elo.numHitEleInner = (iElectron->gsfTrack().isNonnull() ? iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits() : 0.);
        elo.numLostHitEleInner = (iElectron->gsfTrack().isNonnull() ? iElectron->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() : 0.);

        elo.numLostHitEle = (iElectron->gsfTrack().isNonnull() ? iElectron->gsfTrack()->numberOfLostHits() : 0.);
        elo.numValidHitEle = (iElectron->gsfTrack().isNonnull() ? iElectron->gsfTrack()->numberOfValidHits() : 0.);
      	elo.IP3D = iElectron->dB(pat::Electron::PV3D);
	elo.dxy_PV = (iElectron->gsfTrack().isNonnull() ? iElectron->gsfTrack()->dxy(primVertex.position()) : 0.);
	elo.dz_PV = (iElectron->gsfTrack().isNonnull() ? iElectron->gsfTrack()->dz(primVertex.position()) : 0.);

        elo.HoverE = iElectron->hadronicOverEm();
        elo.deltaPhiSuperClusterTrackAtVtx = iElectron->deltaPhiSuperClusterTrackAtVtx();
        elo.deltaEtaSuperClusterTrackAtVtx = iElectron->deltaEtaSuperClusterTrackAtVtx();
        elo.sigmaIetaIeta = iElectron->sigmaIetaIeta();
        elo.sigmaEtaEta = iElectron->sigmaEtaEta();
        elo.ecalIso = iElectron->ecalIso();
        elo.hcalIso = iElectron->hcalIso();
        elo.trackIso = iElectron->trackIso();
        elo.caloIso = iElectron->caloIso();
        elo.ecalEnergy = iElectron->ecalEnergy();
        elo.hcalOverEcal = iElectron->hcalOverEcal();
        elo.eta_SC       = iElectron->superCluster()->eta()   ;
	elo.rawE_SC  = iElectron->superCluster()->rawEnergy(); 
	elo.preshowerE_SC = iElectron->superCluster()->preshowerEnergy();

 



//        elo.EleId95rel = iElectron->electronID("simpleEleId95relIso");
//        elo.EleId90rel = iElectron->electronID("simpleEleId90relIso");
//        elo.EleId85rel = iElectron->electronID("simpleEleId85relIso");
//        elo.EleId80rel = iElectron->electronID("simpleEleId80relIso");
//        elo.EleId70rel = iElectron->electronID("simpleEleId70relIso");
//        elo.EleId60rel = iElectron->electronID("simpleEleId60relIso");

//        elo.EleId95cIso = iElectron->electronID("simpleEleId95cIso");
//        elo.EleId90cIso = iElectron->electronID("simpleEleId90cIso");
//        elo.EleId85cIso = iElectron->electronID("simpleEleId85cIso");
//        elo.EleId80cIso = iElectron->electronID("simpleEleId80cIso");
//        elo.EleId70cIso = iElectron->electronID("simpleEleId70cIso");
//        elo.EleId60cIso = iElectron->electronID("simpleEleId60cIso");

//         elo.CicVeryLoose = iElectron->electronID("eidVeryLoose");
//         elo.CicLoose = iElectron->electronID("eidLoose");
//         elo.CicMedium = iElectron->electronID("eidMedium");
//         elo.CicTight = iElectron->electronID("eidTight");
//         elo.CicSuperTight = iElectron->electronID("eidSuperTight");

//        elo.CicHZZVeryLoose = iElectron->electronID("eidHZZVeryLoose");
//        elo.CicHZZLoose = iElectron->electronID("eidHZZLoose");
//        elo.CicHZZMedium = iElectron->electronID("eidHZZMedium");
//        elo.CicHZZTight = iElectron->electronID("eidHZZTight");
//        elo.CicHZZSuperTight = iElectron->electronID("eidHZZSuperTight");








         Handle<reco::GsfElectronCollection> ElectronsHandleGSF;
         iEvent.getByLabel("gsfElectrons", ElectronsHandleGSF);
         const reco::GsfElectronCollection gsfElectrons_ = *(ElectronsHandleGSF.product());
        

          reco::GsfElectronCollection::const_iterator itr_elec;
       


        elo.Id_mvaTrg = -10;
        elo.Id_mvaNonTrg = -10;
        int num_debug=0;
        int num=-5;
          for (itr_elec = gsfElectrons_.begin(); itr_elec != gsfElectrons_.end(); itr_elec++, num_debug++) {

              if (reco::deltaR( itr_elec->p4() , iElectron->p4()) < 0.001){
		elo.isGsfCtfScPixChargeConsistent= itr_elec->isGsfCtfScPixChargeConsistent() ;
		elo.isGsfScPixChargeConsistent = itr_elec->isGsfScPixChargeConsistent();
		elo.isGsfCtfChargeConsistent = itr_elec->isGsfCtfChargeConsistent();

   edm::Ref<reco::GsfElectronCollection> eleRefgsf(ElectronsHandleGSF, num_debug);
        elo.Id_mvaTrg = (*mvaTrigV0_)[eleRefgsf];
        elo.Id_mvaNonTrg = (*mvaNonTrigV0_)[eleRefgsf];
        num = num_debug;
        break;
}
             }
//cout<<"pt= "<<elo.pt << "      MVAnon=  "<<elo.Id_mvaNonTrg <<"     index= " << index  <<"      num_debug" << num <<endl;

                edm::Ref<pat::ElectronCollection> eleRef(ElectronsHandle, index);
        elo.pfIsoAll = (*isoAllEleMap)[eleRef];
        elo.pfIsoCharged = (*isoChargedEleMap)[eleRef];
        elo.pfIsoNeutral = (*isoNeutralEleMap)[eleRef];
        elo.pfIsoGamma = (*isoGammaEleMap)[eleRef];
        elo.pfIsoPU = (*isoPUEleMap)[eleRef];
//        elo.pfIsoPULow = (*isoPULowEleMap)[eleRef];


                elo.pfIsoAll_NoPFId = (*isoAllEleMap_NoPFId)[eleRef];
        elo.pfIsoCharged_NoPFId = (*isoChargedEleMap_NoPFId)[eleRef];
        elo.pfIsoNeutral_NoPFId = (*isoNeutralEleMap_NoPFId)[eleRef];
        elo.pfIsoGamma_NoPFId = (*isoGammaEleMap_NoPFId)[eleRef];
        elo.pfIsoPU_NoPFId = (*isoPUEleMap_NoPFId)[eleRef];


        eletrack = iElectron->get<GsfTrackRef > ();
        eletranstrack = trackBuilder->build(eletrack);

        TrajectoryStateOnSurface eleTSOS;

        if (useBeamSpot_ == true) {
            eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), BSVertex, eletranstrack.field());
        } else {
            eleTSOS = IPTools::transverseExtrapolate(eletranstrack.impactPointState(), pVertex, eletranstrack.field());
        }


        if (eleTSOS.isValid()) {
            std::pair<bool, Measurement1D> eleIPpair;
            if (useBeamSpot_ == true) {
                eleIPpair = IPTools::signedImpactParameter3D(eletranstrack, eleTSOS.globalDirection(), BSprimVertex);
            } else {
                eleIPpair = IPTools::signedImpactParameter3D(eletranstrack, eleTSOS.globalDirection(), primVertex);
            }
            //
            if (eleIPpair.first) {
                eleSignificance3D = eleIPpair.second.significance();
                eleValue3D = eleIPpair.second.value();
                eleError3D = eleIPpair.second.error();
                elo.SIP = eleSignificance3D;
            }
        }



        //DZ from Ian
        elo.z_expo = 0;
        if (iElectron->gsfTrack().isNonnull()) {
            reco::GsfTransientTrack track(iElectron->gsfTrack(), magField, geomHandle);
            //            FreeTrajectoryState state = track.impactPointTSCP().theState();
            //          states.push_back(state);
            //calculate Z at beamspot
            TransverseImpactPointExtrapolator extrapolator(magField);
            TrajectoryStateOnSurface closestOnTransversePlaneState = extrapolator.extrapolate(track.impactPointState(), GlobalPoint(beamSpot->position().x(), beamSpot->position().y(), 0.0));
           elo.z_expo = closestOnTransversePlaneState.globalPosition().z();
	   // conversion veto
	    bool passconversionveto = false;
	    if( hConversions.isValid()){
	        // this is recommended method
	        passconversionveto = !ConversionTools::hasMatchedConversion( dynamic_cast<reco::GsfElectron const&>(*(iElectron->originalObjectRef())), hConversions, beamSpot->position());
	      }
	    elo.passConversionVeto = passconversionveto;


        }

	
	const pat::TriggerObjectRef trigRef_loose(matchHelper.triggerMatchObject(ElectronsHandle, index, electronMatch_Loose_, iEvent, *triggerEvent));
        elo.hasTrgObject_loose = false;
        elo.TrgObjectEta_loose = -100;
        elo.TrgObjectPt_loose = -100;
        elo.TrgObjectPhi_loose = -100;
        // finally we can fill the histograms
        if (trigRef_loose.isAvailable()) { // check references (necessary!)

	  elo.hasTrgObject_loose = true;
	  elo.TrgObjectEta_loose = trigRef_loose->eta();
	  elo.TrgObjectPt_loose = trigRef_loose->pt();
	  elo.TrgObjectPhi_loose = trigRef_loose->phi();
        }


	const pat::TriggerObjectRef trigRef_medium(matchHelper.triggerMatchObject(ElectronsHandle, index, electronMatch_Medium_, iEvent, *triggerEvent));
	        elo.hasTrgObject_medium = false;
	elo.TrgObjectEta_medium = -100;
	elo.TrgObjectPt_medium = -100;
	        elo.TrgObjectPhi_medium = -100;
	// finally we can fill the histograms
	if (trigRef_medium.isAvailable()) { // check references (necessary!)
	            elo.hasTrgObject_medium = true;
	    elo.TrgObjectEta_medium = trigRef_medium->eta();
	            elo.TrgObjectPt_medium = trigRef_medium->pt();
	    elo.TrgObjectPhi_medium = trigRef_medium->phi();
	  }



        (m->PreSelectedElectrons).push_back(elo);

    }//for loop on electrons



}
