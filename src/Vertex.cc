#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoVertexAnalysis(const edm::Event& iEvent) {
    (m->Vertex).clear();

    (m->runNumber) = 0;
    (m->eventNumber) = 0;
    (m->lumiNumber) = 0;
    m->runNumber = iEvent.id().run();
    m->eventNumber = iEvent.id().event();
    m->lumiNumber = iEvent.id().luminosityBlock();

    m->PUInfo = 0;
    m->PUInfo_Bunch0 = 0;

    Handle<VertexCollection> vertexHandle;
    iEvent.getByLabel(VertexCollection_, vertexHandle);
    const reco::VertexCollection vertexCollection = *(vertexHandle.product());
    VertexCollection::const_iterator iVertex;

    for (iVertex = vertexCollection.begin(); iVertex != vertexCollection.end(); iVertex++) {

        myobject vo;
        vo.Num_Vertex = vertexCollection.size();
        vo.px = iVertex->x();
        vo.py = iVertex->y();
        vo.pz = iVertex->z();
        vo.isFake = iVertex->isFake();
        vo.isValid = iVertex->isValid();
        vo.normalizedChi2 = iVertex->normalizedChi2();
        vo.ndof = iVertex->ndof();

        //new parameter 08/04/2012
//        vo.position_Rho = iVertex->position().Rho();
//        vo.position_rho = iVertex->position().rho();
//        vo.tracksSize = iVertex->tracksSize();

        (m->Vertex).push_back(vo);
    }


    if (IsMC) {
        //    edm::InputTag PileupSrc_ = "addPileupInfo";
        Handle<std::vector< PileupSummaryInfo > > PupInfo;
        iEvent.getByLabel(PileUpInfo_, PupInfo);

        std::vector<PileupSummaryInfo>::const_iterator PVI;
        
        for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
//            cout<< "NumInt= "<<PVI->getPU_NumInteractions()<<"\n";
            
            if(PVI->getBunchCrossing()== 0){  //new from Raman
	      m->PUInfo_Bunch0 = PVI->getPU_NumInteractions();
	      m->PUInfo_true = PVI->getTrueNumInteractions();
            }
m->PUInfo = PVI->getPU_NumInteractions();
            }
        }

//        Handle<std::vector< PileupSummaryInfo > >  PupInfo;
//iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
//
//std::vector<PileupSummaryInfo>::const_iterator PVI;
//
//float Tnpv = -1;
//for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
//
//   int BX = PVI->getBunchCrossing();
//
//   if(BX == 0) {
//     Tnpv = PVI->getTrueNumInteractions();
//     continue;
//   }
//
//}





    }


