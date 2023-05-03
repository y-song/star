#include "TMath.h"
#include "Pythia8/Pythia.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <array>

#include "TRandom3.h"
#include <sstream>
#include <iostream>
#include <cstring>

#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "Pythia8Plugins/FastJet3.h"
#include "Math/ProbFunc.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;

int main()
{

    TFile *f_out = new TFile("$HOME/output/pthat10_test.root", "recreate");

    // New TTree
    UInt_t ntrials, evid, njet, ncharged[100], nneutral[100];
    Float_t xsec, x1, x2, pt3, pt4, eta3, eta4, phi3, phi4, m3, m4, et;
    Float_t e_jet[100], pt_jet[100], eta_jet[100], phi_jet[100], p_jet[100], m_jet[100], q_jet[100];
 
    TTree *tree = new TTree("Tree", "Tree");

    tree->Branch("ntrials", &ntrials, "ntrials/I");
    tree->Branch("evid", &evid, "evid/I");
    tree->Branch("xsec", &xsec, "xsec/F");
    tree->Branch("x1", &x1, "x1/F");
    tree->Branch("x2", &x2, "x2/F");
    tree->Branch("pt3", &pt3, "pt3/F");
    tree->Branch("pt4", &pt4, "pt4/F");
    tree->Branch("eta3", &eta3, "eta3/F");
    tree->Branch("eta4", &eta4, "eta4/F");
    tree->Branch("phi3", &phi3, "phi3/F");
    tree->Branch("phi4", &phi4, "phi4/F");
    tree->Branch("m3", &m3, "m3/F");
    tree->Branch("m4", &m4, "m4/F");
    tree->Branch("et", &et, "et/F");
    
    tree->Branch("njet", &njet, "njet/I");
    tree->Branch("e_jet", e_jet, "e_jet[njet]/F");
    tree->Branch("pt_jet", pt_jet, "pt_jet[njet]/F");
    tree->Branch("eta_jet", eta_jet, "eta_jet[njet]/F");
    tree->Branch("phi_jet", phi_jet, "phi_jet[njet]/F");
    tree->Branch("p_jet", p_jet, "p_jet[njet]/F");
    tree->Branch("m_jet", m_jet, "m_jet[njet]/F");
    tree->Branch("q_jet", q_jet, "q_jet[njet]/F");
    tree->Branch("ncharged", ncharged, "ncharged[njet]/I");
    tree->Branch("nneutral", nneutral, "nneutral[njet]/I");

   // Select FastJet parameters
    const double R_FULL = 0.5; // Jet size.
    const double R_CH = 0.5;   // Jet size.
    const double pTMin = 0.2;  // Min jet pT.
    const double etaMax{4.0};
    const double etaMin{-4.0};
    fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, R_FULL);
    fastjet::Selector jetrap = fastjet::SelectorEtaMin(etaMin + R_FULL) * fastjet::SelectorEtaMax(etaMax - R_FULL); // for full jets
    //fastjet::Selector forward = fastjet::SelectorEtaMin(2.5+R_FULL) * fastjet::SelectorEtaMax(4.0-R_FULL);

    // Initialize Pythia
    Pythia pythia;
    const double sNN{500};
    const double pTHatMin{10.};
    const double pTHatMax{20.};
    string name_type = "pp";
    int idA;
    int idB;
    TRandom3 rand{0};
    TRandom *gaus = new TRandom();
    TRandom *integer = new TRandom();

    if (name_type == "pp")
    {
        idA = 2212;
        idB = 2212;
    }
    else if (name_type == "pAu")
    {
        idA = 2212;
        idB = 1000822080;
    }
    else if (name_type == "AuAu")
    {
        idA = 1000822080;
        idB = 1000822080;
    }

    for (auto str : vector<string>{
             Form("Beams:eCM = %f", sNN),
             "SoftQCD:inelastic = on",
             Form("PhaseSpace:pTHatMin = %f",pTHatMin),
             //Form("PhaseSpace:pTHatMax = %f",pTHatMax),
             Form("Beams:idA = %i", idA), // moving in +z direction
             Form("Beams:idB = %i", idB), // moving in +z direction
             "Random:setSeed = on",
             "ParticleDecays:limitRadius = on",
             Form("ParticleDecays:rMax = %f", 10.),
             Form("Random:seed = %i", rand.Integer(10000000)),
             "Next:numberCount = 10000",
         })
        pythia.readString(str.c_str());

    pythia.init();
    Pythia8::Info &info = pythia.info;

    cout << " Starting the pythia loop. " << endl;

    int n_events{20000};//0000000};
    const double mips_min_p{0.2};

    for (int iEvent{0}; iEvent < n_events; ++iEvent)
    {
        if (!pythia.next())
            continue;

        Event &event = pythia.event;
        std::vector<fastjet::PseudoJet> part_FULL;
       
        evid = iEvent;
        xsec = pythia.info.sigmaGen();
        ntrials = pythia.info.nTried();
 	et = 0.;
	
        // partons 
        x1 = pythia.info.x1();
	x2 = pythia.info.x2();
	pt3 = event[5].pT();
	pt4 = event[6].pT();
	eta3 = event[5].eta();
	eta4 = event[6].eta();
	phi3 = event[5].phi();
	phi4 = event[6].phi();
	m3 = event[5].m();
	m4 = event[6].m();
			
       // particle loop
        for (int i{0}; i < event.size(); ++i)
        {
            if (!event[i].isFinal())
                continue;
            auto &e{event[i]};
            if (!e.isVisible())
                continue;

            double pAbs{e.pAbs()};
            if (pAbs < mips_min_p)
                continue;

            double eta{e.eta()};
            double pt{e.pT()};
            if (fabs(eta) > etaMax || fabs(eta) < etaMin)
                continue;
	    
	    if (eta > 2.5 && eta < 4.0){	    
	    	et += e.eT();
	    }

	    PseudoJet p = event[i];
	    //PseudoJet p{e.px(), e.py(), e.pz(), e.e()};
            part_FULL.push_back(p);
  	}

        // make and use full jets
        fastjet::ClusterSequence clustSeq_FULL(part_FULL, jetdef);
	vector<fastjet::PseudoJet> jetsFULL = sorted_by_pt(jetrap(clustSeq_FULL.inclusive_jets(pTMin)));
	njet = jetsFULL.size();
	if (njet > 100)	continue;
	//jet loop
	for (unsigned int ij = 0; ij < jetsFULL.size(); ij++)
        {
            e_jet[ij] = jetsFULL[ij].e();
            pt_jet[ij] = jetsFULL[ij].perp();
            eta_jet[ij] = jetsFULL[ij].eta();
            phi_jet[ij] = jetsFULL[ij].phi();
            p_jet[ij] = jetsFULL[ij].perp() * TMath::CosH(jetsFULL[ij].eta());
	    m_jet[ij] = jetsFULL[ij].m();
   
	    //jet constituents
	    vector<fastjet::PseudoJet> constituents = jetsFULL[ij].constituents();
	    vector<fastjet::PseudoJet> charged_parts = fastjet::SelectorIsCharged()(constituents);
	    ncharged[ij] = charged_parts.size();
	    nneutral[ij] = constituents.size() - charged_parts.size();
	    q_jet[ij] = 0;

    	    if (charged_parts.size() == 0)	continue;
	    
	    double numerator = 0;
	    for (int i = 0; i < charged_parts.size(); i++){
		
		fastjet::PseudoJet charged_part = charged_parts.at(i);
		double charge = charged_part.user_info<Pythia8::Particle>().charge();
		numerator += charged_part.perp() * charge;
	     }

	     q_jet[ij] = numerator / pt_jet[ij];
	}
	tree -> Fill();
  
    }
    
    cout << " finished PYTHIA loop " << endl;

    f_out->Write();
    f_out->Save();
    f_out->Close();

    return 0;
};