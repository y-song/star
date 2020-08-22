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

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
/* #include "Pythia8Plugins/FastJet3.h" */
#include "Math/ProbFunc.h"

using namespace Pythia8;
using namespace fastjet;
using namespace std;


int main() {
    
    TFile* f_out = new TFile ("$HOME/output/main.root","recreate");

    // New TTree
    UInt_t ntrials, evid, ncharged, nneutral, nnh, nch, n_decay_photon;
    Float_t xsec, x, y, Q2, W2, e_jet, pt_jet, eta_jet, phi_jet, p_jet, theta_jet;
    TTree *tree = new TTree("Tree", "Tree");
    tree->Branch("ntrials", &ntrials, "ntrials/I");
    tree->Branch("evid", &evid, "evid/I");
    tree->Branch("ncharged", &ncharged, "ncharged/I");
    tree->Branch("nneutral", &nneutral, "nneutral/I");
    tree->Branch("nnh", &nnh, "nnh/I");
    tree->Branch("nch", &nch, "nch/I");
    tree->Branch("n_decay_photon", &n_decay_photon, "n_decay_photon/I");
    tree->Branch("xsec", &xsec, "xsec/F");
    tree->Branch("x", &x, "x/F");
    tree->Branch("y", &y, "y/F");
    tree->Branch("Q2", &Q2, "Q2/F");
    tree->Branch("W2", &W2, "W2/F");
    tree->Branch("e_jet", &e_jet, "e_jet/F");
    tree->Branch("pt_jet", &pt_jet, "pt_jet/F");
    tree->Branch("eta_jet", &eta_jet, "eta_jet/F");
    tree->Branch("phi_jet", &phi_jet, "phi_jet/F");
    tree->Branch("p_jet", &p_jet, "p_jet/F");
    tree->Branch("theta_jet", &theta_jet, "theta_jet/F");

    // New TH1D
    const int    nbins_pt { 400 };
    const double lobin_pt { 0. };
    const double hibin_pt {100. * 260100/200 };
    TH1D* kshort { new TH1D("kshort", "K0 short multiplicity in |eta|<0.5;#it{p}_{T};dK/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* kaon   { new TH1D("kaon", "kaon multiplicity in |eta|<0.5;#it{p}_{T};dK^{-}/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* antikaon   { new TH1D("antikaon", "antikaon multiplicity in |eta|<0.5;#it{p}_{T};dK/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* pion   { new TH1D("pion", "pion multiplicity in |eta|<0.5;#it{p}_{T};d#pi/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* antipion   { new TH1D("antipion", "anti-pion multiplicity in |eta|<0.5;#it{p}_{T};d#pi^{-}/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* proton   { new TH1D("proton",  "proton multiplicity in |eta|<0.5;#it{p}_{T};d#p/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) };
    TH1D* pbar   { new TH1D("pbar", "#bar{p} multiplicity in |eta|<0.5;#it{p}_{T};d#bar{p}/d#it{p}_{T}",nbins_pt,lobin_pt,hibin_pt) }; 

    // Select FastJet parameters
    const double R_FULL       = 0.4;    // Jet size.
    const double R_CH         = 0.4;    // Jet size.
    const double pTMin        = 0.2;    // Min jet pT.
    const double etaMax { 10. };
    fastjet::JetDefinition jetdef (fastjet::antikt_algorithm, R_FULL);
    fastjet::Selector jetrap  = fastjet::SelectorAbsRapMax(etaMax)  ; // for full jets
   
    // Initialize Pythia 
    Pythia pythia;
    const double sNN { 260100 };
    const double pTHatMin { 1000. };
    const double pTHatMax { 2000. };
    string name_type = "AuAu";
    int idA; int idB;
    TRandom3 rand{0};   
 
    if      (name_type == "pp" )  { idA = 2212; idB = 2212; }
    else if (name_type == "pAu")  { idA = 2212; idB = 1000822080; }
    else if (name_type == "AuAu") { idA = 1000822080; idB = 1000822080; }

    for (auto str : vector<string>{ 
            Form("Beams:eCM = %f",sNN),
            "HardQCD:all = on",
            //Form("PhaseSpace:pTHatMin = %f",pTHatMin),
            //Form("PhaseSpace:pTHatMax = %f",pTHatMax),
            Form("Beams:idA = %i", idA), // moving in +z direction
            Form("Beams:idB = %i", idB), // moving in +z direction
            "Random:setSeed = on",
            "ParticleDecays:limitRadius = on",
            Form("ParticleDecays:rMax = %f", 10.),
            Form("Random:seed = %i",rand.Integer(10000000)),
    }) pythia.readString(str.c_str());

    pythia.init();
    Pythia8::Info& info = pythia.info ;

    cout << " Starting the pythia loop. " << endl;

    int n_events{1000};
    const double mips_min_p { 0.2 };
    for (int iEvent{0}; iEvent < n_events; ++iEvent) {
        if (!(iEvent % 10000)) {
            cout << "Finished event " << iEvent << endl;
        }
        if (!pythia.next()) continue;

        Event& event = pythia.event;
        std::vector <fastjet::PseudoJet> part_FULL;
        std::vector <fastjet::PseudoJet> part_CH;
        float max_bemc_Et = 0.;

        for (int i {0}; i < event.size(); ++i) {
            if (!event[i].isFinal()) continue;
            auto& e { event[i] }; 
            if (!e.isVisible()) continue;

            double pAbs   { e.pAbs() };
            if (pAbs < mips_min_p) continue;

            double eta { e.eta() };
            double pt  { e.pT()  };
            if ( fabs(eta) > etaMax ) continue;
            if ( fabs(eta) < 0.5 ) {
                int id { e.id() };
                if      (id ==  310 ) kshort  ->Fill(pt);
                else if (id ==  321 ) kaon    ->Fill(pt);
                else if (id == -321 ) antikaon->Fill(pt);
                else if (id ==  2212) proton  ->Fill(pt);
                else if (id == -2212) pbar    ->Fill(pt);
                else if (id ==  211 ) pion    ->Fill(pt);
                else if (id == -211 ) antipion->Fill(pt);
            }

            PseudoJet p {  e.px(), e.py(), e.pz(), e.pAbs() }; // very lazy way to get psuedorapidity
            part_FULL.push_back(p);
            if (e.isCharged()) {
                part_CH.push_back(p);
            } else {
                double et { pt / TMath::CosH(eta) };
                if (fabs(eta)<1 && et > max_bemc_Et) { max_bemc_Et = et; };
            }
        }
        
        // make and use full jets
        fastjet::ClusterSequence clustSeq_FULL(part_FULL, jetdef);
        vector<fastjet::PseudoJet> jetsFULL = sorted_by_pt(jetrap(clustSeq_FULL.inclusive_jets(pTMin)));
        for (unsigned int ij = 0; ij < jetsFULL.size(); ij++) {
		evid = iEvent;
		xsec = pythia.info.sigmaGen();
		ntrials = pythia.info.nTried();
		e_jet = jetsFULL[ij].e();
		pt_jet = jetsFULL[ij].perp();
		eta_jet = jetsFULL[ij].eta();
		phi_jet = jetsFULL[ij].phi();
		p_jet = jetsFULL[ij].perp() * TMath::CosH(jetsFULL[ij].eta());
		theta_jet = jetsFULL[ij].theta();
		
		tree->Fill(); 
	}  

        // make and use charged jets
        fastjet::ClusterSequence clustSeq_CH(part_CH, jetdef);
        vector<fastjet::PseudoJet> jetsCH = sorted_by_pt(jetrap(clustSeq_CH.inclusive_jets(pTMin)));
        // do something with these jets...
    }
    cout << " finished PYTHIA loop " << endl;

    f_out->Write();
    f_out->Save();
    f_out->Close();

    return 0;
};