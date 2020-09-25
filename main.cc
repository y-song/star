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

int main()
{

    TFile *f_out = new TFile("$HOME/output/main.root", "recreate");

    // New TTree
    UInt_t ntrials, evid, ncharged, nneutral, nnh, nch, n_decay_photon, njet, njet_s;
    Float_t xsec, x, y, Q2, W2;
    Float_t e_jet[30], pt_jet[30], eta_jet[30], phi_jet[30], p_jet[30], theta_jet[30], e_jet_s[30], pt_jet_s[30], eta_jet_s[30], phi_jet_s[30], p_jet_s[30],theta_jet_s[30];
 
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
    tree->Branch("njet", &njet, "njet/I");
    tree->Branch("njet_s", &njet_s, "njet_s/I");
    tree->Branch("e_jet", e_jet, "e_jet[njet]/F");
    tree->Branch("pt_jet", pt_jet, "pt_jet[njet]/F");
    tree->Branch("eta_jet", eta_jet, "eta_jet[njet]/F");
    tree->Branch("phi_jet", phi_jet, "phi_jet[njet]/F");
    tree->Branch("p_jet", p_jet, "p_jet[njet]/F");
    tree->Branch("theta_jet", theta_jet, "theta_jet[njet]/F");
    tree->Branch("e_jet_s", e_jet_s, "e_jet_s[njet_s]/F");
    tree->Branch("pt_jet_s", pt_jet_s, "pt_jet_s[njet_s]/F");
    tree->Branch("eta_jet_s", eta_jet_s, "eta_jet_s[njet_s]/F");
    tree->Branch("phi_jet_s", phi_jet_s, "phi_jet_s[njet_s]/F");
    tree->Branch("p_jet_s", p_jet_s, "p_jet_s[njet_s]/F");
    tree->Branch("theta_jet_s", theta_jet_s, "theta_jet_s[njet_s]/F");

    // New TH1D
    const int nbins_pt{100};
    const double lobin_pt{0.};
    const double hibin_pt{10.};
    TH1D *photon_dee{new TH1D("photon dE/E", "photon dE/E;dE/E; count", nbins_pt, -2, 2)};
    TH1D *pion_dee{new TH1D("pi- dE/E", "pion dE/E;dE/E; count", nbins_pt, -2, 2)};

    // Select FastJet parameters
    const double R_FULL = 0.4; // Jet size.
    const double R_CH = 0.4;   // Jet size.
    const double pTMin = 0.2;  // Min jet pT.
    const double etaMax{4.0};
    const double etaMin{2.5};
    fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, R_FULL);
    fastjet::Selector jetrap = fastjet::SelectorEtaMin(etaMin + R_FULL) * fastjet::SelectorEtaMax(etaMax - R_FULL); // for full jets

    // Initialize Pythia
    Pythia pythia;
    const double sNN{200};
    const double pTHatMin{1000.};
    const double pTHatMax{2000.};
    string name_type = "pp";
    int idA;
    int idB;
    TRandom3 rand{0};
    TRandom *gaus = new TRandom();

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
             //Form("PhaseSpace:pTHatMin = %f",pTHatMin),
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

    int n_events{10000};
    const double mips_min_p{0.2};

    for (int iEvent{0}; iEvent < n_events; ++iEvent)
    {
        if (!pythia.next())
            continue;

        Event &event = pythia.event;
        std::vector<fastjet::PseudoJet> part_FULL;
        std::vector<fastjet::PseudoJet> part_SMEAR;
        std::vector<fastjet::PseudoJet> part_CH;
        float max_bemc_Et = 0.;

        evid = iEvent;
        xsec = pythia.info.sigmaGen();
        ntrials = pythia.info.nTried();
 
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

            PseudoJet p{e.px(), e.py(), e.pz(), e.pAbs()};
            part_FULL.push_back(p);

            double ratio, es, pxs, pys, pzs, ps;
            int id{e.id()};

            // smear by resolution of EMCAL and HCAL
            if (id == 111 || id == 22 || (id > 11 && id < 18) || (id < -11 && id > -18))
            {
                es = gaus->Gaus(e.e(), 0.1 * (sqrt(e.e())));
            }
            else
            {
                es = gaus->Gaus(e.e(), 0.5 * sqrt(e.e()) + 0.1 * e.e());
            }

            if (es < e.m0() || es < 0.2)
                continue;
            ps = sqrt(pow(es, 2) - pow(e.m0(), 2));
            if (ps < 0.2)
                continue;

            ratio = ps / e.pAbs();
            if (id == 22)
            {
                //cout << "photon, " << "true e: " << e.e() << ", smeared e: " << es << ", delta e: " << es-e.e() << endl;
                photon_dee->Fill((es - e.e()) / e.e());
            }
            else if (id == -211)
            {
                //cout << "pi-, " << "true e: " << e.e() << ", smeared e: " << es << ", delta e: " << es-e.e() << endl;
                pion_dee->Fill((es - e.e()) / e.e());
            }

            PseudoJet psmeared{e.px() * ratio, e.py() * ratio, e.pz() * ratio, ps}; // very lazy way to get psuedorapidity
            part_SMEAR.push_back(psmeared);
        }

        // make and use full jets
        fastjet::ClusterSequence clustSeq_FULL(part_FULL, jetdef);
        vector<fastjet::PseudoJet> jetsFULL = sorted_by_rapidity(jetrap(clustSeq_FULL.inclusive_jets(pTMin)));
        njet = jetsFULL.size();
	for (unsigned int ij = 0; ij < jetsFULL.size(); ij++)
        {
            e_jet[ij] = jetsFULL[ij].e();
            pt_jet[ij] = jetsFULL[ij].perp();
            eta_jet[ij] = jetsFULL[ij].eta();
            phi_jet[ij] = jetsFULL[ij].phi();
            p_jet[ij] = jetsFULL[ij].perp() * TMath::CosH(jetsFULL[ij].eta());
            theta_jet[ij] = jetsFULL[ij].theta();
        }

        // make and use smeared jets
        fastjet::ClusterSequence clustSeq_SMEAR(part_SMEAR, jetdef);
        vector<fastjet::PseudoJet> jetsSMEAR = sorted_by_rapidity(jetrap(clustSeq_SMEAR.inclusive_jets(pTMin)));
        njet_s = jetsSMEAR.size();
	for (unsigned int ij = 0; ij < jetsSMEAR.size(); ij++)
        {
            e_jet_s[ij] = jetsSMEAR[ij].e();
            pt_jet_s[ij] = jetsSMEAR[ij].perp();
            eta_jet_s[ij] = jetsSMEAR[ij].eta();
            phi_jet_s[ij] = jetsSMEAR[ij].phi();
            p_jet_s[ij] = jetsSMEAR[ij].perp() * TMath::CosH(jetsSMEAR[ij].eta());
            theta_jet_s[ij] = jetsSMEAR[ij].theta();

        }

        // make and use charged jets
        fastjet::ClusterSequence clustSeq_CH(part_CH, jetdef);
        vector<fastjet::PseudoJet> jetsCH = sorted_by_pt(jetrap(clustSeq_CH.inclusive_jets(pTMin)));
        // do something with these jets...
        
	tree -> Fill();
  
    }
    
    cout << " finished PYTHIA loop " << endl;

    f_out->Write();
    f_out->Save();
    f_out->Close();

    return 0;
};
