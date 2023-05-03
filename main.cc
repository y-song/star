#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
//#include "TVectorD.h"
//#include "TStopwatch.h"
//#include "TClonesArray.h"
#include "TObject.h"
#include <array>

#include "TRandom3.h"
//#include <sstream>
#include <iostream>
#include <cstring>

#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
/*#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "Pythia8Plugins/FastJet3.h"
#include "Math/ProbFunc.h"*/

using namespace Pythia8;
//using namespace fastjet;
using namespace std;

int main()
{
    // +++++++++++++++++++++++++++++++++++++ set up output file +++++++++++++++++++++++++++++++++++++
    TFile out("test.root", "RECREATE");
    out.cd();

    // +++++++++++++++++++++++++++++++++++++ set up histograms and track tree +++++++++++++++++++++++++++++++++++++
    Double_t pt, eta, phi, m2;
    Int_t pid, q;

    TTree *t = new TTree("tracks", "tracks");
    t->Branch("pt", &pt, "pt/D");
    t->Branch("eta", &eta, "eta/D");
    t->Branch("phi", &phi, "phi/D");
    t->Branch("m2", &m2, "m2/D");
    t->Branch("pid", &pid, "pid/I");
    t->Branch("q", &q, "q/I");

    // +++++++++++++++++++++++++++++++++++++ initialize Pythia +++++++++++++++++++++++++++++++++++++
    Pythia pythia;
    const double sNN{200};
    const double pTHatMin{0};
    // const double pTHatMax{20.};
    int idA = 2212;
    int idB = 2212;
    TRandom3 rand{0};

    // to do: implement config from rhig308 files, adding in Detroit tune parameters and turning off weak decays
    for (auto str : vector<string>{
             Form("Beams:eCM = %f", sNN),
             "HardQCD:all = on",
             Form("PhaseSpace:pTHatMin = %f", pTHatMin),
             // Form("PhaseSpace:pTHatMax = %f",pTHatMax),
             Form("Beams:idA = %i", idA), // moving in +z direction
             Form("Beams:idB = %i", idB),
             "Random:setSeed = on",
             Form("Random:seed = %i", rand.Integer(10000000)),
             "Next:numberCount = 10000",
         })

        pythia.readString(str.c_str());
    pythia.init();
    Pythia8::Info &info = pythia.info;

    cout << " Starting the pythia loop. " << endl;

    // +++++++++++++++++++++++++++++++++++++ event loop +++++++++++++++++++++++++++++++++++++
    int n_events{20000};
    for (int iEvent{0}; iEvent < n_events; ++iEvent)
    {
        if (!pythia.next())
            continue;
        Event &event = pythia.event;

        // +++++++++++++++++++++++++++++++++++++ particle loop +++++++++++++++++++++++++++++++++++++
        for (int i{0}; i < event.size(); ++i)
        {
            auto &e{event[i]};
            if (!e.isFinal())
                continue;
            if (!e.isVisible())
                continue;
            if (e.pT() < 2.0)
                continue;
            if (e.pT() > 30.0)
                continue;
            if (e.eta() < -1.0 || e.eta() > 1.0)
                continue;

            pt = e.pT();
            eta = e.eta();
            phi = e.phi();
            m2 = e.m()*e.m();
            pid = e.id();
            q = e.chargeType();

            t->Fill();
        } // end particle loop
    } // end event loop

    // +++++++++++++++++++++++++++++++++++++ save output +++++++++++++++++++++++++++++++++++++
    out.Write();
    out.Save();
    out.Close();

    return 0;
};