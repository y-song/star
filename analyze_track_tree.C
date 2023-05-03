void analyze_track_tree()
{
    // create pT bins
    const int nbin = 16;
    const double pt_min = 2.0;
    const double pt_max = 10.0;
    const double pt_interval = 0.5;
    double pt_bin[nbin + 1];
    for (int ibin = 0; ibin < nbin + 1; ibin++)
    {
        pt_bin[ibin] = pt_min + ibin * pt_interval;
    }

    // set up histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(11);

    TObjArray partpid;
    
    for (int ibin = 0; ibin < nbin; ibin++)
    {
        partpid.Add(new TH1F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";fraction;particle PID", 6000, -3000, 3000));
    }

    /*TH1F *h_trackm2 = new TH1F("h_trackm2", ";track M^{2} [GeV^{2}]", 50, -0.5, 1.5);
    TH1F *h_nspi = new TH1F("h_tracknspi", ";track nsigmapi", 50, -5, 5);
    TH1F *h_nsk = new TH1F("h_tracknsk", ";track nsigmak", 50, -5, 5);
    TH1F *h_nsp = new TH1F("h_tracknsp", ";track nsigmap", 50, -5, 5);

    TH2F *h_pt_nsknspidiff = new TH2F("h_pt_nsknspidiff", ";track p_{T} [GeV];track (nsigmaK - nsigmapi)", nbin, pt_min, pt_max, 60, 0, 3);
    TH2F *h_pt_m2 = new TH2F("h_pt_m2", ";track p_{T} [GeV];track M^{2} [GeV^{2}]", 2 * nbin, pt_min, pt_max / 2, 100, -0.5, 1.5);
    TH2F *h_trackm2_tracknspi = new TH2F("h_trackm2_tracknspi", ";track M^{2} [GeV^{2}];track nsigmapi", 75, -1, 2, 100, -5, 5);
    TH2F *h_pt_dedx = new TH2F("h_pt_dedx", ";p_{T} [GeV];dE/dx [keV/cm]", 300, 2, 17, 100, 0, 10);*/

    // TChain *t = new TChain("tracks");
    // t->Add("/gpfs01/star/pwg/youqi/run15/result/tracks_0428/E4*.root");
    TFile *f = new TFile("test.root");
    TTree *t = (TTree *)f->Get("tracks");

    Double_t pt, eta, phi, m2;
    Int_t pid, q;

    t->SetBranchAddress("pt", &pt);
    t->SetBranchAddress("eta", &eta);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("m2", &m2);
    t->SetBranchAddress("pid", &pid);
    t->SetBranchAddress("q", &q);

    int ntrack = t->GetEntries();
    for (int itrack = 0; itrack < ntrack; itrack++)
    {
        t->GetEntry(itrack);
        for (int ibin = 0; ibin < nbin; ibin++)
        {
            if (pt_bin[ibin] < pt && pt < pt_bin[ibin + 1])
            {
                ((TH1F *)partpid.At(ibin))->Fill(pid);
            }
        }
    }

    TCanvas *c = new TCanvas("c", "c");

    /*h_nspi->Draw();
    h_nsk->Draw();
    h_nsp->Draw();
    h_pt_dedx->Draw("colz");
    gPad->SetLogy(1);
    h_trackm2->Draw();
    h_pt_nsknspidiff->Draw("colz");
    h_pt_nsknspidiff->ProfileX()->Draw("same");
    h_trackm2_tracknspi->Draw("colz");
    h_pt_m2->Draw("colz");
    c->SaveAs("plots/pt_m2.png", "png");*/

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        auto dist = ((TH1F *)partpid.At(ibin));
        Double_t integral = dist->GetEntries();
        dist->Scale(1./integral);
        dist->Draw();
        dist->GetYaxis()->SetRangeUser(0., 1.);
        c->SaveAs(Form("plots/pid_%0.1f_pT_%0.1f.png", pt_bin[ibin], pt_bin[ibin + 1]), "png");
    }
}