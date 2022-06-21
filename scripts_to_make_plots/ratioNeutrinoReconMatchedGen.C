#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TLegend.h"

void ratioNeutrinoReconMatchedGen(){
  TH1::AddDirectory(kFALSE);
	gROOT->Reset();

	TFile* f = new TFile("nanoLatino_AZH__part0.root", "READ");
	TTree* tSelected = (TTree*)f->Get("Events");

  Float_t AZH_GeneratorMatchedNeutrinoPZ = 0.;
  Float_t AZH_GeneratorNeutrinoPZ = 0.;

  tSelected->SetBranchAddress("AZH_GeneratorMatchedNeutrinoPZ", &AZH_GeneratorMatchedNeutrinoPZ);
  tSelected->SetBranchAddress("AZH_GeneratorNeutrinoPZ", &AZH_GeneratorNeutrinoPZ);


  TCanvas *c = new TCanvas("c", "canvas", 800, 800);

  TH1F* histo1 = new TH1F("NeutrinoPZ", "NeutrinoPZ", 40,-400, 400);
  TH1F* histo2 = new TH1F("NeutrinoPZ", "NeutrinoPZ", 40,-400, 400);


  TH1F* hratio = new TH1F("AZH_GeneratorNeutrinoPZ", "AZH_GeneratorNeutrinoPZ", 40,-400, 400);


  int nEntries = tSelected->GetEntries();

  for (int i = 0; i < nEntries; i++) {
    tSelected ->GetEntry(i);

    if (i%1000) std::cout << "reading event: " << i << std::endl;
          if (AZH_GeneratorMatchedNeutrinoPZ != -9999) {
              histo1->Fill(AZH_GeneratorMatchedNeutrinoPZ);
              histo2->Fill(AZH_GeneratorNeutrinoPZ);
            }
      }


  histo1->Scale(1/histo1->Integral());
  histo2->Scale(1/histo2->Integral());

  TLegend* legend = new TLegend(0.62,0.68,0.86,0.85);
  legend->SetTextSize(0.020);
  legend->AddEntry(histo1, "gen_matched ", "f");
  legend->AddEntry(histo2, "true ", "f");

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  histo2->SetStats(0);          // No statistics on upper plot
  histo2->Draw("hist");               // Draw h1
  histo1->Draw("hist same");         // Draw h2 on top of h1
  legend->Draw();

  //BOTTOM RATIO PLOT
  c->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();

  TH1F *h3 = (TH1F*)histo1->Clone("h3");
  h3->Divide(histo2);

  h3->SetLineColor(kBlack);
  h3->SetMinimum(-1);  // Define Y ..
  h3->SetMaximum(5); // .. range
  h3->Sumw2();
  h3->SetStats(0);      // No statistics on lower plot
  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  h3->Draw();

  histo1->SetLineColor(kBlue+1);
  histo1->SetLineWidth(3);

  // Y axis h1 plot settings
  histo1->GetYaxis()->SetTitleSize(20);
  histo1->GetYaxis()->SetTitleFont(43);
  histo1->GetYaxis()->SetTitleOffset(1.55);

  // h2 settings
  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(3);

  // Ratio plot (h3) settings
  h3->SetTitle(""); // Remove the ratio title

  // Y axis ratio plot settings
  h3->GetYaxis()->SetTitle("ratio gen_matched/true ");
  h3->GetYaxis()->SetNdivisions(505);
  h3->GetYaxis()->SetTitleSize(20);
  h3->GetYaxis()->SetTitleFont(43);
  h3->GetYaxis()->SetTitleOffset(1.55);
  h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetYaxis()->SetLabelSize(15);

  // X axis ratio plot settings
  h3->GetXaxis()->SetTitle("neutrino pz (GeV)");
  h3->GetXaxis()->SetTitleSize(20);
  h3->GetXaxis()->SetTitleFont(43);
  h3->GetXaxis()->SetTitleOffset(3.05);
  h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetXaxis()->SetLabelSize(15);

  c->Update();
  c->Print("ratioNeutrinoReconMatchedGen.png");
}
