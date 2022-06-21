#include <iostream>
#include<string>
#include<vector>

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TSystem.h"

void leptonicTop() {
  TH1::AddDirectory(kFALSE);
	gROOT->Reset();

    TFile *f1 = new TFile("file.root");

    TH1F* histo = (TH1F*) f1->Get("hist1");

TCanvas *c1 = new TCanvas("c1", "c14", 800,1000);
histo->Fit("gaus","","",100,230.);
histo->GetYaxis()->SetTitleOffset(1.55);
gStyle->SetOptFit(0111);
histo->Draw();

c1->Print("leptonicT.png");


}
