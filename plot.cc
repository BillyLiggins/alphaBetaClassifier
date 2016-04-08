#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>

#include <TGraph2D.h>
#include <string>
#include <TMultiGraph.h>
#include <TH2D.h>
#include <TFrame.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TPad.h>


#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <fstream>


string partflag="";

vector<string> glob( const string& path, const string& start )
{
  vector<string> result;
  TSystemDirectory dir(path.c_str(), path.c_str());
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith( start ) ) {
        stringstream fullPath; fullPath << path << "/" << fname.Data();
        result.push_back(fullPath.str());
      }
    }
  }
  return result;
}


void FillHist(TFile* file,TH1D * hist,TH2D * hist2){

	TTree* Tree = (TTree*) file->Get("output");
	Double_t para, mcEdepQuenched;
	Bool_t  Qfit;
	Int_t pdg1, pdg2;
        Int_t parentpdg1,parentpdg2;

        Tree->SetBranchAddress("pdg1",&pdg1);
        Tree->SetBranchAddress("pdg2",&pdg2);
        Tree->SetBranchAddress("parentpdg1",&parentpdg1);
        Tree->SetBranchAddress("parentpdg2",&parentpdg2);
	
	Tree->SetBranchAddress("berkeleyAlphaBeta",&para);
	Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
	Tree->SetBranchAddress("fitValid",&Qfit);
	Int_t n = (Int_t)Tree->GetEntries();
	Int_t code;

	for( Int_t i =0;i<n;i++){
		Tree->GetEntry(i);
		if (partflag=="Bi"){
			code=11;
		}else if( partflag=="Po"){
			code=1000020040;
		}
		if( Qfit && pdg1==code  ){
			hist->Fill(para);
			hist2->Fill(mcEdepQuenched,para);
		}
	}

}
int main(){
	
	gStyle->SetOptStat(0);
	TH1D *hBi210   = new TH1D("hBi210","berkeleyAlphaBeta",100,-100,100);
	hBi210->SetLineColor(kBlue);hBi210->SetLineWidth(3);
	hBi210->SetFillColor(kBlue);
	TH1D *hPo210   = new TH1D("hPo210","berkeleyAlphaBeta",100,-100,100);
	hPo210->SetLineColor(kRed);hPo210->SetLineWidth(3);
	TH1D *hBlank   = new TH1D("hBlank","berkeleyAlphaBeta",100,-100,100);

	TH2D* compareBi210   = new TH2D("compareBi210","berkeleyAlphaBeta",23,0,2.3,260,-160,100);
	compareBi210->SetLineColor(kBlue);compareBi210->SetLineWidth(3);
	compareBi210->SetMarkerColor(kBlue);
	compareBi210->SetFillColor(kBlue);
	TH2D* comparePo210   = new TH2D("comparePo210","berkeleyAlphaBeta",23,0,2.3,260,-160,100);
	comparePo210->SetLineColor(kRed);comparePo210->SetLineWidth(3);
	comparePo210->SetMarkerColor(kRed);
	comparePo210->SetFillColor(kRed);
	//TH2D *hBlank2   = new TH2D("hBlank2","berkeleyAlphaBeta",100,-100,100);
	
	TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
	TLegend* t2 = new TLegend( 0.11, 0.11, 0.21, 0.21 );
	
	//===================================Bi210 & Po210=================================================
	//----------------These are the classifiers from the production 5_3 files---------------------------
	//vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi210","SolarBi210");
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po210","SolarPo210");
	//----------------Files Recoordinated with 1 MeV e- and 5 MeV alpha---------------------------
	//vector<string> biFileList= glob("output","Bi");
	//vector<string> poFileList= glob("output","Po");
	//----------------Files Recoordinated with 5 MeV e- and 50 MeV alpha---------------------------
	//vector<string> biFileList= glob("output2","Bi");
	//vector<string> poFileList= glob("output2","Po");

	//===================================Bi212 & Po212=================================================
	//----------------These are the classifiers from the production 5_3 files---------------------------
	//vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi212","SolarBi212");
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po212","SolarPo212");
	//----------------Files Recoordinated with 1 MeV e- and 5 MeV alpha---------------------------
	//vector<string> biFileList= glob("output","Bi");
	//vector<string> poFileList= glob("output","Po");
	//----------------Files Recoordinated with 5 MeV e- and 50 MeV alpha---------------------------
	//vector<string> biFileList= glob("output2","Bi");
	//vector<string> poFileList= glob("output2","Po");

	//===================================Bi214 & Po214=================================================
	//----------------These are the classifiers from the production 5_3 files---------------------------
	//vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi214","SolarBi214");
	//vector<string> poFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Po214","SolarPo214");
	//----------------Files Recoordinated with 1 MeV e- and 5 MeV alpha---------------------------
	//vector<string> biFileList= glob("output","Bi");
	//vector<string> poFileList= glob("output","Po");
	//----------------Files Recoordinated with 5 MeV e- and 50 MeV alpha---------------------------
	//vector<string> biFileList= glob("output2","Bi");
	//vector<string> poFileList= glob("output2","Po");
	//===================================Bi210 & Alphas=================================================
	//----------------These are the classifiers from the production 5_3 files---------------------------
	vector<string> biFileList= glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi212","SolarBi212");
	vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

	for( int i=0; i<biFileList.size(); i++ ){
		TFile * file= TFile::Open(biFileList[i].c_str());	
		partflag="Bi";
		FillHist( file, hBi210,compareBi210);
	  }

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		partflag="Po";
		FillHist( file, hPo210,comparePo210);
	  }

	TCanvas * c1 = new TCanvas();
	hBlank->SetMaximum(0.3);
	hBlank->SetTitle("Standard Classifiers From Production 5.3");
	//hBlank->SetTitle("Classifers Recoordinated with 1 MeV #beta & 5 MeV #alpha");
	//hBlank->SetTitle("Classifers Recoordinated with 5 MeV #beta & 50 MeV #alpha");
	hBlank->Draw();
	cout<<"# of entries in hBi210 = "<< hBi210->GetEntries()<<endl;
	cout<<"# of entries in hPo210 = "<< hPo210->GetEntries()<<endl;
	hBi210->DrawNormalized("same");
	hPo210->DrawNormalized("same");

	t1->AddEntry( hBi210, "Bi 212","f");
	t1->AddEntry( hPo210, "Po 212","f");
	t1->Draw();
	//c1->Print("BerkeleyAlphaBeta_NormalWithRAT5_3.png");
	//c1->Print("output2.png");
	//c1->Print("production_5_3.png");
	c1->Print("StandardClassifersFromProduction5_3_212.png");
	//c1->Print("StandardClassifersFromProduction5_3.png");
	//c1->Print("RecoodinatedClassifiers_1MeV_beta_5MeV_alpha.png");
	//c1->Print("RecoodinatedClassifiers_5MeV_beta_50MeV_alpha.png");

	TCanvas * c2 = new TCanvas();
	cout<<"# of entries in hBi210 = "<< hBi210->GetEntries()<<endl;
	cout<<"# of entries in hPo210 = "<< hPo210->GetEntries()<<endl;
	compareBi210->SetTitle("Standard Classifiers From Production 5.3");
	compareBi210->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
	compareBi210->GetYaxis()->SetTitle("BerkeleyAlphaBeta");
	compareBi210->Draw();
	comparePo210->Draw("same");

	t2->AddEntry( compareBi210, "Bi 210","f");
	t2->AddEntry( comparePo210, "Po 210","f");
	t2->Draw();
	c2->Print("compareBerkeleyAlphaBetaVsMCEnergy.png");


	TFile fileout("plot.root","RECREATE");
	fileout.cd();
	hBi210->Write();
	hPo210->Write();
	fileout.Close();
	return 0;

}

