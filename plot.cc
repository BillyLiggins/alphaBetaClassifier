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


void FillHist(TFile* file,TH1D * hist){

	TTree* Tree = (TTree*) file->Get("output");
	Double_t para;

	Tree->SetBranchAddress("berkeleyAlphaBeta",&para);
	Int_t n = (Int_t)Tree->GetEntries();

	for( Int_t i =0;i<n;i++){
	Tree->GetEntry(i);	
	hist->Fill(para);
	}

}
int main(){
	
	gStyle->SetOptStat(0);
	TH1D *hBi210   = new TH1D("hBi210","berkeleyAlphaBeta",100,-50,50);
	hBi210->SetLineColor(kBlue);hBi210->SetLineWidth(3);
	hBi210->SetFillColor(kBlue);
	TH1D *hPo210   = new TH1D("hPo210","berkeleyAlphaBeta",100,-50,50);
	hPo210->SetLineColor(kRed);hPo210->SetLineWidth(3);
	
	TLegend* t1 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
	vector<string> biFileList= glob("output","Bi");
	vector<string> poFileList= glob("output","Po");
	for( int i=0; i<biFileList.size(); i++ ){
		TFile * file= TFile::Open(biFileList[i].c_str());	
		FillHist( file, hBi210);
	  }

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		FillHist( file, hPo210);
	  }

	hBi210->Draw();
	hPo210->Draw("same");

	t1->AddEntry( hBi210, "Bi 210","f");
	t1->AddEntry( hPo210, "Po 210","f");
	t1->Draw();

	TFile fileout("plot.root","RECREATE");
	fileout.cd();
	hBi210->Write();
	hPo210->Write();
	fileout.Close();
	return 0;

}
