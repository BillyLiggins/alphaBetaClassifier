// This file has a main function that goes and find cuts for the berkeleyAlphaBeta across energy and radius (independently). 
// You were trying to get a radial fit done.
// TODOs:
//
// 	1) combine the radial and energy varibles in to one TH3D 
// 	2) Clean it 
//

#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>

#include <TGraph2D.h>
#include <TMath.h>
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
#include <numeric>

#define SSTR( x ) static_cast< std::ostringstream & >( \
                     ( std::ostringstream() << std::dec << x ) ).str()


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


void FillHist(TFile* file,TH1D * hist){

	TTree* Tree = (TTree*) file->Get("output");
	Double_t para, mcEdepQuenched,posr;
	Bool_t  Qfit;
	Int_t pdg1, pdg2;
        Int_t parentpdg1,parentpdg2;

        Tree->SetBranchAddress("pdg1",&pdg1);
        Tree->SetBranchAddress("pdg2",&pdg2);
        Tree->SetBranchAddress("parentpdg1",&parentpdg1);
        Tree->SetBranchAddress("parentpdg2",&parentpdg2);
	
	Tree->SetBranchAddress("berkeleyAlphaBeta",&para);
	Tree->SetBranchAddress("mcEdepQuenched",&mcEdepQuenched);
	Tree->SetBranchAddress("posr",&posr);
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
			hist->Fill(mcEdepQuenched);
		}
	}

}


void findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values){
	/* This function should take a 2d histrogram split it into energy 
	 * strips and then find a cut value which retain ~99% of the signal.
	 */
	//for(double energy =0; energy<3.5;energy+=0.1){

	double startingCut=-1000;
	double step=0.1;
	double cutValue=startingCut;
	double accept=0;
	double mistagged=0;
	double rejection=0;
	double minBin=compareBi210->GetXaxis()->GetXmin();
	double maxBin=compareBi210->GetXaxis()->GetXmax();
	double sliceWidth=compareBi210->GetXaxis()->GetBinWidth(1);
	cout<<"minBin = "<<minBin<<endl;
	cout<<"maxBin = "<<maxBin<<endl;
	cout<<"sliceWidth = "<<sliceWidth<<endl;
	double threshold=0.995;
	TAxis* bi_x=compareBi210->GetXaxis();
	TAxis* bi_y=compareBi210->GetYaxis();
	TAxis* po_x=comparePo210->GetXaxis();
	TAxis* po_y=comparePo210->GetYaxis();

	TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(4.));
	TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(4.));
	double bi_numEV_complete= compareBi210->GetEntries();
	double po_numEV_complete= comparePo210->GetEntries();
	cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl;
	cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl;

	for(double energy =minBin; energy<maxBin;energy+=sliceWidth){
	/* for(double energy =0; energy<2.5;energy+=0.1){ */
		cout<<"Energy = "<< energy<<endl;

		//These histograms contain the sliced projections
		TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth));
		TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth));

		//These two numbers are the num of events before cuts
		double bi_numEV=slice_bi->GetEntries();
		double po_numEV=slice_po->GetEntries();

		TAxis* bi_slice_x=slice_bi->GetXaxis();
		TAxis* po_slice_x=slice_po->GetXaxis();

		while(accept<threshold){

			//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step'
			double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
			//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
			cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl;
			//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
			accept=bi_remain/bi_numEV;	
			cout<<"Accepted fraction = "<<accept<<endl;
		}
		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));

		rejection= po_numEV/mistagged;
		cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl;
		Rejection_values.push_back(rejection);
		cutValues.push_back(cutValue);
		energyValues.push_back(energy+ sliceWidth/2);

		if(false){
		TCanvas* c1 = new TCanvas();
		c1->cd();
		//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
		TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
		cutLine->SetLineWidth(2);cutLine->SetLineStyle(4);
		slice_bi->Draw();
		slice_po->Draw("same");
		cutLine->Draw("same");
		}
		cutValue=startingCut;
		accept=0;
	}

}

void findCutsRadial(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& rejection_values){
	/* This function should take a 2d histrogram split it into energy 
	 * strips and then find a cut value which retain ~99% of the signal.
	 */

	double startingCut=-1000;
	double step=0.1;
	double step_radial=500;
	double cutValue=startingCut;
	double accept=0;
	double mistagged=0;
	double rejection=0;
	double sliceWidth=100;
	double threshold=0.9;
	TAxis* bi_x=compareBi210->GetXaxis();
	TAxis* bi_y=compareBi210->GetYaxis();
	TAxis* po_x=comparePo210->GetXaxis();
	TAxis* po_y=comparePo210->GetYaxis();

	TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(8000.));
	TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(8000.));
	double bi_numEV_complete= compareBi210->GetEntries();
	double po_numEV_complete= comparePo210->GetEntries();
	cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl;
	cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl;

	for(double energy =0; energy<6000;energy+=step_radial){
		cout<<"radius = "<< energy<<endl;

		//TH1D* slice_bi=compareBi210->ProjectionY("slice_bi",bi_x->FindBin(0.5),bi_x->FindBin(1.0));
		TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth));
		TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth));

		double bi_numEV=slice_bi->GetEntries();
		double po_numEV=slice_po->GetEntries();
		cout<<"Number of entries in bi slice = "<< bi_numEV<<endl;
		cout<<"Number of entries in po slice = "<< po_numEV<<endl;
		TAxis* bi_slice_x=slice_bi->GetXaxis();
		TAxis* po_slice_x=slice_po->GetXaxis();

		while(accept<threshold){

			double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step));
			//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step));
			//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl;
			accept=bi_remain/bi_numEV;	
		}
		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue));

		//This rejection is a cheat to get around div0.
		if(po_numEV>1 && mistagged>1){
		rejection= po_numEV/(mistagged+0.000001);
		}else{
		rejection= 0;
		}
		cout<<mistagged<<" events were mistagged. Rejection =  "<<rejection<<"."<<endl;
		rejection_values.push_back(rejection);
		cutValues.push_back(cutValue);
		energyValues.push_back(energy+step_radial/2);
		if(false){
		TCanvas* c1 = new TCanvas();
		c1->cd();
		//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() );
		TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() );
		cutLine->SetLineWidth(2);cutLine->SetLineStyle(4);
		slice_bi->Draw();
		slice_po->Draw("same");
		cutLine->Draw("same");
		}
		cutValue=startingCut;
		accept=0;
	}

}

int main(){
//====================================================================	
	gStyle->SetOptStat(0);

	TH1D *hist_blanck_beta= new TH1D("hist_blanck_beta","mcEdepQuenched",1000,0,2.6);
	TH1D *hist_blanck_alpha= new TH1D("hist_blanck_alpha","mcEdepQuenched",1000,0,2.6);
	
	TLegend* t2 = new TLegend( 0.6, 0.7, 0.89, 0.88 );
	

	vector<string>  biFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
	vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

	TCanvas * c2 = new TCanvas();
	c2->cd();
	hist_blanck_beta->SetTitle("#beta across mcEdepQuenched");
	hist_blanck_beta->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");

	hist_blanck_beta->Draw();

	for( int i=0; i<biFileList.size(); i++ ){
		TH1D *hist_beta= new TH1D("hist_beta","#beta across mcEdepQuenched",200,0,2.6);
		hist_beta->SetLineColor(kBlue);//hist_beta->SetLineWidth(3);
		/* hist_beta->SetFillColor(kBlue); */

		TFile * file= TFile::Open(biFileList[i].c_str());	
		partflag="Bi";
		FillHist( file,hist_beta);
		if(i==0) hist_beta->Draw();
		hist_beta->Draw("same");
	  }

	c2->Print("betaEnergies.png");

	TCanvas * c3 = new TCanvas();
	c3->cd();
	hist_blanck_alpha->SetTitle("#alpha across mcEdepQuenched");
	hist_blanck_alpha->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");

	hist_blanck_alpha->Draw();
	for( int i=0; i<poFileList.size(); i++ ){
		TH1D *hist_alpha= new TH1D("hist_alpha","#alpha across mcEdepQuenched",200,0,2.6);
		hist_alpha->SetLineColor(kBlue);
		/* hist_alpha->SetLineWidth(3); */
		/* hist_alpha->SetFillColor(kBlue); */

		TFile * file= TFile::Open(poFileList[i].c_str());	
		partflag="Po";
		FillHist( file,hist_alpha);
		if(i==0) hist_alpha->Draw();
		hist_alpha->Draw("same");
	  }

	c3->Print("alphaEnergies.png");


	return 0;

}

