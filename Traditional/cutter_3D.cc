#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>

#include <TGraph2D.h>
#include <TMath.h>
#include <string>
#include <TMultiGraph.h>
#include <TH2D.h>
#include <TH3D.h>
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


void FillHist(TFile* file,TH3D* hist){

	TTree* Tree = (TTree*) file->Get("output");
	Double_t berkeleyAlphaBeta, mcEdepQuenched,posr;
	Bool_t  Qfit;
	Int_t pdg1, pdg2;
        Int_t parentpdg1,parentpdg2;

        Tree->SetBranchAddress("pdg1",&pdg1);
        Tree->SetBranchAddress("pdg2",&pdg2);
        Tree->SetBranchAddress("parentpdg1",&parentpdg1);
        Tree->SetBranchAddress("parentpdg2",&parentpdg2);
	
	Tree->SetBranchAddress("berkeleyAlphaBeta",&berkeleyAlphaBeta);
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
			// x-axis- MC energy
			// y-axis- posr
			// z-axis- berkeleyAlphaBeta
			hist->Fill(mcEdepQuenched,posr,berkeleyAlphaBeta);
		}
	}

}


/* void findCutsEnergy(TH3D *ele_hist,TH3D *alpha_hist,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values){ */
/* 	/1* This function should take a 2d histrogram split it into energy */ 
/* 	 * strips and then find a cut value which retain ~99% of the signal. */
/* 	 *1/ */
/* 	//for(double energy =0; energy<3.5;energy+=0.1){ */

/* 	double startingCut=-1000; */
/* 	double step=0.1; */
/* 	double cutValue=startingCut; */
/* 	double accept=0; */
/* 	double mistagged=0; */
/* 	double rejection=0; */
/* 	double sliceWidth=0.1; */
/* 	double threshold=0.995; */
/* 	TAxis* bi_x=compareBi210->GetXaxis(); */
/* 	TAxis* bi_y=compareBi210->GetYaxis(); */
/* 	TAxis* po_x=comparePo210->GetXaxis(); */
/* 	TAxis* po_y=comparePo210->GetYaxis(); */

/* 	TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(4.)); */
/* 	TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(4.)); */
/* 	double bi_numEV_complete= compareBi210->GetEntries(); */
/* 	double po_numEV_complete= comparePo210->GetEntries(); */
/* 	cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl; */
/* 	cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl; */

/* 	for(double energy =0; energy<2.5;energy+=0.1){ */
/* 		cout<<"Energy = "<< energy<<endl; */

/* 		//These histograms contain the sliced projections */
/* 		TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth)); */
/* 		TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth)); */

/* 		//These two numbers are the num of events before cuts */
/* 		double bi_numEV=slice_bi->GetEntries(); */
/* 		double po_numEV=slice_po->GetEntries(); */

/* 		TAxis* bi_slice_x=slice_bi->GetXaxis(); */
/* 		TAxis* po_slice_x=slice_po->GetXaxis(); */

/* 		while(accept<threshold){ */

/* 			//This loop finds the cut value to a certain threshold is achevied, it does so by incrementing by 'step' */
/* 			double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step)); */
/* 			//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step)); */
/* 			cout<<"Number of entries remaining in bi slice = "<< bi_remain<<endl; */
/* 			//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl; */
/* 			accept=bi_remain/bi_numEV; */	
/* 			cout<<"Accepted fraction = "<<accept<<endl; */
/* 		} */
/* 		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue)); */

/* 		rejection= po_numEV/mistagged; */
/* 		cout<<mistagged<<" events were mistagged. Rejection= "<<rejection<<"."<<endl; */
/* 		Rejection_values.push_back(rejection); */
/* 		cutValues.push_back(cutValue); */
/* 		energyValues.push_back(energy+0.05); */

/* 		if(false){ */
/* 		TCanvas* c1 = new TCanvas(); */
/* 		c1->cd(); */
/* 		//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() ); */
/* 		TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() ); */
/* 		cutLine->SetLineWidth(2);cutLine->SetLineStyle(4); */
/* 		slice_bi->Draw(); */
/* 		slice_po->Draw("same"); */
/* 		cutLine->Draw("same"); */
/* 		} */
/* 		cutValue=startingCut; */
/* 		accept=0; */
/* 	} */

/* } */

/* void findCutsRadial(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& rejection_values){ */
/* 	/1* This function should take a 2d histrogram split it into energy */ 
/* 	 * strips and then find a cut value which retain ~99% of the signal. */
/* 	 *1/ */

/* 	double startingCut=-1000; */
/* 	double step=0.1; */
/* 	double step_radial=500; */
/* 	double cutValue=startingCut; */
/* 	double accept=0; */
/* 	double mistagged=0; */
/* 	double rejection=0; */
/* 	double sliceWidth=100; */
/* 	double threshold=0.9; */
/* 	TAxis* bi_x=compareBi210->GetXaxis(); */
/* 	TAxis* bi_y=compareBi210->GetYaxis(); */
/* 	TAxis* po_x=comparePo210->GetXaxis(); */
/* 	TAxis* po_y=comparePo210->GetYaxis(); */

/* 	TH1D* complete_bi=compareBi210->ProjectionY("complete_bi",bi_x->FindBin(0.),bi_x->FindBin(8000.)); */
/* 	TH1D* complete_po=comparePo210->ProjectionY("complete_po",po_x->FindBin(0.),po_x->FindBin(8000.)); */
/* 	double bi_numEV_complete= compareBi210->GetEntries(); */
/* 	double po_numEV_complete= comparePo210->GetEntries(); */
/* 	cout<<"Number of entries in complete bi = "<< bi_numEV_complete<<endl; */
/* 	cout<<"Number of entries in complete po = "<< po_numEV_complete<<endl; */

/* 	for(double energy =0; energy<6000;energy+=step_radial){ */
/* 		cout<<"radius = "<< energy<<endl; */

/* 		//TH1D* slice_bi=compareBi210->ProjectionY("slice_bi",bi_x->FindBin(0.5),bi_x->FindBin(1.0)); */
/* 		TH1D* slice_bi=compareBi210->ProjectionY(SSTR(energy).c_str(),bi_x->FindBin(energy),bi_x->FindBin(energy+sliceWidth)); */
/* 		TH1D* slice_po=comparePo210->ProjectionY(("Po"+SSTR(energy)).c_str(),po_x->FindBin(energy),po_x->FindBin(energy+sliceWidth)); */

/* 		double bi_numEV=slice_bi->GetEntries(); */
/* 		double po_numEV=slice_po->GetEntries(); */
/* 		cout<<"Number of entries in bi slice = "<< bi_numEV<<endl; */
/* 		cout<<"Number of entries in po slice = "<< po_numEV<<endl; */
/* 		TAxis* bi_slice_x=slice_bi->GetXaxis(); */
/* 		TAxis* po_slice_x=slice_po->GetXaxis(); */

/* 		while(accept<threshold){ */

/* 			double bi_remain=slice_bi->Integral(bi_slice_x->FindBin(startingCut),bi_slice_x->FindBin(cutValue+=step)); */
/* 			//double po_remain=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(startingCut+step)); */
/* 			//cout<<"Number of entries remaining in po slice = "<< po_remain<<endl; */
/* 			accept=bi_remain/bi_numEV; */	
/* 		} */
/* 		mistagged=slice_po->Integral(po_slice_x->FindBin(startingCut),po_slice_x->FindBin(cutValue)); */

/* 		//This rejection is a cheat to get around div0. */
/* 		if(po_numEV>1 && mistagged>1){ */
/* 		rejection= po_numEV/(mistagged+0.000001); */
/* 		}else{ */
/* 		rejection= 0; */
/* 		} */
/* 		cout<<mistagged<<" events were mistagged. Rejection =  "<<rejection<<"."<<endl; */
/* 		rejection_values.push_back(rejection); */
/* 		cutValues.push_back(cutValue); */
/* 		energyValues.push_back(energy+step_radial/2); */
/* 		if(false){ */
/* 		TCanvas* c1 = new TCanvas(); */
/* 		c1->cd(); */
/* 		//TLine* cutLine = new TLine( cutValue, slice_bi->GetYaxis()->GetXmin(), cutValue, slice_bi->GetYaxis()->GetXmax() ); */
/* 		TLine* cutLine = new TLine( cutValue, slice_bi->GetMinimum(), cutValue, slice_bi->GetMaximum() ); */
/* 		cutLine->SetLineWidth(2);cutLine->SetLineStyle(4); */
/* 		slice_bi->Draw(); */
/* 		slice_po->Draw("same"); */
/* 		cutLine->Draw("same"); */
/* 		} */
/* 		cutValue=startingCut; */
/* 		accept=0; */
/* 	} */

/* } */

int main(){
//====================================================================	
	gStyle->SetOptStat(0);

	TH3D* ele_hist= new TH3D("ele_hist","berkeleyAlphaBeta",18,0,3.6,60,0,6000,40,-200,200);
	ele_hist->SetLineColor(kBlue);ele_hist->SetLineWidth(3);
	ele_hist->SetMarkerColor(kBlue);
	ele_hist->SetFillColor(kBlue);
	TH3D* alpha_hist= new TH3D("alpha_hist","berkeleyAlphaBeta",18,0,3.6,60,0,6000,40,-200,200);
	alpha_hist->SetLineColor(kRed);alpha_hist->SetLineWidth(3);
	alpha_hist->SetMarkerColor(kRed);
	alpha_hist->SetFillColor(kRed);

	
	TLegend* t1 = new TLegend( 0.11, 0.11, 0.21, 0.21 );
	
	vector<string>  biFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
	vector<string> poFileList= glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

	for( int i=0; i<biFileList.size(); i++ ){
		TFile * file= TFile::Open(biFileList[i].c_str());	
		partflag="Bi";
		FillHist( file, ele_hist);
	  }

	for( int i=0; i<poFileList.size(); i++ ){
		TFile * file= TFile::Open(poFileList[i].c_str());	
		partflag="Po";
		FillHist( file, alpha_hist);
	  }

	/* vector<double> bi_cuts, energyValues,Rejection_energy; */
	/* vector<double> radial_cuts, RadiusValues,Rejection_radius; */

	/* findCutsEnergy(compareBi210,comparePo210,energyValues,bi_cuts,Rejection_energy); */
	/* findCutsRadial(compareEle_radial,compareAlpha_radial,RadiusValues,radial_cuts,Rejection_radius); */
/* // ----Finding Rejection Powers--------- */
	

/* //------Drawing and saving--------------- */


	/* TGraph* cutGraph = new TGraph(bi_cuts.size(),&energyValues[0],&bi_cuts[0]); */
	/* TF1 *f = new TF1("f", "[3]*x*x*x +[2]*x*x +[1]*x +[0]"); */
	/* cutGraph->Fit(f); */

	/* TGraph* cutGraph_radial = new TGraph(radial_cuts.size(),&RadiusValues[0],&radial_cuts[0]); */
	/* TF1 *f_Rad = new TF1("f_Rad", "[1]*x +[0]"); */
	/* cutGraph_radial->Fit(f_Rad); */



	TCanvas * c2 = new TCanvas();
	cout<<"# of entries in ele_hist= "<< ele_hist->GetEntries()<<endl;
	cout<<"# of entries in alpha_hist = "<< alpha_hist->GetEntries()<<endl;
	ele_hist->SetTitle("Classifiers values across truth quenched energy and reconstruted radius");
	ele_hist->GetXaxis()->SetTitle("mcEdepQuenched (MeV)");
	ele_hist->GetYaxis()->SetTitle("posr(mm)");
	ele_hist->GetZaxis()->SetTitle("BerkeleyAlphaBeta");
	ele_hist->Draw();
	alpha_hist->Draw("same");

	t1->AddEntry(ele_hist, "beta","f");
	t1->Draw();

	//c2->Print("compareBerkeleyAlphaBetaVsMCEnergy_withCuts.png");
	TFile fileout("3D_Plot.root","RECREATE");
	fileout.cd();
	alpha_hist->Write();
	ele_hist->Write();
	fileout.Close();

	return 0;

}

