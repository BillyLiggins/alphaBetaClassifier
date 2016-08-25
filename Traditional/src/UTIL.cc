#include "UTIL.h"

//OXSX libaries

// #include <BinnedPdf.h>        
// #include <Rand.h>
// #include <ROOTNtuple.h>        
// #include <BinnedNLLH.h>        
// #include <GridSearch.h>
// #include <MetropolisHastings.h> 
// #include <Minuit.h> 
// #include <Histogram.h>
// #include <PdfConverter.h>        
// #include <CompositePdf.h>        
// #include <Convolution.h>        
// #include <Gaussian.h>        
// #include <DataSetGenerator.h>        
// #include <OXSXDataSet.h>        
// #include <BoolCut.h>        
// #include <BoxCut.h>        

//ROOT libaries

#include <TString.h>
#include <TH1D.h>        
#include <TH2D.h>        
#include <TFrame.h> 
#include <TGraph.h> 
#include <TGraphErrors.h> 
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
#include <THStack.h> 
#include <TPad.h> 
#include <TPaveStats.h> 
#include <TAttFill.h>

//Standard libaries

#include <iomanip>
#include <sstream>
#include <vector>        
#include <string>        
#include <math.h>	

std::vector<std::string> UTIL::glob( const std::string& path, const std::string& start )
{
				std::vector<std::string> result;
				TSystemDirectory dir(path.c_str(), path.c_str());
				TList *files = dir.GetListOfFiles();
				if (files) {
								TSystemFile *file;
								TString fname;
								TIter next(files);
								while ((file=(TSystemFile*)next())) {
												fname = file->GetName();
												if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith( start ) ) {
												// if (!file->IsDirectory()  && fname.BeginsWith( start ) ) {
																std::stringstream fullPath; 
																fullPath << path << "/" << fname.Data();
																result.push_back(fullPath.str());
												}
								}
				}
				return result;
}


TH1D* UTIL::diffHist(TH1D * h1,TH1D * h2){
				/*This Finds the fractional difference between two histograms, (see code for details). 
				 */


				double minBin=h1->GetXaxis()->GetXmin();
				double maxBin=h1->GetXaxis()->GetXmax();
				double sliceWidth=h1->GetXaxis()->GetBinWidth(1);
				double numOfBins=h1->GetNbinsX();
				std::cout<<"minBin = "<<minBin<<std::endl;
				std::cout<<"maxBin = "<<maxBin<<std::endl;
				std::cout<<"sliceWidth = "<<sliceWidth<<std::endl;
				std::cout<<"number of bins = "<<numOfBins<<std::endl;

				TH1D* rhist = new TH1D("rhist","",numOfBins,minBin,maxBin);
				for(double i=0;i<numOfBins;i++){
								double h1cont=h1->GetBinContent(i);
								double h2cont=h2->GetBinContent(i);
								double weight;
								if (h1cont!=0 && h1cont-h2cont!=0) {
												weight= (h1cont-h2cont)/h1cont;
								} else {
												weight= 0;
								}
								rhist->SetBinContent(i,weight);
								// std::cout << "weight = "<<weight << std::endl;

				}
				// TCanvas * c1 =new TCanvas();
				// rhist->Draw();
				// c1->Print("temp.png");

				return rhist;


}

double UTIL::find2n2bRate(double percentage,double fidRad)
{
				//This function returns the 2n2b rate given a % loading and fid radius.
				double Vratio= pow(fidRad,3)/pow(6000,3);	
				std::cout << "Volume ratio = "<<Vratio << std::endl;
				double massOfLABPPO=(7.8e5)*Vratio;
				std::cout << "Mass of LABPPO = "<<massOfLABPPO << std::endl;
				double rawTeMass= percentage*massOfLABPPO/100;
				std::cout << "Raw mass of Te = "<<rawTeMass << std::endl;
				double _130TeMass= 0.3408*rawTeMass;
				std::cout << "130Te mass = "<<_130TeMass << std::endl;
				double _130TeHalflife= 7e20;
				std::cout << "130Te halflife = "<<_130TeHalflife << std::endl;
				double decayrate=log(2)/_130TeHalflife;
				std::cout << "decay rate = "<<decayrate << std::endl;
				double molarMassTe=129.906;
				std::cout << "130Te molar mass = "<<molarMassTe << std::endl;
				double molesOfTe=molarMassTe*_130TeMass;
				std::cout << "Number of moles of 130Te = "<<molesOfTe << std::endl;
				double numberOfParticles=molesOfTe*molarMassTe;
				std::cout << "num of particles = "<<numberOfParticles<< std::endl;
				double numOfEvents=numberOfParticles*decayrate*6.022e23;
				std::cout << "Number of events  = "<<numOfEvents  << std::endl;
				std::cout << "check this function again."  << std::endl;

				return numOfEvents;
}
//After chatting to Jeanne she said you may need to apply the decaying forumla.


double UTIL::scaleForTime(double yearRate,double runtime){
				//Should be given the runtime in days and yearly rates (events per year).
				return yearRate*runtime/365;

}

