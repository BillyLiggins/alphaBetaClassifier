#include "UTIL.h"
#include "Cutter.h"

//ROOT libaries

#include <TString.h>
#include <TH1D.h>        
#include <TH2D.h>        
#include <TStyle.h> 
#include <TPaveStats.h> 
#include <TAttFill.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TMultiGraph.h>
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
#include <TPad.h>

//Standard libaries

#include <iomanip>
#include <sstream>
#include <vector>        
#include <string>        
#include <math.h>	
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>


int main(){

				ofstream File_bi;
				ofstream File_po;

				File_bi.open("Classifier_data_bi.dat");
				File_po.open("Classifier_data_po.dat");

				File_bi << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" <<  std::endl;
				File_po << "mcEdepQuenched,"<<"posr," << "BerekelyAlphaBeta" << std::endl;

				UTIL* util = new UTIL();

				// std::vector<std::string> betaFileList= util->glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output_electron/ntuple","electron");
				// std::vector<std::string> alphaFileList= util->glob("/data/snoplus/liggins/year1/fitting/fitting/alphaSims/output/ntuple","alpha");

				std::vector<std::string> betaFileList= util->glob("/home/billy/workspace/PhD/testData/beta/output/ntuple","electron");
				std::vector<std::string> alphaFileList= util->glob("/home/billy/workspace/PhD/testData/alpha/output/ntuple","alpha");

				Cutter* alpha = new Cutter("alpha");
				Cutter* beta = new Cutter("beta");

				alpha->SetHist();
				beta->SetHist();

				for( int i=0; i<betaFileList.size(); i++ ){
								TFile * file= TFile::Open(betaFileList[i].c_str());	
								std::cout<<"Loaded file "<<betaFileList[i]<<std::endl;
								beta->FillHist(file,File_bi);
				}

				beta->PrintHist();

				for( int i=0; i<alphaFileList.size(); i++ ){
								TFile * file= TFile::Open(alphaFileList[i].c_str());	
								std::cout<<"Loaded file "<<alphaFileList[i]<<std::endl;
								alpha->FillHist( file, File_po);
				}
				alpha->PrintHist();
				File_bi.close();
				File_po.close();

				return 0;
}
