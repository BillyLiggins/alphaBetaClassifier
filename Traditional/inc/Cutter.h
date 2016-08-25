#ifndef __Cutter__
#define __Cutter__

//ROOT libaries
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

//Standard libaries
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>

class Cutter{
public:
			 	Cutter(){
				}

			 	Cutter(std::string PID): PID(PID){
				}

void FillHist(TFile* file ,ofstream& outputfile);

void SetHist();

void PrintHist();



private:
				TH2D * BabVsEnergy;
				// ofstream& outputfile;
				std::string PID;

};



#endif
