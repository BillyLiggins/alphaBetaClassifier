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

// void findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values,vector<double>& Rejection_errors,vector<double>& energy_errors);
void findCutsEnergy();

private:
				TH2D * BabVsEnergy;
				// ofstream& outputfile;
				std::string PID;
				vector<double> energyValues,cutValues, Rejection_values, Rejection_errors, energy_errors;

};



#endif
