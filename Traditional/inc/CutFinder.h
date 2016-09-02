#ifndef __CutFinder__
#define __CutFinder__

//Your libaries
#include "UTIL.h"
#include "Cutter.h"

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

class CutFinder{
public:
			 	CutFinder(){
				}

			 	CutFinder(Cutter* acceptor,Cutter* rejector): acceptor(acceptor), rejector(rejector){
				}

				~CutFinder(){
				}
				void FindBoundary();
				// void FindCutValue();

				double GetThreshold(){
								return threshold;
				}
				void SetThreshold(double value){
								threshold=value;
				}
				
				void Scan();
				void Holder();
				void Plotting();


				const std::vector<double> & GetRejectionValuesVector() const { return Rejection_values;}
				const std::vector<double> & GetRejectionErrorVector() const { return Rejection_errors;}
				const std::vector<double> & GetEnergyValuesVector() const { return energyValues;}
				const std::vector<double> & GetEnergyErrorVector() const { return energy_errors;}
				const std::vector<double> & GetMistaggedValueVector() const { return mistagged;}
				const std::vector<double> & GetMistaggedErrorVector() const { return mistagged_error;}

				const	std::vector<double> & GetCutValuesVector() const { return cutValues;}

				void  EnterRejectionValue(double entry) { Rejection_values.push_back(entry);}
				void  EnterRejectionError(double entry) { Rejection_errors.push_back(entry);}
				void  EnterEnergyValues(double entry) { energyValues.push_back(entry);}
				void  EnterEnergyError(double entry) { energy_errors.push_back(entry);}
				void  EnterMistaggedValue(double entry) { mistagged.push_back(entry);}
				void  EnterMistaggedError(double entry) { mistagged_error.push_back(entry);}
				void  EnterCutValues(double entry) { cutValues.push_back(entry);}

				void findCutsEnergy();


private:

				Cutter* acceptor;
				Cutter* rejector;
				double acceptor_entries;
				double rejector_entries;
				double threshold;
				std::vector<double> energyValues,cutValues, Rejection_values, Rejection_errors, energy_errors, mistagged,mistagged_error;
};



#endif
