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
								Cutter():  numberOfEntries(0),remainingAfterCut(0){
												;}

								Cutter(std::string PID): PID(PID),  numberOfEntries(0), remainingAfterCut(0) {
												;}

								~Cutter(){
												;}

								void FillHist(TFile* file ,ofstream& outputfile);
								void FillHist(TFile* file );
								void FillCutter(std::string folder,std::string fileStart);
								void FillCutter(std::string folder,std::string fileStart,ofstream&  File_bi);


								void ApplyCut(TFile * file );
								void ApplyBoundary(std::string folder,std::string fileStart);

								void SetHist();

								void SetGradient(double value){ gradient= value;}
								void SetIntercept(double value){ intercept= value;}
								void SetRemainingPercentage(double value){ percentageRemaining = value;}

								double GetGradient(){ return gradient;}
								double GetIntercept(){ return intercept;}
								double GetRemainingPercentage(){ return percentageRemaining;}

								double GetRemainingPercentageError(){ 
												percentageRemainingError= percentageRemaining*sqrt((remainingAfterCut+numberOfEntries)/(remainingAfterCut*numberOfEntries));
												return percentageRemainingError;}

								void SetHistLimits(double Ebins,double ELow,double EHigh,double BabBins, double BabLow, double BabHigh);

								void PrintHist();

								void SetRadialCut(double rad){ radialCut= rad;}

								TH2D* GetHist();


								std::string GetPID() const {return PID;}
								const std::vector<double> & GetRejectionValuesVector() const { return Rejection_values;}
								const std::vector<double> & GetRejectionErrorVector() const { return Rejection_errors;}
								const std::vector<double> & GetEnergyValuesVector() const { return energyValues;}
								const std::vector<double> & GetEnergyErrorVector() const { return energy_errors;}
								const std::vector<double> & GetMistaggedValueVector() const { return mistagged;}
								const std::vector<double> & GetMistaggedErrorVector() const { return mistagged_error;}

								std::vector<double> GetCutValuesVector() const { return cutValues;}

								void  EnterRejectionValue(double entry) { Rejection_values.push_back(entry);}
								void  EnterRejectionError(double entry) { Rejection_errors.push_back(entry);}
								void  EnterEnergyValues(double entry) { energyValues.push_back(entry);}
								void  EnterEnergyError(double entry) { energy_errors.push_back(entry);}
								void  EnterMistaggedValue(double entry) { mistagged.push_back(entry);}
								void  EnterMistaggedError(double entry) { mistagged_error.push_back(entry);}
								void  EnterCutValues(double entry) { cutValues.push_back(entry);}

								void findCutsEnergy(std::string folder, std::string fileStart);


								
								void FindNEntries(std::string folder, std::string fileStart);
								double GetNumberOfEntries(){return numberOfEntries;}
								double GetRemainingAfterCut(){return remainingAfterCut;}
								double SetRemainingAfterCut(double value){remainingAfterCut = value;}

								void FindRejection(std::string folder,std::string fileStart){

				private:
								TH2D * BabVsEnergy;
								TH2D * BabVsEnergyPassCut;
								double radialCut;
								// ofstream& outputfile;
								double numberOfEntries;
								double gradient;
								double intercept;
								std::string PID;
								std::vector<double> energyValues, cutValues, Rejection_values, Rejection_errors, energy_errors, mistagged, mistagged_error;
								std::vector<double> rejection, rejector_errors;
								double remainingAfterCut;
								double percentageRemaining;
								double percentageRemainingError;

};



#endif
