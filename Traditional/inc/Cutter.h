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
				;}

			 	Cutter(std::string PID): PID(PID){
				;}

			 	~Cutter(){
				;}

void FillHist(TFile* file ,ofstream& outputfile);
void FillHist(TFile* file );

void ApplyCut(TFile * file );

void SetHist();

void SetGradient(double value){ gradient= value;}
void SetIntercept(double value){ intercept= value;}
				
double GetGradient(){ return gradient;}
double GetIntercept(){ return intercept;}

void SetHistLimits(double Ebins,double ELow,double EHigh,double BabBins, double BabLow, double BabHigh);

void PrintHist();

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

// void findCutsEnergy(TH2D *compareBi210,TH2D *comparePo210,vector<double>& energyValues, vector<double>& cutValues,vector<double>& Rejection_values,vector<double>& Rejection_errors,vector<double>& energy_errors);
void findCutsEnergy();


double GetNumberOfEntries(){return numberOfEntries;}
double GetRemainingAfterCut(){return remainingAfterCut;}

private:
				TH2D * BabVsEnergy;
				// ofstream& outputfile;
				double numberOfEntries;
				double gradient;
				double intercept;
				std::string PID;
				std::vector<double> energyValues,cutValues, Rejection_values, Rejection_errors, energy_errors, mistagged,mistagged_error;
				double remainingAfterCut;

};



#endif
