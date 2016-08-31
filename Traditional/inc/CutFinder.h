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
				
				void scan();


private:

				Cutter* acceptor;
				Cutter* rejector;
				double acceptor_entries;
				double rejector_entries;
				double threshold;
};



#endif
