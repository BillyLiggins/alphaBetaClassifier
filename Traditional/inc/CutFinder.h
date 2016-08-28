#ifndef __CutFinder__
#define __CutFinder__

//Your libaries
#include "Cutter.h"
#include "Util.h"

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

				void findBoundary();
				void FindCutValue();

				


private:

				Cutter* acceptor;
				Cutter* rejector;
				double acceptor_entries;
				double rejector_entries;
};



#endif
