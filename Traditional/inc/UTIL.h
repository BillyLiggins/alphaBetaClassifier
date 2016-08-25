#ifndef __UTIL__
#define __UTIL__

//ROOT libaries
#include <TH1D.h>

//Standard libaries
#include <string>
#include <vector>
#include <iostream>

class UTIL{
public:
			 	UTIL(){
				}

				std::vector<std::string> glob( const std::string& path, const std::string& start );

				TH1D* diffHist(TH1D * h1,TH1D * h2);

				double find2n2bRate(double percentage=0.5,double fidRad=6000);
				//After chatting to Jeanne she said you may need to apply the decaying forumla.

				double scaleForTime(double yearRate,double runtime);


private:

				double dummy;

};



#endif
