// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/OutCrossingResults.h,v $

#ifndef OutCrossingResults_h
#define OutCrossingResults_h

#include <Vector.h>
#include <TimePoints.h>
#include <ReliabilityDomain.h>
#include <Vector.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class OutCrossingResults
{
	public:
		OutCrossingResults(int passednumLsf, int passednumFrag, 
						   int passednRV, 
						   char* passedfileName, 
						   bool print=false);
		OutCrossingResults(OutCrossingResults& passedobject);

		~OutCrossingResults();

//		void clear(int,int, double, double, int, int);
		void clear(int, int, int);
		void setnumAna(int,int);
		void setLsf(int,int);
		void setnumAnaIncSens(int,int);
		void setnumLinSearch(int,int);
		void setxDesPoints(int,Vector&);
		void setuDesPoints(int,Vector&);
		void setDesAlpha(int,Vector&);
		void setHfunc(int,Vector&);
		void setthresholdVal(int,double);
		void setbeta(int,double);
		void setpf(int,double);
		void setnu(int,double);
		void settime(int,double);
		void setcheck1(int,double);
		void setcheck2(int,double);
		void setcheck1_init(int,double);
		void setcheck2_init(int,double);
		void setnumSteps(int,int);
		void setiresult(int,int);

//		int getIdLsf(void){return IdLsf;}
//		int getIdFragility(void){return IdFragility;}
//		double getthreshold(void){return threshold;}
		int getnumPoints(void){return numPoints;}
		int getnumRV(void){ return numRV; } 
		int getLsf(int);
		int getnumLsf(void);
		int getnumFragility(void);
		int getnumAna(int);
		int getnumAnaIncSens(int);
		int getnumLinSearch(int);
		int getnumSteps(int);
		int getiresult(int);
		Vector& getxDesPoints(int);
		Vector& getuDesPoints(int);
		Vector& getDesAlpha(int);
		Vector& getHfunc(int);
		double getthresholdVal(int);
		double getbeta(int);
		double getpf(int);
		double getnu(int);
		double getTime(int);
		double getcheck1(int);
		double getcheck2(int);
		double getcheck1_init(int);
		double getcheck2_init(int);
		void printResult(ofstream& outputFile);
		char* getfileName() { return fileName; }

//		int* getstepVector(){return numSteps;}
//		int* getnumAnaVector(){ return numAna;}
//		int* getnumAnaIncSensVector(){ return numAnaIncSens;}
//		int* getnumLinSearchVector(){ return numLinSearch;}
//		int* getiresultVector(){ return iresult;}

		void outtoFile(void);
		void readfromFile(void);
		void printResults(ofstream& outputFile);
        void printSinglePoint(ofstream& outputFile,int ipt);
        void printSingleTitle(ofstream& outputFile,int passedlsf);
	protected: 
    
	private:
//		void setIdLsf(int passedValue){IdLsf=passedValue;}
//		void setIdFragility(int passedValue){IdFragility=passedValue;}
//		void setthreshold(double passedValue){IdFragility=passedValue;}
//		void setfactorFragility(double passedValue){factorFragility=passedValue;}
//		void setnumPoints(int passedValue){numPoints=passedValue;}
//		void setnumRV(int passedValue){numRV=passedValue;}
		void allocate(int,int);
//		int IdLsf;  
//		int numLsf;  
//		int numthreshold;  
//		int IdFragility;
//		double factorFragility;
//		double threshold;
		int numPoints;
		int numLsf;
		int numFragility;
		int numRV;
		int* numAna;
		int* numAnaIncSens;
		int* numLinSearch;
		int* numSteps;
		int* iresult;
		int* Lsf;  
		Vector** xDesPoints;
		Vector** uDesPoints;
		Vector** alphaPoints;
		Vector** hfuncs;
		Vector* thresholdVal;
		Vector* betaPoints;
		Vector* pfPoints;
		Vector* nuPoints;
		Vector* resTime;
		Vector* check1;
		Vector* check2;
		Vector* check1_init;
		Vector* check2_init;

		bool print;
		ofstream output;
		char* fileName;

};

#endif
