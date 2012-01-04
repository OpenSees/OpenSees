

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/FOSeriesSimulation.h,v $


#ifndef FOSeriesSimulation_h
#define FOSeriesSimulation_h

#include <Vector.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <RandomNumberGenerator.h>
#include <NormalRV.h>
#include <GeneralRandGenerator.h>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class FOSeriesSimulation
{
  public:
	FOSeriesSimulation(int passedMaxSim=2000,
					   int passedInterval=100,
					   double passedEps=0.01,
					   bool passedtowside=true,
					   int passedanalysis=2,
					   bool passedprint=false);
	FOSeriesSimulation(int passedNrv,
					   int passedNcomp,
					   int passeMaxSim,
					   int passedCheckInterval,
					   double passedEps,
					   Vector* passedBetaVec,
					   Vector** passedAlphaVec,
                       Vector** passeduDesVec,
					   int passedType,
				       bool passedtwoSide,
					   bool print=false);
	~FOSeriesSimulation();

	void setNrv(int);
	void setNcomp(int);
	void setMaxSim(int);
	void setCheckInterVal(int);
	void setEps(double);
	void setBetaVec(Vector*);
	void setAlphaVec(Vector**);
	void setuDesVec(Vector**);
	void setAnalysisType(int);
	void setTwoSide(bool);
	double getpfres(){return pfres;}
	double getcvar(){return cvar;}
	double getbetares(){return betares;}
	int getnumSimulation(){ return numSimulation;}

	int selectComp();
	int analyze(void);
	void analyze0(void);
	void analyze1(void);
	void analyze2(void);
	void inputcheck();
	int makeWeights();
	double hFunc(bool);


  protected: 
    
  private:
	int Nrv; // Number of Random Variables;
	int Ncomp; // Number of Componenet;;
	Vector* BetaVec; // Vector of beta ( size = nComp );
	Vector** AlphaVec;  // alpha Vectors ( ncomp Vectors of size nRV)
	Vector** uDesVec;  // design point Vectors ( ncomp Vectors of size nRV)
	int MaxSim; // Maximum Number of Simulation
	int CheckInterVal; // checking cov with this interval
	double Eps; // convergence criteria for cov;
	bool twoSide;
	int analysisType;
	RandomNumberGenerator* theRandomNumberGenerator;
	ofstream output;
	Vector* Weight;
	Vector* ytrial;
	Vector* ajVec;
	Vector* PfVec;
	Vector* WLim;
	double pfsum;
	double pfres;
	double cvar;
	double betares;
	int numSimulation;
	bool print;

};

#endif
