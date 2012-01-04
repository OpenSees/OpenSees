#include <FOSeriesSimulation.h>
FOSeriesSimulation::FOSeriesSimulation(int passedMaxSim,
									   int passedInterval,
									   double passedEps,
									   bool passedtowside,
									   int passedanalysis,
									   bool passedprint)
{
	print=passedprint;
	if(print){
		output.open("FOSeries.txt", ios::out);
		output << "\n";
		output << "FOSeries::FOSeries\n";
		output << "\n";
		output.flush();
	}
	Nrv=0; // Number of Random Variables;
	Ncomp=0; // Number of Componenet;;
	BetaVec=0; // Vector of beta ( size = nComp );
	AlphaVec=0;  // alpha Vectors ( ncomp Vectors of size nRV)
	uDesVec=0;
	MaxSim=passedMaxSim; // Maximum Number of Simulation
	CheckInterVal=passedInterval; // checking cov with this interval
	Eps=passedEps; // convergence criteria for cov;
	twoSide=passedtowside;
	analysisType=passedanalysis;
    theRandomNumberGenerator = new GeneralRandGenerator();
	if (theRandomNumberGenerator==0) {
		opserr << "FOSeriesSimulation::FOSeriesSimulation() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}

	ytrial=0;
	PfVec=0;
	Weight=0;
	WLim=0;
	ajVec=0;
} 
FOSeriesSimulation::FOSeriesSimulation
					(int passedNrv,
					 int passedNcomp,
					 int passedMaxSim,
					 int passedCheckInterval,
					 double passedEps,
					 Vector* passedBetaVec,
					 Vector** passedAlphaVec,
					 Vector** passeduDesVec,
					 int passedType,
				     bool passedtwoSide,
					 bool passedprint)
{
	print=passedprint;
	if(print){
		output.open("FOSeries.txt", ios::out);
		output << "\n";
		output << "FOSeries::FOSeries\n";
		output << "\n";
		output.flush();
	}
	Nrv=passedNrv; 
	Ncomp=passedNcomp;
	BetaVec=passedBetaVec;
	AlphaVec=passedAlphaVec; 
	uDesVec=passeduDesVec; 
	MaxSim=passedMaxSim; 
	CheckInterVal=passedCheckInterval; 
	Eps=passedEps; 
	twoSide=passedtwoSide;
	analysisType=passedType;
    theRandomNumberGenerator = new GeneralRandGenerator();
	if (theRandomNumberGenerator==0) {
		opserr << "FOSeriesSimulation::FOSeriesSimulation() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}

	ytrial=0;
	PfVec=0;
	Weight=0;
	WLim=0;
	ajVec=0;
} 
FOSeriesSimulation::~FOSeriesSimulation()
{
	if(ytrial!=0){delete ytrial;ytrial=0;}
	if(PfVec != 0) { delete PfVec; PfVec=0; } 
	if(Weight != 0) { delete Weight; Weight=0; } 
	if(WLim != 0) { delete WLim; WLim=0; } 
	if(ajVec != 0) { delete ajVec; ajVec=0; } 
} 
void FOSeriesSimulation::setNrv(int passedNrv)
{
	Nrv=passedNrv;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setNrv\n";
		output<<"Nrv is set to "<<Nrv<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setNcomp(int passedNcomp)
{
	Ncomp=passedNcomp;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setNcomp\n";
		output<<"Ncomp is set to "<<Ncomp<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setMaxSim(int passedMaxSim)
{
	MaxSim=passedMaxSim;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setMaxSim\n";
		output<<"Maxsim is set to "<<MaxSim<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setCheckInterVal(int passedValue)
{
	CheckInterVal=passedValue;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setCheckInterVal\n";
		output<<"CheckInterVal is set to "<<CheckInterVal<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setEps(double passedEps)
{
	Eps=passedEps;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setEps\n";
		output<<"Eps is set to "<<Eps<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setBetaVec(Vector* passedVec)
{
	BetaVec=passedVec;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setBetaVec\n";
		output<<"Size of Vector"<<BetaVec->Size();
		for(int i=0;i<BetaVec->Size();i++){
			output<<setw(5)<<i;
			output<<setw(15)<<setprecision(5)<<(*BetaVec)(i);
			output<<"\n";
		}
		output.flush();
	}
}
void FOSeriesSimulation::setAlphaVec(Vector** passedVec)
{
	AlphaVec=passedVec;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setAlphaVec\n";
		output<<"Size of Vector"<<AlphaVec[0]->Size()<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setuDesVec(Vector** passedVec)
{
	uDesVec=passedVec;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setuDesVec\n";
		output<<"Size of Vector"<<uDesVec[0]->Size()<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setAnalysisType(int passedType)
{
	analysisType=passedType;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setAnalysisType\n";
		output<<"AnalysisType"<<analysisType<<"\n";
		output.flush();
	}
}
void FOSeriesSimulation::setTwoSide(bool passedTwoSide)
{
	twoSide=passedTwoSide;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::setTwoSide\n";
		output<<"TwoSide"<<twoSide<<"\n";
		output.flush();
	}
}
int FOSeriesSimulation::selectComp()
{
	int i;
	double Uniform=theRandomNumberGenerator->generate_singleUniformNumber();
	int ifunc=Ncomp-1;
	for(i=0;i<Ncomp;i++){
	  if(Uniform<=(*WLim)(i)) break;
	}
	ifunc=i;
	if(print){
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		output<<"FOSeriesSimulation::selectComp\n";
		output<<"uniform RV"<<Uniform<<"\n";
		output<<"selevted component"<<ifunc<<"\n";
		output.flush();
	}

	return ifunc;
}
int FOSeriesSimulation::analyze(void)
{
// check input
	this->inputcheck();
	if(ytrial!=0){delete ytrial;ytrial=0;}
	ytrial= new Vector(Nrv);
	if (ytrial==0) {
		opserr << "FOSeriesSimulation::analyze1() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}

	if(analysisType==0){
		this->analyze0();
	}else if(analysisType==1){
		this->analyze1();
	}else if(analysisType==2){
		this->analyze2();
	}
	if(Weight!=0){ delete Weight; Weight=0;}
	if(PfVec!=0){ delete PfVec; PfVec=0;}
	if(WLim!=0){ delete WLim; WLim=0;}
	if(ytrial!=0){ delete ytrial; ytrial=0;}
	if(ajVec!=0){ delete ajVec; ajVec=0;}

	return 0;
}
void FOSeriesSimulation::analyze0(void)
{
	int i,j,k,numFail;
	double sum,sp,rnm,dev,xpp;
	int iconv=0;
	double coef=1.0;
	double cum=0.0;
	double cum2=0.0;
	int n1=CheckInterVal;

	static NormalRV aStdNormRV(1,0.0,1.0);

	opserr.setFloatField(SCIENTIFIC);
	opserr<<"\n";
	opserr<<"===== crude FO series simulation =====\n";
	opserr<<"\n";
	opserr<<"Number of Components.................."<<Ncomp<<"\n";
	opserr<<"Number of Random Variables............"<<Nrv<<"\n";
	opserr<<"\n";
	opserr<<"Maximum number of simulations........."<<MaxSim<<"\n";
	opserr<<"Cecking Interval......................"<<CheckInterVal<<"\n";
	opserr<<"Convergence criteria.................."<<Eps<<"\n";
	opserr<<"\n";
	opserr<<"TwoSide..............................."<<twoSide<<"\n";
	opserr<<"\n";
	opserr<<"      nsim";
	opserr<<" prob. estimate";
	opserr<<"      rel index";
	opserr<<"          c.o.v";
	opserr<<"\n";


	for(i=0;i<MaxSim;i++){
		for(j=0;j<Nrv;j++) (*ytrial)(j)=theRandomNumberGenerator->
							generate_singleStdNormalNumber();
		numFail=0;
		for(j=0;j<Ncomp;j++){
			sum=0.0;
			for(k=0;k<Nrv;k++)sum-=(*AlphaVec[j])(k)*(*ytrial)(k);
			if(!twoSide){
				sum+=(*BetaVec)(j);
				if(sum<=0.0) numFail++;
			}else{
				if(fabs(sum)>=(*BetaVec)(j)) numFail++;
			}
		}
		if(numFail!=0)xpp=1.0;
		else xpp=0.0;
		cum+=xpp;
		cum2+=xpp*xpp;

		if(i>n1) {
		  rnm=(double)i;
		  pfres=cum/rnm;
		  sp=cum2-cum*cum/rnm;
		  sp=sp/(rnm*(rnm-1));
		  if(sp>0.0) {
			  dev=sqrt(sp);
              cvar=dev/pfres;
			  betares = -aStdNormRV.getInverseCDFvalue(pfres);
              if(cvar<=Eps) iconv=1;
			  opserr.width(10);
			  opserr<<i;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<pfres;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<betares;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<cvar;
			  opserr<<"\n";
		  }else{
			  opserr.width(10);
			  opserr<<i;
			  opserr<<"\n";
		  }
		  n1+=CheckInterVal;
		}
		if(iconv==1) break;
	}
	if( iconv==1 ) numSimulation=i;
	else numSimulation=-i;
}
void FOSeriesSimulation::analyze1(void)
{
	opserr.setFloatField(SCIENTIFIC);
	opserr<<"\n";
	opserr<<"===== importance sampling FO series simulation =====\n";
	opserr<<"\n";
	opserr<<"Number of Components.................."<<Ncomp<<"\n";
	opserr<<"Number of Random Variables............"<<Nrv<<"\n";
	opserr<<"\n";
	opserr<<"Maximum number of simulations........."<<MaxSim<<"\n";
	opserr<<"Cecking Interval......................"<<CheckInterVal<<"\n";
	opserr<<"Convergence criteria.................."<<Eps<<"\n";
	opserr<<"\n";
	opserr<<"TwoSide..............................."<<twoSide<<"\n";
	opserr<<"\n";
	opserr<<"      nsim";
	opserr<<" prob. estimate";
	opserr<<"      rel index";
	opserr<<"          c.o.v";
	opserr<<"\n";

	this->makeWeights();

	int iconv=0;
	int i,j,iComp,numFail;
	double coef,sum,rnm,sp,dev,xpp,sample;
	double cum=0.0;
	double cum2=0.0;
	int n1=CheckInterVal;

	static NormalRV aStdNormRV(1, 0.0, 1.0);

	coef=this->hFunc(true);
	sample=0.0;
	for(i=0;i<MaxSim;i++){
	    iComp=this->selectComp();
		for(j=0;j<Nrv;j++) {
			(*ytrial)(j)=theRandomNumberGenerator->generate_singleStdNormalNumber();
			(*ytrial)(j)+=(*uDesVec[iComp])(j);
		}
		coef=this->hFunc(false);

		numFail=0;
		for(j=0; j<Ncomp; j++){
			sum=0.0;
			for(int k=0;k<Nrv;k++) sum-=(*AlphaVec[j])(k)*(*ytrial)(k);
			if(!twoSide){
				sum+=(*BetaVec)(j);
				if(sum<=0.0) numFail++;
			}else{
				if(fabs(sum)>=(*BetaVec)(j)) numFail++;
			}
		}
		if(numFail!=0){
			xpp=1.0*coef;
			sample+=1.0;
		}else{
			xpp=0.0;
		}
		cum+=xpp;
		cum2+=xpp*xpp;

		if(i>n1) {
		  rnm=(double)i;
		  pfres=cum/rnm;
		  sp=cum2-cum*cum/rnm;
		  sp=sp/(rnm*(rnm-1));
		  if(sp>0.0) {
			  dev=sqrt(sp);
              cvar=dev/pfres;
			  betares = -aStdNormRV.getInverseCDFvalue(pfres);
              if(cvar<=Eps) iconv=1;
			  opserr.width(10);
			  opserr<<i;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<pfres;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<betares;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<cvar;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<sample/rnm;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<coef;
			  opserr<<"\n";
		  }else{
			  opserr.width(10);
			  opserr<<i;
			  opserr<<"\n";
		  }
		  n1+=CheckInterVal;
		}
		if(iconv==1) break;
	}
	if( iconv==1 ) numSimulation=i;
	else numSimulation=-i;

}
void FOSeriesSimulation::analyze2(void)
{
	char outputString[100];

	opserr.setFloatField(SCIENTIFIC);
	opserr<<"\n";
	opserr<<"===== importance sampling FO series simulation =====\n";
	opserr<<"\n";
	opserr<<"Number of Components.................."<<Ncomp<<"\n";
	opserr<<"Number of Random Variables............"<<Nrv<<"\n";
	opserr<<"\n";
	opserr<<"Maximum number of simulations........."<<MaxSim<<"\n";
	opserr<<"Cecking Interval......................"<<CheckInterVal<<"\n";
	opserr<<"Convergence criteria.................."<<Eps<<"\n";
	opserr<<"\n";
	opserr<<"TwoSide..............................."<<twoSide<<"\n";
	opserr<<"\n";
	opserr<<"      nsim";
	opserr<<" prob. estimate";
	opserr<<"      rel index";
	opserr<<"          c.o.v";
	opserr<<"\n";

	this->makeWeights();
	if(twoSide) pfsum*=2.0;

	int iconv=0;
	int i,j,k,iComp,numFail;
	double sum,rnm,sp,dev,xpp,sample;
	double cum=0.0;
	double cum2=0.0;
	int n1=CheckInterVal;
	double uniform,ptemp,qtemp,ytemp,positive;

	static NormalRV aStdNormRV(1, 0.0, 1.0);

	sample=0.0;
	for(i=0;i<MaxSim;i++){
	    iComp=this->selectComp();
		positive=1.0;
		if(twoSide){
			uniform=theRandomNumberGenerator->generate_singleUniformNumber();
			if(uniform>=0.5) positive=1.0;
			else positive=-1.0;
		}
		sum=0.0;
		for(j=0;j<Nrv;j++){
			(*ytrial)(j)=theRandomNumberGenerator->generate_singleStdNormalNumber();
			sum+=positive*(*AlphaVec[iComp])(j)*(*ytrial)(j);
		}
		for(j=0; j<Nrv; j++)(*ytrial)(j)-=sum*positive*(*AlphaVec[iComp])(j);
		uniform=theRandomNumberGenerator->generate_singleUniformNumber();
		ptemp=1.0+(*PfVec)(iComp)*(uniform-1.0);
        qtemp=1.0-ptemp;
		ytemp=aStdNormRV.getInverseCDFvalue(ptemp);
		for(j=0;j<Nrv;j++) (*ytrial)(j)+=ytemp*(*AlphaVec[iComp])(j);
		numFail=0;
		for(int j=0;j<Ncomp;j++){
			sum=0.0;
			for(k=0;k<Nrv;k++)sum-=(*AlphaVec[j])(k)*(*ytrial)(k);
			if(!twoSide){
				sum+=(*BetaVec)(j);
				if(sum<=0.0) numFail++;
			}else{
				if(fabs(sum)>=(*BetaVec)(j)) numFail++;
			}
		}
		if(numFail!=0) sample++;
		xpp=pfsum/(double)numFail;

		cum+=xpp;
		cum2+=xpp*xpp;

		if(i>n1) {
		  rnm=(double)i;
		  pfres=cum/rnm;
		  sp=cum2-cum*cum/rnm;
  		  if(fabs(sp)<=1.0e-3) sp=0.0;
		  sp=sp/(rnm*(rnm-1));
		  if(sp>=0.0) {
	 		  dev=sqrt(sp);
              cvar=dev/pfres;
			  if(pfres>=1.0) pfres=0.99999999999999;
			  betares = -aStdNormRV.getInverseCDFvalue(pfres);
              if(cvar<=Eps) iconv=1;

			  sprintf(outputString,"%10d%15.5e%15.5e%15.5e%15.5e" 
				  ,i,pfres,betares,cvar,sample/rnm);
			  opserr << outputString << endln;

/*			  opserr.width(10);
			  opserr<<i;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<pfres;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<betares;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<cvar;
			  opserr.width(15);
			  opserr.setPrecision(5);
			  opserr<<sample/rnm;
			  opserr<<"\n";*/
		  }else{
			  opserr.width(10);
			  opserr<<i;
			  opserr<<"\n";
		  }
		  n1+=CheckInterVal;
		}
		if(iconv==1) break;
	}
	if( iconv==1 ) numSimulation=i;
	else numSimulation=-i;
}
int FOSeriesSimulation::makeWeights()
{
	if(PfVec != 0) { delete PfVec; PfVec=0; } 
	PfVec = new Vector(Ncomp);
	if (PfVec==0) {
		opserr << "FOSeriesSimulation::makeWeights() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}

	static NormalRV aStdNormRV(1, 0.0, 1.0);

	double beta,pf;
	pfsum=0.0;
	for(int i=0; i<Ncomp; i++){
		beta=(*BetaVec)(i);
		pf=1.0 - aStdNormRV.getCDFvalue(beta);
		(*PfVec)(i)=pf;
		pfsum+=pf;
	}

	if(Weight != 0) { delete Weight; Weight=0; } 
	if(WLim != 0) { delete WLim; WLim=0; } 
	Weight = new Vector(Ncomp);
	WLim = new Vector(Ncomp);
	if (Weight==0) {
		opserr << "FOSeriesSimulation::makeWeights() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}
	if (WLim==0) {
		opserr << "FOSeriesSimulation::makeWeights() - " << endln
			<< " out of memory while instantiating internal objects." << endln;
		exit(-1);
	}

	(*Weight)(0)=(*PfVec)(0)/pfsum;
    (*WLim)(0)=(*Weight)(0);	 
	for(int i=1;i<Ncomp; i++){
	    (*Weight)(i)=(*PfVec)(i)/pfsum;
        (*WLim)(i)=(*WLim)(i-1)+(*Weight)(i);
	}

	if(print){
		output<<"FOSeriesSimulation::makeWeight\n";
 		output<<"num of Comp"<<Ncomp<<"\n";
		output<<"      comp";
		output<<"           beta";
		output<<"             pf";
		output<<"         weight";
		output<<"           wlim\n";
		output<<
		output.setf(ios::right);
		output.setf(ios::scientific, ios::floatfield);
		for(int i=0;i<Ncomp;i++){
			output<<setw(10)<<i;
			output<<setw(15)<<setprecision(5)<<(*BetaVec)(i);
			output<<setw(15)<<setprecision(5)<<(*PfVec)(i);
			output<<setw(15)<<setprecision(5)<<(*Weight)(i);
			output<<setw(15)<<setprecision(5)<<(*WLim)(i);
			output<<"\n";
		}
		output.flush();
	}

	return 0;
}

double FOSeriesSimulation::hFunc(bool preComp)
{
	int i,j;
	double sum,udum,sum0,coef;

	if(ajVec != 0) { delete ajVec; ajVec=0; } 
	ajVec = new Vector(Ncomp);

	if(preComp){
		for(i=0;i<Ncomp;i++){
			sum=0.0;
			for(j=0;i<Nrv;j++){
				udum=(*uDesVec[i])(j);
				sum+=udum*udum;
			}
			(*ajVec)(i)=-sum;
		}
		coef=0.0;
	}else{
//		sum0=0.0;
//		for(j=0;j<Nrv;j++) sum0+=(*ytrial)(j)*(*ytrial)(j)
		sum0=0.0;
		coef=0.0;
		for(i=0;i<Ncomp;i++){
			sum=0.0;
			for(j=0;j<Nrv;j++) sum+=(*uDesVec[i])(j)*(*ytrial)(j);
            coef+=(*Weight)(i)*exp((*ajVec)(i)+sum+sum0);
		}
		coef=1.0/coef;
	}
//	    COEF=COEF*(SIGMA**FLOAT(NUM))
	return coef;
}
void FOSeriesSimulation::inputcheck(void)
{
	if(Nrv==0){
		opserr<< "Numbero of Random Variablles : Nrv is zero\n";
		opserr<< "check the data\n";
		exit(-1);
	}
	if(Ncomp==0){
		opserr<< "Numbero of Components : Ncomp is zero\n";
		opserr<< "check the data\n";
		exit(-1);
	}
	if(BetaVec==0){
		opserr<< "Beta vector is not set yet\n";
		opserr<< "check the data\n";
		exit(-1);
		if(BetaVec->Size()<Ncomp){ 
			opserr<< "Size of Beta is not consistent\n";
			opserr<< "check the data\n";
			opserr<< "Size of Beta"<<BetaVec->Size()<<"\n";
			opserr<< "Ncomp"<<Ncomp<<"\n";
			exit(-1);
		}
	}

	if(AlphaVec==0){
		opserr<< "Alpha vectors are not set yet\n";
		opserr<< "check the data\n";
		exit(-1);
		if(AlphaVec[0]->Size()<Nrv){ 
			opserr<< "Size of AlphaVec is not consistent\n";
			opserr<< "check the data\n";
			opserr<< "Size of Alpha"<<AlphaVec[0]->Size()<<"\n";
			opserr<< "Nrv"<<Nrv<<"\n";
			exit(-1);
		}
	}

	if(uDesVec==0){
		opserr<< "uDes vectors are not set yet\n";
		opserr<< "check the data\n";
		exit(-1);
		if(uDesVec[0]->Size()<Nrv){ 
			opserr<< "Size of uDesVec is not consistent\n";
			opserr<< "check the data\n";
			opserr<< "Size of uDes"<<uDesVec[0]->Size()<<"\n";
			opserr<< "Nrv"<<Nrv<<"\n";
			exit(-1);
		}
	}

	if(MaxSim==0){
		opserr<< "MaxSim is not set yet\n";
		opserr<< "check the data\n";
		exit(-1);
	}
	
	if(CheckInterVal==0){
		opserr<< "CheckInterVal not set yet\n";
		opserr<< "check the data\n";
		exit(-1);
	}
	if(Eps==0.0){
		opserr<< "Eps is not set yet\n";
		opserr<< "check the data\n";
		exit(-1);
	}
	if(analysisType<0){
		opserr<< "analysisType is not set yet\n";
		opserr<< "check the data\n";
		exit(-1);
	}
}
