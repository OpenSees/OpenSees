// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/OutCrossingResults.cpp,v $

#include <OutCrossingResults.h>
OutCrossingResults::OutCrossingResults(int passednumLsf, int passednumFrag, 
									   int passednumRV, 
									   char* passedfileName, 
									   bool passedprint)
{
//	IdLsf=0;
//	IdFragility=0;
//	threshold=0.0;
	numLsf=0;
	numFragility=0;
	numPoints=0;
	numRV=0;
	Lsf=0;
	numAna=0;
	numAnaIncSens=0;
	numLinSearch=0;
	numSteps=0;
	iresult=0;
	xDesPoints=0;
	uDesPoints=0;
	hfuncs=0;
	alphaPoints=0;
	thresholdVal=0;
	betaPoints=0;
	pfPoints=0;
	nuPoints=0;
	resTime=0;
	check1=0;
	check2=0;
	check1_init=0;
	check2_init=0;

	print=passedprint;
	numLsf=passednumLsf;
	numFragility=passednumFrag;
	numPoints= numLsf*numFragility;
	numRV=passednumRV;

	if(print){
		output.open("OutCrossingResults.txt", ios::out);
		output << "\n";
		output << "OutCrossingResults::OutCrossingResults\n";
		output << "\n";
		output << "numPoints " << numPoints << "\n"; 
		output << "numRV " << numRV << "\n"; 
		output.flush();
	}
	allocate(numPoints,numRV);
	if(passedfileName==0){
		fileName="designpoints.bin";
	}else{
		fileName=passedfileName;
	}
}
OutCrossingResults::OutCrossingResults(OutCrossingResults& passedobject) 
{
	numLsf=passedobject.getnumLsf();
	numFragility=passedobject.getnumFragility();
	numPoints=passedobject.getnumPoints();
	numRV=passedobject.getnumRV();
	allocate(numPoints,numRV);

//	IdLsf=passedobject.getIdLsf();  
//	IdFragility=passedobject.getIdFragility();
	
//	threshold=passedobject.getthreshold();

	for(int i=0; i<numPoints; i++){
		Lsf[i]=passedobject.getLsf(i);
		numAna[i]=passedobject.getnumAna(i);
		numAnaIncSens[i]=passedobject.getnumAnaIncSens(i);
		numLinSearch[i]=passedobject.getnumLinSearch(i);
		numSteps[i]=passedobject.getnumSteps(i);
		iresult[i]=passedobject.getiresult(i);
		(*betaPoints)(i)=passedobject.getbeta(i);
		(*thresholdVal)(i)=passedobject.getthresholdVal(i);
		(*pfPoints)(i)=passedobject.getpf(i);
		(*nuPoints)(i)=passedobject.getnu(i);
		(*resTime)(i)=passedobject.getTime(i);
		(*check1)(i)=passedobject.getcheck1(i);
		(*check2)(i)=passedobject.getcheck2(i);
		(*check1_init)(i)=passedobject.getcheck1_init(i);
		(*check2_init)(i)=passedobject.getcheck2_init(i);
		Vector& tmpxDesPoint=passedobject.getxDesPoints(i);
		Vector& tmpuDesPoint=passedobject.getuDesPoints(i);
		Vector& tmpDesAlpha=passedobject.getDesAlpha(i);
		Vector& tmpHfunc=passedobject.getHfunc(i);
		for(int j=0; j<numRV; j++){
			(*xDesPoints[i])(j)=tmpxDesPoint(j);
			(*uDesPoints[i])(j)=tmpuDesPoint(j);
			(*alphaPoints[i])(j)=tmpDesAlpha(j);
			(*hfuncs[i])(j)=tmpHfunc(j);
		}
	}
	print=false;
	fileName=passedobject.getfileName();

}
OutCrossingResults::~OutCrossingResults()
{
	if(Lsf!=0){ delete [] Lsf; Lsf=0; }
	if(numAna!=0){ delete [] numAna; numAna=0; }
	if(numAnaIncSens!=0){ delete [] numAnaIncSens; numAnaIncSens=0;}
	if(numLinSearch!=0){delete [] numLinSearch; numLinSearch=0;}
	if(numSteps!=0){ delete [] numSteps; numSteps=0; }
	if(iresult!=0){ delete [] iresult; iresult=0; }
	if(xDesPoints!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete xDesPoints[i]; xDesPoints[i]=0; }
		delete [] xDesPoints;
		xDesPoints = 0;
	}
	if(uDesPoints!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete uDesPoints[i]; uDesPoints[i]=0; } 
		delete [] uDesPoints;
		uDesPoints = 0;
	}
	if(alphaPoints!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete alphaPoints[i]; alphaPoints[i]=0; }
		delete [] alphaPoints;
		alphaPoints=0;
	}
	if(hfuncs!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete hfuncs[i]; hfuncs[i]=0; } 
		delete [] hfuncs;
		hfuncs = 0;
	}
	if(thresholdVal!=0){ delete thresholdVal; thresholdVal=0 ;}
	if(betaPoints!=0){ delete betaPoints; betaPoints=0 ;}
	if(pfPoints!=0){ delete pfPoints;pfPoints=0;}
	if(nuPoints!=0){ delete nuPoints; nuPoints=0;}
	if(resTime!=0){delete resTime;resTime=0;}
	if(check1!=0){delete check1;check1=0;}
	if(check2!=0){delete check2;check2=0;}
	if(check1_init!=0){delete check1_init;check1_init=0;}
	if(check2_init!=0){delete check2_init;check2_init=0;}
}
//void OutCrossingResults::clear(int passedLsf,int passeFragility, 
//							   double passedthreshold, double passedfactor,
//							   int passednumPoints, int passednumRV)
void OutCrossingResults::clear(int passednumLsf, int passednumFrag, int passednumRV)
{
//	IdLsf=passedLsf;  
//	IdFragility=passeFragility;
//	threshold=passedthreshold;
//	factorFragility=passedfactor;
	int passednumPoints=passednumLsf*passednumFrag;
	if(numLsf!=passednumLsf||numFragility!=passednumFrag||
		passednumPoints!=numPoints||passednumRV!=numRV){
		opserr<< "!!!!!Warning!!!!!\n";
		opserr<< "Inconsistent numbers for theOutCrossingResults \n";
		opserr<< "passednumLsf "<<passednumLsf;
		opserr<< "numLsf "<<numLsf<<"\n";
		opserr<< "passednumFrag "<<passednumFrag;
		opserr<< "numFragility "<<numFragility<<"\n";
		opserr<< "passednumPoints "<<passednumPoints;
		opserr<< "numPoints "<<numPoints<<"\n";
		opserr<< "passednumRV"<<passednumRV;
		opserr<< "numRV "<<numRV<<"\n";
		numPoints=passednumPoints;
		numRV=passednumRV;
		this->allocate(numPoints,numRV);
	}
	for(int i=0; i<numPoints; i++){
		Lsf[i]=0;
		numAna[i]=0;
		numAnaIncSens[i]=0;
		numLinSearch[i]=0;
		numSteps[i]=0;
		iresult[i]=0;
		(*thresholdVal)(i)=0.0;
		(*betaPoints)(i)=0.0;
		(*pfPoints)(i)=0.0;
		(*nuPoints)(i)=0.0;
		(*resTime)(i)=1.0e10;
		(*check1)(i)=0.0;
		(*check2)(i)=0.0;
		(*check1_init)(i)=0.0;
		(*check2_init)(i)=0.0;
		for( int j=0; j<numRV; j++){
			(*xDesPoints[i])(j)=0.0;
			(*uDesPoints[i])(j)=0.0;
			(*alphaPoints[i])(j)=0.0;	
			(*hfuncs[i])(j)=0.0;	
		}
	}
}
void OutCrossingResults::allocate(int numPoints, int numRV)
{
	if(Lsf!=0){ delete [] Lsf; Lsf=0; }
	if(numAna!=0){ delete [] numAna; numAna=0; }
	if(numAnaIncSens!=0){ delete [] numAnaIncSens; numAnaIncSens=0;}
	if(numLinSearch!=0){delete [] numLinSearch; numLinSearch=0;}
	if(numSteps!=0){ delete [] numSteps; numSteps=0; }
	if(iresult!=0){ delete [] iresult; iresult=0; }
	if(xDesPoints!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete xDesPoints[i]; xDesPoints[i]=0; }
		delete [] xDesPoints;
		xDesPoints = 0;
	}
	if(uDesPoints!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete uDesPoints[i]; uDesPoints[i]=0; } 
		delete [] uDesPoints;
		uDesPoints = 0;
	}
	if(alphaPoints!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete alphaPoints[i]; alphaPoints[i]=0; }
		delete [] alphaPoints;
		alphaPoints=0;
	}
	if(hfuncs!=0) {
		for(int i=0;i<numPoints; i++)
		{ delete hfuncs[i]; hfuncs[i]=0; }
		delete [] hfuncs;
		hfuncs=0;
	}
	if(thresholdVal!=0){ delete thresholdVal; thresholdVal=0 ;}
	if(betaPoints!=0){ delete betaPoints; betaPoints=0 ;}
	if(pfPoints!=0){ delete pfPoints;pfPoints=0;}
	if(nuPoints!=0){ delete nuPoints; nuPoints=0;}
	if(resTime!=0){delete resTime;resTime=0;}
	if(check1!=0){delete check1;check1=0;}
	if(check2!=0){delete check2;check2=0;}
	if(check1_init!=0){delete check1_init;check1_init=0;}
	if(check2_init!=0){delete check2_init;check2_init=0;}

	Lsf= new int[numPoints];
	if( Lsf == 0){ opserr << "Fail to allocate Lsf\n";
					  opserr << "OutCrossingResults::allocate\n";
					  exit(1);}
	numAna= new int[numPoints];
	if( numAna == 0){ opserr << "Fail to allocate numAna\n";
					  opserr << "OutCrossingResults::allocate\n";
					  exit(1);}

	numAnaIncSens = new int[numPoints];
	if( numAnaIncSens == 0){ opserr << "Fail to allocate numAnaIncSens\n";
							 opserr << "OutCrossingResults::allocate\n";
							 exit(1); }

	numLinSearch = new int[numPoints];
	if( numLinSearch == 0){	opserr << "Fail to allocate numLinSearch\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1); }

	numSteps = new int[numPoints];
	if( numSteps == 0){ opserr << "Fail to allocate numLinSearch\n";
						opserr << "OutCrossingResults::allocate\n";
						exit(1); }

	iresult = new int[numPoints];
	if( iresult == 0){	opserr << "Fail to allocate iresult\n";
						opserr << "OutCrossingResults::allocate\n";
						exit(1); }

	xDesPoints = new Vector*[numPoints];
	if( xDesPoints == 0){ opserr << "Fail to allocate xDesPoints\n";
						  opserr << "OutCrossingResults::allocate\n";
						  exit(1);  }

	uDesPoints = new Vector*[numPoints];
	if( uDesPoints == 0){	opserr << "Fail to allocate xDesPoints\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1); 	}

	alphaPoints = new Vector*[numPoints];
	if( alphaPoints == 0){	opserr << "Fail to allocate alphaPoints\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1);
						}
	hfuncs = new Vector*[numPoints];
	if( hfuncs == 0){	opserr << "Fail to allocate hfuncs\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1);
						}

	for( int i=0; i<numPoints; i++){
		xDesPoints[i]= new Vector(numRV);
		if( xDesPoints[i] == 0){	opserr << "Fail to allocate xDesPoints[i]\n";
									opserr << "OutCrossingResults::allocate\n";
									exit(1); }

		uDesPoints[i]= new Vector(numRV);
		if( uDesPoints[i] == 0){	opserr << "Fail to allocate uDesPoints[i]\n";
									opserr << "OutCrossingResults::allocate\n";
									exit(1); }

		alphaPoints[i]= new Vector(numRV);
		if( alphaPoints[i] == 0){	opserr << "Fail to allocate alphaPoints[i]\n";
									opserr << "OutCrossingResults::allocate\n";
									exit(1); }
		hfuncs[i]= new Vector(numRV);
		if( hfuncs[i] == 0){	opserr << "Fail to allocate hfuncs[i]\n";
									opserr << "OutCrossingResults::allocate\n";
									exit(1); }
	}

	thresholdVal = new Vector(numPoints);
	if( thresholdVal == 0){	opserr << "Fail to allocate thresholdVal[i]\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1); }

	betaPoints = new Vector(numPoints);
	if( betaPoints == 0){	opserr << "Fail to allocate betaPoints[i]\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1); }

	pfPoints = new Vector(numPoints);
	if( pfPoints == 0){		opserr << "Fail to allocate pfPoints\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1);  }

	nuPoints = new Vector(numPoints);
	if( nuPoints == 0){		opserr << "Fail to allocate nuPoints\n";
							opserr << "OutCrossingResults::allocate\n";
							exit(1); }

	resTime = new Vector(numPoints);
	if( resTime == 0){	opserr << "Fail to allocate resTime\n";
						opserr << "OutCrossingResults::allocate\n";
						exit(1); 	}

	check1 = new Vector(numPoints);
	if( check1 == 0){	opserr << "Fail to allocate check1\n";
						opserr << "OutCrossingResults::allocate\n";
						exit(1);	}

	check2 = new Vector(numPoints);
	if( check2 == 0){	opserr << "Fail to allocate check1\n";
						opserr << "OutCrossingResults::allocate\n";
						exit(1);	}
	check1_init = new Vector(numPoints);
	if( check1_init == 0){	opserr << "Fail to allocate check1\n";
						opserr << "OutCrossingResults::allocate\n";
						exit(1);	}

	check2_init = new Vector(numPoints);
	if( check2_init == 0){	opserr << "Fail to allocate check1\n";
						opserr << "OutCrossingResults::allocate\n";
						exit(1);	}
}
void OutCrossingResults::setLsf(int indicator,int number)
{
	if(indicator<=numPoints){
		Lsf[indicator-1]=number;
	}
	else{
		opserr << "Fail to put number in Lsf\n";
		opserr << "OutCrossingResults::setnumAna\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}

void OutCrossingResults::setnumAna(int indicator,int number)
{
	if(indicator<=numPoints){
		numAna[indicator-1]=number;
	}
	else{
		opserr << "Fail to put number in numnAna\n";
		opserr << "OutCrossingResults::setnumAna\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setnumAnaIncSens(int indicator,int number)
{
	if(indicator<=numPoints){
		numAnaIncSens[indicator-1]=number;
	}
	else{
		opserr << "Fail to put number in numAnaIncSens\n";
		opserr << "OutCrossingResults::setnumAnaIncSens\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setnumLinSearch(int indicator,int number)
{
	if(indicator<=numPoints){
		numLinSearch[indicator-1]=number;
	}
	else{
		opserr << "Fail to put number in numAnaIncSens\n";
		opserr << "OutCrossingResults::setnumLinSearch\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setnumSteps(int indicator,int number)
{
	if(indicator<=numPoints){
		numSteps[indicator-1]=number;
	}
	else{
		opserr << "Fail to put number in numSteps\n";
		opserr << "OutCrossingResults::setnumSteps\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setiresult(int indicator,int number)
{
	if(indicator<=numPoints){
		iresult[indicator-1]=number;
	}
	else{
		opserr << "Fail to put number in numSteps\n";
		opserr << "OutCrossingResults::setnumSteps\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setxDesPoints(int indicator, Vector& design)
{
	if(indicator<=numPoints && design.Size()==numRV ){
		for(int i=0; i< numRV; i++){
			(*xDesPoints[indicator-1])(i)=design(i);
		}
	}
	else{
		opserr << "Fail to put vector in DesPoints\n";
		opserr << "OutCrossingResults::setDesPoints\n";
		opserr << "indicator "<< indicator << "\n";
		opserr << "size of vector "<< design.Size() << "\n";
		exit(1);
	}
}
void OutCrossingResults::setuDesPoints(int indicator, Vector& design)
{
	if(indicator<=numPoints && design.Size()==numRV ){
		for(int i=0; i< numRV; i++){
			(*uDesPoints[indicator-1])(i)=design(i);
		}
	}
	else{
		opserr << "Fail to put vector in DesPoints\n";
		opserr << "OutCrossingResults::setDesPoints\n";
		opserr << "indicator "<< indicator << "\n";
		opserr << "size of vector "<< design.Size() << "\n";
		exit(1);
	}
}
void OutCrossingResults::setDesAlpha(int indicator, Vector& design)
{
	if(indicator<=numPoints && design.Size()==numRV ){
		for(int i=0; i< numRV; i++){
			(*alphaPoints[indicator-1])(i)=design(i);
		}
	}
	else{
		opserr << "Fail to put vector in DesAlpha\n";
		opserr << "OutCrossingResults::setDesAlpha\n";
		opserr << "indicator "<< indicator << "\n";
		opserr << "size of vector "<< design.Size() << "\n";
		exit(1);
	}
}
void OutCrossingResults::setHfunc(int indicator, Vector& hfunc)
{
	if(indicator<=numPoints && hfunc.Size()==numRV ){
		for(int i=0; i< numRV; i++){
			(*hfuncs[indicator-1])(i)=hfunc(i);
		}
	}
	else{
		opserr << "Fail to put vector in hfuncs\n";
		opserr << "OutCrossingResults::setHfunc\n";
		opserr << "indicator "<< indicator << "\n";
		opserr << "size of vector "<< hfunc.Size() << "\n";
		exit(1);
	}
}
void OutCrossingResults::setbeta(int indicator, double value)
{
	if(indicator<=numPoints){
		(*betaPoints)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setbeta\n";
		opserr << "OutCrossingResults::setbeta\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setthresholdVal(int indicator, double value){

	if(indicator<=numPoints){
		(*thresholdVal)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setbeta\n";
		opserr << "OutCrossingResults::setbeta\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setpf(int indicator, double value)
{
	if(indicator<=numPoints){
		(*pfPoints)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setpf\n";
		opserr << "OutCrossingResults::setpf\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setnu(int indicator, double value)
{
	if(indicator<=numPoints){
		(*nuPoints)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setnu\n";
		opserr << "OutCrossingResults::setnu\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::settime(int indicator, double value)
{
	if(indicator<=numPoints){
		(*resTime)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setTime\n";
		opserr << "OutCrossingResults::setTime\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setcheck1(int indicator, double value)
{
	if(indicator<=numPoints){
		(*check1)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setTime\n";
		opserr << "OutCrossingResults::setTime\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setcheck2(int indicator, double value)
{
	if(indicator<=numPoints){
		(*check2)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setTime\n";
		opserr << "OutCrossingResults::setTime\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setcheck1_init(int indicator, double value)
{
	if(indicator<=numPoints){
		(*check1_init)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setTime\n";
		opserr << "OutCrossingResults::setTime\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
void OutCrossingResults::setcheck2_init(int indicator, double value)
{
	if(indicator<=numPoints){
		(*check2_init)(indicator-1)=value;
	}
	else{
		opserr << "Fail to put value in setTime\n";
		opserr << "OutCrossingResults::setTime\n";
		opserr << "indicator "<< indicator << "\n";
		exit(1);
	}
}
int OutCrossingResults::getnumLsf(void){ return numLsf; }
int OutCrossingResults::getnumFragility(void){ return numFragility; }
int OutCrossingResults::getLsf(int iord){ return Lsf[iord-1]; }
int OutCrossingResults::getnumAna(int iord){ return numAna[iord-1]; }
int OutCrossingResults::getnumAnaIncSens(int iord){	return numAnaIncSens[iord-1];}
int OutCrossingResults::getnumLinSearch(int iord){ return numLinSearch[iord-1]; }
int OutCrossingResults::getnumSteps(int iord){ return numSteps[iord-1]; }
Vector& OutCrossingResults::getxDesPoints(int iord){ return *xDesPoints[iord-1]; }
Vector& OutCrossingResults::getuDesPoints(int iord){ return *uDesPoints[iord-1]; }
Vector& OutCrossingResults::getDesAlpha(int iord){ return *alphaPoints[iord-1]; }
Vector& OutCrossingResults::getHfunc(int iord){ return *hfuncs[iord-1]; }
double OutCrossingResults::getthresholdVal(int iord){ return (*thresholdVal)(iord-1); }
double OutCrossingResults::getbeta(int iord){ return (*betaPoints)(iord-1); }
double OutCrossingResults::getpf(int iord){ return (*pfPoints)(iord-1); }
double OutCrossingResults::getnu(int iord){	return (*nuPoints)(iord-1); }
double OutCrossingResults::getTime(int iord){ return (*resTime)(iord-1); }
double OutCrossingResults::getcheck1(int iord){ return (*check1)(iord-1); }
double OutCrossingResults::getcheck2(int iord){ return (*check2)(iord-1); }
double OutCrossingResults::getcheck1_init(int iord){ return (*check1_init)(iord-1); }
double OutCrossingResults::getcheck2_init(int iord){ return (*check2_init)(iord-1); }
int OutCrossingResults::getiresult(int iord){ return iresult[iord-1]; }
void OutCrossingResults::outtoFile()
{
//	int templsf;
//	int tempiFrag;

/*	ifstream ifile;
	ifile.open( fileName, ios::in | ios::binary );  
	int nsize=17+3*(numPoints+1);

	if(ifile!=NULL){
		int ipt=0;
		while(!ifile.eof()){	
			ifile.read((char*)&templsf, sizeof(int) );
			ifile.read((char*)&tempiFrag, sizeof(int) );
			if(templsf==IdLsf&&tempiFrag==IdFragility){
				opserr<<"!! warning !!\n";
				opserr<<" lsf "<<IdLsf<<" ifrag "<<IdFragility<<"\n";
				opserr<<" is already in the file\n";
			}
			ipt+=nsize;
			ifile.seekg(ipt);
		}
		ifile.close();
	}
*/
	int np=numPoints;
	ofstream ofile;
	ofile.open( fileName, ios::out | ios::binary |ios::trunc);  
//	ofile.write((char*)&IdLsf, sizeof(int) );				//01
//	ofile.write((char*)&IdFragility, sizeof(int) );			//02
//	ofile.write((char*)&threshold, sizeof(double) );		//03
//	ofile.write((char*)&factorFragility, sizeof(double) );	//04
	ofile.write((char*)&numPoints, sizeof(int) );			//05
	ofile.write((char*)&numLsf, sizeof(int) );			//05
	ofile.write((char*)&numFragility, sizeof(int) );			//05
 	ofile.write((char*)&numRV, sizeof(int) );	
	int itemp;
	double dtemp;//06
	for(int i=0; i<np; i++){itemp=numSteps[i]; ofile.write((char*)&itemp, sizeof(int));}
//	ofile.write((char*)&numSteps, sizeof(int)*np );			//07
	for(int i=0; i<np; i++){itemp=Lsf[i]; ofile.write((char*)&itemp, sizeof(int));}
//	ofile.write((char*)&Lsf, sizeof(int)*np );			//08
	for(int i=0; i<np; i++){itemp=numAna[i]; ofile.write((char*)&itemp, sizeof(int));}
//	ofile.write((char*)&numAna, sizeof(int)*np );			//08
	for(int i=0; i<np; i++){itemp=numAnaIncSens[i]; ofile.write((char*)&itemp, sizeof(int));}
//	ofile.write((char*)&numAnaIncSens, sizeof(int)*np );	//09	
	for(int i=0; i<np; i++){itemp=numLinSearch[i]; ofile.write((char*)&itemp, sizeof(int));}
//	ofile.write((char*)&numLinSearch, sizeof(int)*np );		//10
	for(int i=0; i<np; i++){itemp=iresult[i]; ofile.write((char*)&itemp, sizeof(int));}
//	ofile.write((char*)&iresult, sizeof(int)*np );			//11
	for(int i=0; i<np; i++){dtemp=(*thresholdVal)(i); ofile.write((char*)&dtemp, sizeof(double));}
//	ofile.write((char*)&((*thresholdVal)(0)), sizeof(double)*np );			//12
	for(int i=0; i<np; i++){dtemp=(*betaPoints)(i); ofile.write((char*)&dtemp, sizeof(double));}
//	ofile.write((char*)&((*betaPoints)(0)), sizeof(double)*np );			//12
	for(int i=0; i<np; i++){dtemp=(*pfPoints)(i); ofile.write((char*)&dtemp, sizeof(double));}
//	ofile.write((char*)&((*pfPoints)(0)), sizeof(double)*np );			//13
	for(int i=0; i<np; i++){dtemp=(*nuPoints)(i); ofile.write((char*)&dtemp, sizeof(double));}
//	ofile.write((char*)&((*nuPoints)(0)), sizeof(double)*np );			//14
	for(int i=0; i<np; i++){dtemp=(*resTime)(i); ofile.write((char*)&dtemp, sizeof(double));}
//	ofile.write((char*)&((*resTime)(0)), sizeof(double)*np );			//15
	for(int i=0; i<np; i++){dtemp=(*check1)(i); ofile.write((char*)&dtemp, sizeof(double));}
//	ofile.write((char*)&((*check1)(0)), sizeof(double)*np );			//16
	for(int i=0; i<np; i++){dtemp=(*check2)(i); ofile.write((char*)&dtemp, sizeof(double));}
//	ofile.write((char*)&((*check2)(0)), sizeof(double)*np );			//17

	for(int i=0; i<numPoints; i++){
		for(int j=0; j<numRV; j++)
		{dtemp=(*xDesPoints[i])(j); ofile.write((char*)&dtemp, sizeof(double));}
//		ofile.write((char*)&(*xDesPoints[i])(0), sizeof(double)*numRV );
		for(int j=0; j<numRV; j++)
		{dtemp=(*uDesPoints[i])(j); ofile.write((char*)&dtemp, sizeof(double));}
//		ofile.write((char*)&(*uDesPoints[i])(0), sizeof(double)*numRV );
		for(int j=0; j<numRV; j++)
		{dtemp=(*alphaPoints[i])(j); ofile.write((char*)&dtemp, sizeof(double));}
//		ofile.write((char*)&(*alphaPoints[i])(0), sizeof(double)*numRV );
		for(int j=0; j<numRV; j++)
		{dtemp=(*hfuncs[i])(j); ofile.write((char*)&dtemp, sizeof(double));}
//		ofile.write((char*)&(*hfuncs[i])(0), sizeof(double)*numRV );
	}
	ofile.close();
}
void OutCrossingResults::readfromFile(void)
{
	/// read in from the file /////
	int tempnumPoints;
	int tempnumRV;
	int tempnumLsf;
	int tempnumFragility;

	ifstream ifile;
	ifile.open( fileName, ios::in | ios::binary );  
	if(!ifile){
		opserr<< "No file " << fileName <<" exists\n";
		exit(-1);
	}
	ifile.read((char*)&tempnumPoints, sizeof(int) );
	ifile.read((char*)&tempnumLsf, sizeof(int) );
	ifile.read((char*)&tempnumFragility, sizeof(int) );
	ifile.read((char*)&tempnumRV, sizeof(int) );
	numPoints=tempnumPoints;
	numLsf=tempnumLsf;
	numFragility=tempnumFragility;
	numRV=tempnumRV;
	this->allocate(numPoints,numRV);

	int tempint;
	double tempdbl;
	int np=numPoints;
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempint, sizeof(int) );numSteps[i]=tempint;}		//07
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempint, sizeof(int) );Lsf[i]=tempint;}		//07
//	ifile.read((char*)&Lsf, sizeof(int)*np );			//08
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempint, sizeof(int) );numAna[i]=tempint;}		//07
//	ifile.read((char*)&numAna, sizeof(int)*np );			//08
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempint, sizeof(int) );numAnaIncSens[i]=tempint;}		//07
//	ifile.read((char*)&numAnaIncSens, sizeof(int)*np );	//09	
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempint, sizeof(int) );numLinSearch[i]=tempint;}		//07
//	ifile.read((char*)&numLinSearch, sizeof(int)*np );		//10
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempint, sizeof(int) );iresult[i]=tempint;}		//07
//	ifile.read((char*)&iresult, sizeof(int)*np );			//11
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempdbl, sizeof(double) );(*thresholdVal)(i)=tempdbl;}		//07
//	double& vvv=(*thresholdVal)(0);
//	ifile.read((char*)&vvv, sizeof(double)*np );			//12
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempdbl, sizeof(double) );(*betaPoints)(i)=tempdbl;}		//07
//	vvv=(*betaPoints)(0);
//	ifile.read((char*)&vvv, sizeof(double)*np );			//12
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempdbl, sizeof(double) );(*pfPoints)(i)=tempdbl;}		//07
//	vvv=(*pfPoints)(0);
//	ifile.read((char*)&vvv, sizeof(double)*np );			//13
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempdbl, sizeof(double) );(*nuPoints)(i)=tempdbl;}		//07
//	vvv=(*nuPoints)(0);
//	ifile.read((char*)&vvv, sizeof(double)*np );			//14
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempdbl, sizeof(double) );(*resTime)(i)=tempdbl;}		//07
//	vvv=(*resTime)(0);
//	ifile.read((char*)&vvv, sizeof(double)*np );			//15
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempdbl, sizeof(double) );(*check1)(i)=tempdbl;}		//07
//	vvv=(*check1)(0);
//	ifile.read((char*)&vvv, sizeof(double)*np );			//16
	for(int i=0; i<np; i++)
	{ifile.read((char*)&tempdbl, sizeof(double) );(*check2)(i)=tempdbl;}		//07
//	vvv=(*check2)(0);
//	ifile.read((char*)&vvv, sizeof(double)*np );			//17

	for(int i=0; i<numPoints; i++){
		for(int j=0; j<numRV; j++)
		{ifile.read((char*)&tempdbl, sizeof(double) );(*xDesPoints[i])(j)=tempdbl;}		//07
//		vvv=(*xDesPoints[i])(0);
//		ifile.read((char*)&vvv, sizeof(double)*numRV );
		for(int j=0; j<numRV; j++)
		{ifile.read((char*)&tempdbl, sizeof(double) );(*uDesPoints[i])(j)=tempdbl;}		//07
//		vvv=(*uDesPoints[i])(0);
//		ifile.read((char*)&vvv, sizeof(double)*numRV );
		for(int j=0; j<numRV; j++)
		{ifile.read((char*)&tempdbl, sizeof(double) );(*alphaPoints[i])(j)=tempdbl;}		//07
//		vvv=(*alphaPoints[i])(0);
//		ifile.read((char*)&vvv, sizeof(double)*numRV );
		for(int j=0; j<numRV; j++)
		{ifile.read((char*)&tempdbl, sizeof(double) );(*hfuncs[i])(j)=tempdbl;}		//07
//		vvv=(*hfuncs[i])(0);
//		ifile.read((char*)&vvv, sizeof(double)*numRV );
	}
	ifile.close();

}
void OutCrossingResults::printResults(ofstream& outputFile)
{

	outputFile<<"\n";
	outputFile<<"  Lsf"<<"      threshold";
	outputFile<<" step"<<"      Time"<<"           beta";
	outputFile<<"             pf"<<"             nu";
	outputFile<<"    numAna"<<" numIncSns"<<" numLinSer";
	outputFile<<"         check1"<<"         check2"<<"   iresult\n";
	outputFile.setf(ios::right);
	outputFile.setf(ios::scientific, ios::floatfield);
	outputFile.flush();

	for(int i=0;i<numPoints;i++){
		int istep=numSteps[i];
		if( istep==0) continue;
		double time=(*resTime)(i);
		outputFile << setw(5) << Lsf[i];
		outputFile << setw(15) << setprecision(5) 
				   << (*thresholdVal)(i);
		outputFile << setw(5) << istep;
		outputFile << setw(10) << setprecision(2) << time;
		outputFile << setw(15) << setprecision(5) 
				   << (*betaPoints)(i);
		outputFile << setw(15) << setprecision(5) 
				   << (*pfPoints)(i);
		outputFile << setw(15) << setprecision(5) 
				   << (*nuPoints)(i);
		outputFile << setw(10) << numAna[i];
		outputFile << setw(10) << numAnaIncSens[i];
		outputFile << setw(10) << numLinSearch[i];
		outputFile << setw(15) << setprecision(5) 
				   << (*check1)(i);
		outputFile << setw(15) << setprecision(5) 
				   << (*check2)(i);
		outputFile << setw(10) << iresult[i];
		outputFile << setw(15) << setprecision(5) 
				   << (*check1_init)(i);
		outputFile << setw(15) << setprecision(5) 
				   << (*check2_init)(i);
		outputFile << "\n";
	}
	outputFile.flush();
}
void OutCrossingResults::printSinglePoint(ofstream& outputFile,
										  int ipt)
{
//	if(passedlsf!=IdLsf){
//		opserr<<"fatal in OutCrossingResults::printSinglePoint\n";
//		opserr<<"passedlsf "<<passedlsf;
//		opserr<<"IdLsf "<<IdLsf;
//		exit(-1);
//	}
	outputFile.setf(ios::right);
	int istep=numSteps[ipt-1];
	if( istep==0) {
		opserr<<"No results found for point"<<ipt<<"\n";
		opserr<<"OutCrossingResults::printSinglePoint\n";
	}
	double time=(*resTime)(ipt-1);

	outputFile << setw(5) <<Lsf[ipt-1];
//	outputFile.setf(ios::fixed, ios::floatfield);
//	outputFile << setw(10)<<setprecision(2)<<factorFragility;
	outputFile.setf(ios::scientific, ios::floatfield);
	outputFile << setw(15)<<setprecision(5)<<(*thresholdVal)(ipt-1);
	outputFile << setw(5) <<istep;
	outputFile.setf(ios::fixed, ios::floatfield);
	outputFile << setw(10)<<setprecision(2)<<time;
	outputFile.setf(ios::scientific, ios::floatfield);
	outputFile << setw(15) << setprecision(5) << (*betaPoints)(ipt-1);
	outputFile << setw(15) << setprecision(5) << (*pfPoints)(ipt-1);
	outputFile << setw(15) << setprecision(5) << (*nuPoints)(ipt-1);
	outputFile << setw(10) << numAna[ipt-1];
	outputFile << setw(10) << numAnaIncSens[ipt-1];
	outputFile << setw(10) << numLinSearch[ipt-1];
	outputFile << setw(15) << setprecision(5) << (*check1)(ipt-1);
	outputFile << setw(15) << setprecision(5) << (*check2)(ipt-1);
	outputFile << setw(10) << iresult[ipt-1];
	outputFile << setw(15) << setprecision(5) << (*check1_init)(ipt-1);
	outputFile << setw(15) << setprecision(5) << (*check2_init)(ipt-1);
	outputFile << "\n";
//	outputFile << " DESIGN POINT INFORMATION \n";
//	outputFile << "\n";
//	outputFile << "      step";
//	outputFile << "                xdes";
//	outputFile << "                udes";
//	outputFile << "               alpha";
//	int nSize=xDesPoints[ipt]->Size();
//	outputFile.setf(ios::scientific, ios::floatfield);
//  for( int i=0; i<nSize; i++){
//    outputFile << setw(10) << i;
//	outputFile << setw(20) << setprecision(10) << (*xDesPoints[ipt])(i);
//	outputFile << setw(20) << setprecision(10) << (*uDesPoints[ipt])(i);
//	outputFile << setw(20) << setprecision(10) << (*alphaPoints[ipt])(i);
//	outputFile << "\n";
//	}
	outputFile.flush();
}
void OutCrossingResults::printSingleTitle(ofstream& outputFile,int passedlsf)
{
	outputFile << "\n";
	outputFile << "-------------------------------------------------------\n";
	outputFile << "----- Design Point Results for Lsf "<<passedlsf<<"\n";
	outputFile << "-------------------------------------------------------\n";
	outputFile << "\n";
	outputFile << "  Lsf";
	outputFile << "      Threshold";
	outputFile << " step";
	outputFile << "      time";
	outputFile << "           beta";
	outputFile << "             pf";
	outputFile << "             nu";
	outputFile << "    numAna";
	outputFile << "   IncSens";
	outputFile << "   LinSear";
	outputFile << "         check1";
	outputFile << "         check2";
	outputFile << "    iresult";
	outputFile << "\n";
}
/*
int OutCrossingResults::ordering()
{
	ordered=true;

	if((*resTime)(0)!=0.0){ 
		for(int i=0;i<numPoints-1;i++){
			double timei=(*resTime)(i);
			if(timei==0.0) timei=1.0e10;
			for( int j=i+1; j<numPoints; j++){
				double timej=(*resTime)(j);
				if(timej==0.0) timej=1.0e10;
				if(timej<timei){
					(*resTime)(i)=timei;
					(*resTime)(j)=timej;
					int ntmp=numAna[i];
					numAna[i]=numAna[j];
					numAna[j]=ntmp;
					ntmp=numAnaIncSens[i];
					numAnaIncSens[i]=numAnaIncSens[j];
					numAnaIncSens[j]=ntmp;
					ntmp=numLinSearch[i];
					numLinSearch[i]=numLinSearch[j];
					numLinSearch[j]=ntmp;
					ntmp=numSteps[i];
					numSteps[i]=numSteps[j];
					numSteps[j]=ntmp;
					ntmp=iresult[i];
					iresult[i]=iresult[j];
					iresult[j]=ntmp;
					double tmp=(*betaPoints)(i);
					(*betaPoints)(i)=(*betaPoints)(j);
					(*betaPoints)(j)=tmp;
					tmp=(*pfPoints)(i);
					(*pfPoints)(i)=(*pfPoints)(j);
					(*pfPoints)(j)=tmp;
					tmp=(*nuPoints)(i);
					(*nuPoints)(i)=(*nuPoints)(j);
					(*nuPoints)(j)=tmp;
					tmp=(*check1)(i);
					(*check1)(i)=(*check1)(j);
					(*check1)(j)=tmp;
					tmp=(*check2)(i);
					(*check2)(i)=(*check2)(j);
					(*check2)(j)=tmp;
					Vector* vtmp=xDesPoints[i];
					xDesPoints[i]=xDesPoints[j];
					xDesPoints[j]=vtmp;
					vtmp=uDesPoints[i];
					uDesPoints[i]=uDesPoints[j];
					uDesPoints[j]=vtmp;
					vtmp=alphaPoints[i];
					alphaPoints[i]=alphaPoints[j];
					alphaPoints[j]=vtmp;

//					Vector vtmp=(*xDesPoints[i]);
//					(*xDesPoints[i])=(*xDesPoints[j]);
//					(*xDesPoints[j])=vtmp;
//					vtmp=(*uDesPoints[i]);
//					(*uDesPoints[i])=(*uDesPoints[j]);
//					(*uDesPoints[j])=vtmp;
//					vtmp=(*alphaPoints[i]);
//					(*alphaPoints[i])=(*alphaPoints[j]);
//					(*alphaPoints[j])=vtmp;
				}
			}
		}
		numPoints++;
		return numPoints; 
	}else{
		for(int i=0;i<numPoints;i++){
			int inext=i+1;
			(*resTime)(i)=(*resTime)(inext);
			numAna[i]=numAna[inext];
			numAnaIncSens[i]=numAnaIncSens[inext];
			numLinSearch[i]=numLinSearch[inext];
			numSteps[i]=numSteps[inext];
			iresult[i]=iresult[inext];
			(*betaPoints)(i)=(*betaPoints)(inext);
			(*pfPoints)(i)=(*pfPoints)(inext);
			(*nuPoints)(i)=(*nuPoints)(inext);
			(*check1)(i)=(*check1)(inext);
			(*check2)(i)=(*check2)(inext);
			xDesPoints[i]=xDesPoints[inext];
			uDesPoints[i]=uDesPoints[inext];
			alphaPoints[i]=alphaPoints[inext];
		}
		(*resTime)(numPoints)=0.0;
		numAna[numPoints]=0;
		numAnaIncSens[numPoints]=0;
		numLinSearch[numPoints]=0;
		numSteps[numPoints]=0;
		iresult[numPoints]=0;
		(*betaPoints)(numPoints)=0.0;
		(*pfPoints)(numPoints)=0.0;
		(*nuPoints)(numPoints)=0.0;
		(*check1)(numPoints)=0.0;
		(*check2)(numPoints)=0.0;
		xDesPoints[numPoints]->Zero();
		uDesPoints[numPoints]->Zero();
		alphaPoints[numPoints]->Zero();
		return numPoints;
	}
}
*/
