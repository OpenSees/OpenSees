//////////////////////////////////////////////////////////////////////
//  This file contains the constructor, destructor, and member	    //  
//  functions for the ShallowFoundationGen class.                   //
//  								    //
//  Written by: Prishati Raychowdhury (prishati@gmail.com)	    //
//              Graduate Student, UC San Diego			    //
//              November, 2007					    //
//////////////////////////////////////////////////////////////////////

/***************************************************************************************************
** Copyright @ 2007 The Regents of the University of California (The Regents). All Rights Reserved.
**
** The Regents grants permission, without fee and without a written license agreement, for (a) use, 
** reproduction, modification, and distribution of this software and its documentation by educational, 
** research, and non-profit entities for noncommercial purposes only; and (b) use, reproduction and 
** modification of this software by other entities for internal purposes only. The above copyright 
** notice, this paragraph and the following three paragraphs must appear in all copies and modifications 
** of the software and/or documentation.
**
** Permission to incorporate this software into products for commercial distribution may be obtained 
** by contacting the University of California 
** Office of Technology Licensing 
** 2150 Shattuck Avenue #510, 
** Berkeley, CA 94720-1620, 
** (510) 643-7201.
**
** This software program and documentation are copyrighted by The Regents of the University of California. 
** The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The 
** end-user understands that the program was developed for research purposes and is advised not to rely 
** exclusively on the program for any reason.
**
** IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 
** CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
** DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. REGENTS GRANTS 
** NO EXPRESS OR IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN INDIVIDUAL 
** CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES PROJECT AT THE UNIVERSITY OF CALIFORNIA, BERKELEY 
** TO BENEFIT THE END USER.
**
** REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
** OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
** IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, 
** SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
****************************************************************************************************/

#include "ShallowFoundationGen.h"
#include <stdlib.h>
#include <string>
using namespace std;

#include <elementAPI.h>
#include <sstream>

int OPS_ShallowFoundationGen()
{
    if (OPS_GetNumRemainingInputArgs() < 4) {
	opserr << "WARNING ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType?";
	opserr << "Must have 4 arguments." << endln;
	return -1;
    }

    int tags[2];
    int num = 2;
    if (OPS_GetIntInput(&num,tags) < 0) {
	opserr<<"WARNING: invalid integer input\n";
	return -1;
    }

    const char* filename = OPS_GetString();

    int ftype;
    num = 1;
    if (OPS_GetIntInput(&num,&ftype) < 0) {
	opserr<<"WARNING: invalid integer input\n";
	return -1;
    }

    std::stringstream ss;
    ss << tags[0] << " " << tags[1] << " " << ftype;
    std::string id, cnode, foundtype;
    ss >> id >> cnode >> foundtype;

    ShallowFoundationGen gen;
    gen.GetShallowFoundation(id.c_str(), cnode.c_str(), filename, foundtype.c_str());

    return 0;
}

////////////////////////////////////////////////////////////////////////
// Constructor initializes global variables to zero
ShallowFoundationGen::ShallowFoundationGen()
{	
}

/////////////////////////////////////////////////////////////////////////
// Destructor deletes dynamically allocated arrays
ShallowFoundationGen::~ShallowFoundationGen()
{
}

///////////////////////////////////////////////////////////////////////////////////////
// Function to get inputfile and writing output tcl source file
void ShallowFoundationGen::GetShallowFoundation(const char *FoundationID, const char *ConnectingNode, const char *InputDataFile, const char *FoundationCondition)
{
	 int FoundationTag;
	 int ConNode;
	 int FootingCondition;

	 FoundationTag=atoi(FoundationID);
	 ConNode=atoi(ConnectingNode);
	 FootingCondition=atoi(FoundationCondition);

	
    // Define local variables
	string str1("Foundation_");
    string str2(FoundationID);
    string str3(".tcl");
	string str4;
	str4= (str1+str2+str3);
    const char *OutFile;
    OutFile=str4.c_str();
   
	ofstream ShallowFoundationOut(OutFile,ios::out); // Initialize output stream

	// Write headers for output file
	ShallowFoundationOut << "######################################################################################" << endln;
	ShallowFoundationOut << "#                                                                                    #" << endln;
	ShallowFoundationOut << "# This is an intermediate file generated by the command ShallowFoundationGen.        #" << endln; 
	ShallowFoundationOut << "# Source it after the ShallowFoundationGen command.                                  #" << endln;
	ShallowFoundationOut << "# Use this file to check shallow foundation nodes, elements,  fixity details         #" << endln;
	ShallowFoundationOut << "# ShallowFoundationGen.cpp is developed by Prishati Raychowdhury (UCSD)              #" << endln;
	ShallowFoundationOut << "#                                                                                    #" << endln;
	ShallowFoundationOut << "######################################################################################" << endln << endln;
	ShallowFoundationOut << endln;
	ShallowFoundationOut << " # Foundation Tag ="  << FoundationTag <<endln;
	ShallowFoundationOut << " # Foundation Base Condition Tag ="  << FootingCondition <<endln;
    ShallowFoundationOut << endln;

	//Opening input file for reading	
	ifstream inFile(InputDataFile, ios::in);
	
	if(!inFile)
	{
		opserr << "File " << InputDataFile << "does not exist.  Must exit." << endln;
		exit(-1);
	}
	
   string tmpString; //defining temporary string

    int SoilType; double cSoil; double PhiSoil; double GammaSoil; double GSoil; double NuSoil;
	double CradSoil; double TpSoil; int CapTag; double Wg;
	double Qult; double Pult; double Tult; double Kv; double Kx; double Lfoot;
	double Bfoot; double Hfoot; double Dfoot; double Efoot; 
	double Rk; double Re; double leRatio; double beta;
	char charBuffer[256];

//individually read lines from the input file, into tmpString class
	while( !inFile.eof() ) //while not at end of file
	{
		inFile >> tmpString;
		if(tmpString[0] == '#') 
		{
		   // Skip comment lines
		}
		//Reading soil property data
		else if( tmpString.compare("SoilProp")==0 ) 	
		{
			inFile >> SoilType;
			inFile >> cSoil;
			inFile >> PhiSoil;
			inFile >> GammaSoil;
			inFile >> GSoil;
			inFile >> NuSoil;
			inFile >> CradSoil;
			inFile >> TpSoil; 
			CapTag =0;  //setting tag for further input
		}

       	//Reading soil capacity properties
		else if( tmpString.compare("CapSoil")==0 )
		{
			inFile >> Qult;
			inFile >> Pult;
			inFile >> Tult;
			inFile >> Kv;
			inFile >> Kx;
			CapTag = 1; //setting tag for not calculating Qult
		}
		    	
		//Reading Footing properties
		else if( tmpString.compare("FootProp")==0 )
		{
			inFile >> Lfoot;
			inFile >> Bfoot;
			inFile >> Hfoot;
			inFile >> Dfoot;
			inFile >> Efoot;
			inFile >> Wg;
			inFile >> beta;
		}

		//Reading Footing Mesh properties
		else if( tmpString.compare("MeshProp")==0 )
		{
			inFile >> Rk;
			inFile >> Re;
			inFile >> leRatio;
		}

		inFile.getline(charBuffer, 256); // purge the line
	}
	inFile.close(); //Closing input file

	//Calculate capacities (Qult, Pult and Tult); if not given


if(CapTag == 0) 
	{
double pi; double phirad; double tph; double Nphi; double q; double Nq; double Nc; double Ngamma;
double Fcs; double Fqs; double Fgammas; double Fcd; double Fqd; double Fgammad; double Fci; 
double Fqi; double Fgammai; double qu;  
double Kp; double delta; 

//-------- CALCULATION OF Qult (VERTICAL BEARING CAPACITY)
pi= 3.14159265;
phirad = PhiSoil*pi/180;
tph= tan(phirad);
Nphi=pow (tan(pi/4+phirad/2),2);
q= GammaSoil*Dfoot;

// Bearing capacity factors (Meyerhof, 1963)
Nq= Nphi*exp(pi*tph);   
if (PhiSoil == 0)
{
Nc = 5.7;  //using table 3.2 (Das, 2006)
}
else
{
Nc= (Nq-1)/tph;   //(Meyerhof, 1963)
}
Ngamma= (Nq-1)*tan(1.4*phirad);  //(Meyerhof, 1963)
// Shape factors (Meyerhof, 1963)
if (PhiSoil == 0)
{
Fcs= 1+ 0.2*(Bfoot/Lfoot);
Fqs= 1;
Fgammas= Fqs;
}
else
{
Fcs= 1+ 0.2*(Bfoot/Lfoot)*Nphi;
Fqs= 1+ 0.1*(Bfoot/Lfoot)*Nphi;
Fgammas= Fqs;
}
// Depth factors (Meyerhof, 1963)
if (PhiSoil == 0)
{
Fcd= 1+ 0.2*(Dfoot/Bfoot);
Fqd= 1;
Fgammad= Fqd;
}
else 
{
Fcd= 1+ 0.2*(Dfoot/Bfoot)*(sqrt(Nphi));
Fqd= 1+ 0.1*(Dfoot/Bfoot)*(sqrt(Nphi));
Fgammad= Fqd;
}
// Inclination factors (Meyerhof, 1963)
Fci= pow ((1.0-beta/90.0),2);
Fqi= Fci;

if (PhiSoil == 0) 
{
Fgammai= 1.0;
}
else
{
Fgammai= pow ((1.0-beta/PhiSoil),2); 
}
  
//Ultimate bearing capacity
qu= cSoil*Nc*Fcs*Fcd*Fci + q*Nq*Fqs*Fqd*Fqi + 0.5*GammaSoil*Bfoot*Ngamma*Fgammas*Fgammad*Fgammai;
// Ultimate load capacity
Qult = qu*Lfoot*Bfoot; 
//----------- END CALCULATION OF BEARING CAPACITY


//-------- CALCULATION OF Pult (LATERAL PASSIVE CAPACITY)
Kp=Nphi;                                    //following Rankin's passive earth pressure theory (Rankin, 1857)
Pult=(0.5*Kp*GammaSoil*Dfoot*Dfoot+2*cSoil*Dfoot*sqrt(Kp))*Lfoot;  //passive capacity of the footing 
//-------- END OF Pult CALCULATION

//-------- CALCULATION OF Tult (LATERAL SLIDING CAPACITY)
delta=0.667*phirad;                //assuming that for concrete footing, the friction angle delta=2/3 of phi
Tult=(cSoil+Wg*tan(delta))*Lfoot*Bfoot; //sliding capacity of the footing 
//-------- END OF Tult CALCULATION


//-------- CALCULATION OF STIFFNESS
//Stiffness calculation (Gazetas 1991)
// Horizontal soil stiffness
Kx=((GSoil*Lfoot)/(2-NuSoil))*(2+2.5*(pow((Bfoot/Lfoot),0.85)));
// Vertical soil stiffness
Kv=((GSoil*Lfoot)/(1-NuSoil))*(0.73+1.54*(pow((Bfoot/Lfoot),0.75)));
//-------- END OF STIFFNESS CALCULATION
	}

//------MESH GENERATION FOR SHALLOW FOUNDATION
	double Lmid; double LMidEleTrial; double LmidEle; double LendEle; int midNode; 
	double nodeTotal; int ndivMid; int ndivEnd; int i; double midNodeX; 
	double midNodeY; double LEndMidFac; int ndiv; 
	double Aelefoot; double Ielefoot;  
	double qmid; double qend; double qendext;
	double kvmid; double kvend; double kvMidSpring; double kvEndSpring; double kvEndExt;
	double Afoot; double Ifoot; 

LEndMidFac = 0.5;              // ratio of length of end element to mid region element
Lmid  =Lfoot*(1-2*Re);         //length of mid region
LMidEleTrial = leRatio*Lfoot;       //Trial element length at mid region
ndivMid = (int(Lmid/(2*LMidEleTrial)+0.5))*2;              //Calculating total number of elements excluding end elements
ndivEnd= (int(Re*Lfoot/(2*LEndMidFac*LMidEleTrial) + 0.5))*2;  //Calculating total number of elements at each end region
ndiv = ndivEnd*2+ndivMid;      //Calculating total number of elements including end elements (even number forced)
LmidEle = Lmid/(ndivMid);      //element length at mid region 
LendEle = Re*Lfoot/(ndivEnd);      //element length at end region 
nodeTotal = ndiv+1;            //Total number of nodes
midNode = ndiv/2+1;            //mid node

//Properties for each footing and spring elements
Afoot=Lfoot*Bfoot;                 //area of entire footing
Ifoot=Bfoot*pow(Lfoot,3)/12.0;     //inertia of each footing element

Aelefoot=(Bfoot*Hfoot);                      //area of each footing element -mid and end
Ielefoot=Bfoot*pow(Hfoot,3)/(12.0);          //inertia of each footing element - mid and end


qmid=Qult*(LmidEle/Lfoot);    //vertical capacity of each spring at mid region (unit=load/length)
qend=Qult*(LendEle/Lfoot);    //vertical capacity of each spring at end region (unit=load/length)
qendext=qend/2.0;             //vertical capacity of each spring at extreme two ends (unit=load/length)
kvmid=Kv/(Lmid+(Lfoot-Lmid)*Rk);   //stiffness intensity at mid region
kvend=kvmid*Rk;                     //stiffness intensity at end region  
kvMidSpring=kvmid*LmidEle;        //component stiffness of mid region springs
kvEndSpring=kvend*LendEle;		   //component stiffness of end region springs
kvEndExt=kvend*LendEle*0.5;       //component stiffness of two extreme end springs

//nodes at the footing level starts from FoundationTag*1000 + 1
//nodes at the spring level starts from FoundationTag*10,000 + 1

double Xcoord; double Ycoord; 

// calculating node coordinates
midNodeX = 0.0;  //Coordinate of mid node of foundation along X
midNodeY = 0.0;  //Coordinate of mid node of foundation along Y

//Writing node and coordinate information
ShallowFoundationOut << " #node  " << " $NodeTag " << " $Xcoord " << " $Ycoord " <<endln;

for(i=1; i<=nodeTotal; i++)
{
	if((i >= 1) && (i <= ndivEnd+1)) {
	Xcoord = midNodeX - Lfoot/2 + (i-1)*LendEle; // for end nodes in left side 
    }
	else if((i >= nodeTotal- ndivEnd) && (i <= nodeTotal)) {
	Xcoord = midNodeX + Lfoot/2 +(i-nodeTotal)*LendEle; //for end nodes in right side 
    }
    else if((i > ndivEnd+1) && (i < nodeTotal-ndivEnd)){
	Xcoord = midNodeX - Lfoot/2 + LendEle*ndivEnd + (i-ndivEnd-1)*LmidEle; // for middle nodes 
    }
	Ycoord = midNodeY; // Y coordinate 

	ShallowFoundationOut << " node  " << FoundationTag*1000+i << "  " << Xcoord << " " << Ycoord <<endln;
	ShallowFoundationOut << " node  " << FoundationTag*100000+i << " " << Xcoord << " " << Ycoord <<endln;
	
	// For horizontal springs for foot conditions 3 
	if ((i==nodeTotal) && (FootingCondition == 3)) {
	ShallowFoundationOut << " node  " << FoundationTag*100000+i+1 << " " << Xcoord << " " << Ycoord <<endln;
	}
	// For horizontal springs for foot conditions 5 
	if ((i==nodeTotal) && (FootingCondition == 5)) {
	ShallowFoundationOut << " node  " << FoundationTag*100000+i+1 << " " << Xcoord << " " << Ycoord <<endln;
	ShallowFoundationOut << " node  " << FoundationTag*100000+i+2 << " " << Xcoord << " " << Ycoord <<endln;
	}
}




//EqualDOF command
ShallowFoundationOut <<	endln << " #equalDOF $rNodeTag $cNodeTag $dof1 $dof2 $dof3" <<endln;
ShallowFoundationOut <<	" equalDOF "<< ConNode << "  " << FoundationTag*1000+midNode<< " 1 2 3 "<< endln;

//ShallowFoundationOut << endln <<" set midNode_"<<ConNode << "   " << FoundationTag*1000+midNode <<endln;

//-------END OF MESH GENERATION 


//--------- Z50 AND X50 CALCULATION
double kvfactor; double kxpfactor; double kxtfactor; double z50mid; double z50end; 
double z50endext; double xp50; double xt50; double Cd;

if 	(SoilType == 1) //soilType=clay
{
kvfactor = 0.525;    //for QzSimple2 material
kxpfactor = 8.0;     //for PySimple2 material
kxtfactor = 0.708;   //for TzSimple2 material
} 
else {               //soilType=sand
kvfactor = 1.39;     //for QzSimple2 material
kxpfactor = 0.542;   //for PySimple2 material
kxtfactor = 2.05;    //for TzSimple2 material
}
z50mid=kvfactor*qmid/kvMidSpring;
z50end=kvfactor*qend/kvEndSpring;
z50endext=kvfactor*qendext/kvEndExt;
xp50=kxpfactor*Pult/Kx;
xt50=kxtfactor*Tult/Kx;
Cd=0.1;
//--------- END OF Z50 AND X50 CALCULATION


//------MATERIAL GENERATION 
/*
/Material tag for foundation materials start from 21, assuming that the structure has used maximum 20 material types
Creating FIVE different foundation conditions (detailed description in the document)
Condition 1: fixed base case
Condition 2: elastic vertical springs; sliding restrained
Condition 3: elastic vertical sand lateral springs
Condition 4: nonlinear qz vertical springs; sliding restrained 
Condition 5: nonlinear vertical and lateral springs
*/

ShallowFoundationOut << endln << " #Materials for shallow foundation"<<endln; 

//writing foundation conditions 
if(FootingCondition == 1)  //fixed base case
	{
ShallowFoundationOut << endln <<" #fix " << " $ConnectingNode " << " " << 1 << " " << 1 << " "<< 1 << endln; 
ShallowFoundationOut << " fix " << ConNode << "  " << 1 << " " << 1 << " " << 1 << endln; 

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $KvendExtreme " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+1 <<"  " <<  kvEndExt << endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $Kvend " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+2 <<"  " <<  kvEndSpring << endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $Kvmid " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+3 <<"  "<<  kvMidSpring << endln;

	}

if(FootingCondition == 2)  //elastic vertical springs; sliding restrained
	{
ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $KvendExtreme " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+1 <<"  " <<  kvEndExt << endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $Kvend " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+2 <<"  " <<  kvEndSpring << endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $Kvmid " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+3 <<"  " <<  kvMidSpring << endln;

ShallowFoundationOut << endln <<" #fix " << " $midNode " << " " << 1 << " " << 0 << " "<< 0 << endln; 
ShallowFoundationOut << " fix " << ConNode << " " << 1 << " " << 0 << " " << 0 << endln; 
	}

if(FootingCondition == 3)    //elastic vertical sand lateral springs
	{

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $KvendExtreme " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+1 <<"  " <<  kvEndExt << endln;		

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $Kvend " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+2 <<"  " <<  kvEndSpring << endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $Kvmid " << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+3 <<"  " <<  kvMidSpring << endln;
	
ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " Elastic " << " $matTag " << " $Kx" << endln;
ShallowFoundationOut << " uniaxialMaterial " << " Elastic " << FoundationTag*100+4 <<"  " <<  Kx << endln;
	}

if(FootingCondition == 4)    //nonlinear qz vertical springs; sliding restrained
	{		

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " QzSimple2 " << " $matTag " << " $SoilType " << " $Qult-end-extreme " << " $z50-end "  << " <TpSoil> " << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " QzSimple2 " << FoundationTag*100+1 <<"  " << SoilType << " " << qendext << " " << z50endext << " " << TpSoil << " " << CradSoil <<endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " QzSimple2 " << " $matTag " << " $SoilType " << " $Qult-end " << " $z50-end "  << " <TpSoil> " << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " QzSimple2 " << FoundationTag*100+2 <<"  " << SoilType << " " << qend << " " << z50end << " " << TpSoil << " " << CradSoil <<endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " QzSimple2 " << " $matTag " << " $SoilType " << " $Qult-mid " << " $z50-mid "  << " <TpSoil> " << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " QzSimple2 " << FoundationTag*100+3 <<"  " << SoilType << " " << qmid<< " " << z50mid << " " << TpSoil << " " << CradSoil <<endln;

ShallowFoundationOut << endln <<" #fix " << " $midNode " << " " << 1 << " " << 0 << " "<< 0 << endln; 
ShallowFoundationOut <<" fix " << ConNode<< " " << 1 << " " << 0 << " " << 0 << endln; 	
}

if(FootingCondition == 5)    //nonlinear vertical and lateral springs
	{

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " QzSimple2 " << " $matTag " << " $SoilType " << " $Qult-end-extreme " << " $z50-end "  << " <TpSoil> " << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " QzSimple2 " << FoundationTag*100+1 <<"  " << SoilType << " " << qendext << " " << z50endext << " " << TpSoil << " " << CradSoil <<endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " QzSimple2 " << " $matTag " << " $SoilType " << " $Qult-end " << " $z50-end "  << " <TpSoil> " << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " QzSimple2 " << FoundationTag*100+2 <<"  " << SoilType << " " << qend << " " << z50end << " " << TpSoil << " " << CradSoil <<endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " QzSimple2 " << " $matTag " << " $SoilType " << " $Qult-mid " << " $z50-mid "  << " <TpSoil> " << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " QzSimple2 " << FoundationTag*100+3 <<"  " << SoilType << " " << qmid<< " " << z50mid << " " << TpSoil << " " << CradSoil <<endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " PySimple2 " << " $matTag " << " $SoilType " << " $Pp " << " $xp50 "  << " Cd " << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " PySimple2 " << FoundationTag*100+5 <<"  " << SoilType << " " << Pult << " " <<xp50 << " " << Cd << " " << CradSoil <<endln;

ShallowFoundationOut <<	endln << " #uniaxialMaterial " << " TzSimple2 " << " $matTag " << " $SoilType " << " $Tult " << " $xt50 "  << " <CradSoil> "<<endln;
ShallowFoundationOut << " uniaxialMaterial " << " TzSimple2 " << FoundationTag*100+6 <<"  " << SoilType << " " << Tult << " " <<xt50 << " " << Cd << " " << CradSoil <<endln;
	}

//--------END OF MATERIAL GENERATION

//--------ELEMENT GENERATION
//Vertical Spring element connectivity 
ShallowFoundationOut << endln <<" #Vertical spring element connectivity"<<endln; 
ShallowFoundationOut << " #element  " << " zeroLength " << " $eleTag " << " $iNode " << " $jNode " << " -mat" << "$matTag " << " -dir " << " $dir " <<endln;
for(i=1; i<=nodeTotal; i++)
{
	if (i > ndivEnd && i <= ndivEnd+ndivMid+1)
	{
ShallowFoundationOut << " element  " << " zeroLength " << FoundationTag*100000+i << "  " <<FoundationTag*100000+i << "  " << FoundationTag*1000+i << " -mat " << FoundationTag*100+3 <<"  " << " -dir " << " 2 " <<endln;
	}
	else if (i ==1 || i == nodeTotal)
	{
ShallowFoundationOut << " element  " << " zeroLength " << FoundationTag*100000+i << "  " <<FoundationTag*100000+i << "  " << FoundationTag*1000+i << " -mat " << FoundationTag*100+1 <<"  " << " -dir " << " 2 " <<endln;
	} 
	else
	{
ShallowFoundationOut << " element  " << " zeroLength " << FoundationTag*100000+i << "  " <<FoundationTag*100000+i << "  " << FoundationTag*1000+i << " -mat " << FoundationTag*100+2 <<"  " << " -dir " << " 2 " <<endln;
	}
}


//Horizontal Spring element connectivity 

if(FootingCondition == 3)    //elastic lateral spring
	{
ShallowFoundationOut << endln <<" #Horizontal spring element connectivity"<<endln; 
ShallowFoundationOut << " #element  " << " zeroLength " << " $eleTag " << " $iNode " << " $jNode " << " -mat" << "$matTag " << " -dir " << " $dir " <<endln;
ShallowFoundationOut << " element  " << " zeroLength " << FoundationTag*100000+nodeTotal+1 << "  " <<FoundationTag*1000+nodeTotal << "  " << FoundationTag*100000+nodeTotal+1 << " -mat " << FoundationTag*100+4 <<"  " << " -dir " << " 1 " <<endln;
}

if(FootingCondition == 5)    //Nonlinear lateral springs
	{
ShallowFoundationOut << endln <<" #Horizontal spring element connectivity"<<endln; 
ShallowFoundationOut << " #element  " << " zeroLength " << " $eleTag " << " $iNode " << " $jNode " << " -mat" << "$matTag " << " -dir " << " $dir " <<endln;
ShallowFoundationOut << " element  " << " zeroLength " << FoundationTag*100000+nodeTotal+1 << "  " <<FoundationTag*1000+nodeTotal << "  " << FoundationTag*100000+nodeTotal+1 << " -mat " << FoundationTag*100+5 <<"  " << " -dir " << " 1 " <<endln;
ShallowFoundationOut << " element  " << " zeroLength " << FoundationTag*100000+nodeTotal+2 << "  " <<FoundationTag*1000+nodeTotal << "  " << FoundationTag*100000+nodeTotal+2 << " -mat " << FoundationTag*100+6 <<"  " << " -dir " << " 1 " <<endln;
}

// Foundation element connectivity
double transfTag;
transfTag =10*FoundationTag;  //geometric transfer tag (for foundation)
ShallowFoundationOut << endln <<" # geomTransf Linear $transfTag <-jntOffset $dXi $dYi $dXj $dYj>"<<endln;
ShallowFoundationOut <<" geomTransf Linear  "<< transfTag <<endln;

   

ShallowFoundationOut << endln <<" #foundation element connectivity"<<endln; 
ShallowFoundationOut << " #element  " << " elasticBeamColumn " << " $eleTag " << " $iNode " << " $jNode " << " $A " << " $E " << " $Iz " << " $transfTag " <<endln;

for(i=1; i<nodeTotal; i++){
//writing Foundation connectivity
	if (i > ndivEnd && i <= ndivEnd+ndivMid)
	{
ShallowFoundationOut << " element " <<"elasticBeamColumn " << FoundationTag*1000+i << " " <<FoundationTag*1000+i << "  " << FoundationTag*1000+i+1 << " " << Aelefoot << " " << Efoot << " " <<  Ielefoot << " " << transfTag <<endln;
	}
	else
	{
ShallowFoundationOut << " element " <<"elasticBeamColumn " << FoundationTag*1000+i << " " <<FoundationTag*1000+i << "  " << FoundationTag*1000+i+1 << " " << Aelefoot << " " << Efoot << " " <<  Ielefoot << " " << transfTag <<endln;
	}
}

//writing Spring Fixity
ShallowFoundationOut << endln <<" #fixity "<<endln; 
for(i=1; i<=nodeTotal; i++)
{
ShallowFoundationOut << " fix  " <<FoundationTag*100000+i << " 1 1 1"<<endln;
}

if(FootingCondition == 3) {
ShallowFoundationOut << " fix  " <<FoundationTag*100000+nodeTotal+1 << " 1 1 1"<<endln;
}

if(FootingCondition == 5) {
ShallowFoundationOut << " fix  " <<FoundationTag*100000+nodeTotal+1 << " 1 1 1"<<endln;
ShallowFoundationOut << " fix  " <<FoundationTag*100000+nodeTotal+2 << " 1 1 1"<<endln;
}

//--------END OF ELEMENT GENERATION

//Setting required nodes and elements 
ShallowFoundationOut <<	endln << " set endFootNodeL_"<<FoundationTag<< "   "<<FoundationTag*1000+1<<endln;
ShallowFoundationOut << " set endFootNodeR_"<<FoundationTag<< "   "<<FoundationTag*1000+nodeTotal <<endln;
ShallowFoundationOut << " set endSprEleL_"<<FoundationTag<< "   "<<FoundationTag*100000+1<<endln;
ShallowFoundationOut << " set endSprEleR_"<<FoundationTag<< "   "<<FoundationTag*100000+nodeTotal<<endln;
ShallowFoundationOut << " set midSprEle_"<<FoundationTag<< "   "<<FoundationTag*100000+midNode<<endln;
	ShallowFoundationOut.close();

	return;
}



