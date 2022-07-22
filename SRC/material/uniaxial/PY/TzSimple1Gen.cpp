/*
Copyright @ 2002 The Regents of the University of California. All Rights Reserved.

The Regents grants permission, without fee and without a written license agreement, for (a) use, reproduction, modification, and distribution of this software and its documentation by educational, research, and non-profit entities for noncommercial purposes only; and (b) use, reproduction and modification of this software by other entities for internal purposes only. The above copyright notice, this paragraph and the following three paragraphs must appear in all copies and modifications of the software and/or documentation.

Permission to incorporate this software into commercial products may be obtained by contacting the University of California. [ Office of Technology Licensing, 2150 Shattuck Avenue #150 Berkeley, CA 94720-1620, (510) 643-7201]

This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/


////////////////////////////////////////////////////////////////////////
//  This file contains the constructor, destructor, and member		  //  
//  functions for the TzSimple1Gen class.  The purpose of the		  //
//  class is to create TzSimple1 materials associated with			  //
//  pre-defined	zeroLength elements, beam column elements, and		  //
//	nodes.															  //
//																	  //
//  Written by: Scott Brandenberg									  //
//              Graduate Student, UC Davis							  //
//              December 2, 2003									  //	
//				Now Associate Professor, UCLA (sjbrandenberg@ucla.edu //
////////////////////////////////////////////////////////////////////////

//$Revision: 1.6 $
//$Date: 2008/07/03 18:41:05 $
//$Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PY/TzSimple1Gen.cpp,v $

#include "TzSimple1Gen.h"
#include <stdlib.h>
#include <string.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Default Constructor
TzSimple1Gen::TzSimple1Gen()
{
	NumNodes = 0;
	NumTzEle = 0;
	NumPileEle = 0;
	NPile = 0;
	NumLayer = 0;
	NumMtLoadSp = 0;
	NumLoad = 0;
	NumSp = 0;
	NumMt = 0;
	NumMat = 0;
	p = 0.0;
	zground = 0.0;
	TULT = 0.0;
	Z50 = 0.0;
	ru = 0.0;
	ca = 0.0;
	depth = 0.0;
	stress = 0.0;
	delta = 0.0;
	b = 0.0;
	Sa = 0.0;
}

/////////////////////////////////////////////////////////////////////////
// Destructor deletes dynamically allocated arrays
TzSimple1Gen::~TzSimple1Gen()
{
	delete[] Nodex;
	delete[] Nodey;
	delete[] NodeNum;
	delete[] TzEleNum;
	delete[] TzNode1;
	delete[] TzNode2;
	delete[] TzMat;
	delete[] TzDir;
	delete[] PileEleNum;
	delete[] PileNode1;
	delete[] PileNode2;
	delete[] gamma_t;
	delete[] gamma_b;
	delete[] z_t;
	delete[] z_b;
	delete[] p_t;
	delete[] p_b;
	delete[] c_t;
	delete[] c_b;
	delete[] ca_t;
	delete[] ca_b;
	delete[] delta_t;
	delete[] delta_b;
	delete[] zLoad_t;
	delete[] zLoad_b;
	delete[] load_val_t;
	delete[] load_val_b;
	delete[] zSp_t;
	delete[] zSp_b;
	delete[] sp_val_t;
	delete[] sp_val_b;
	delete[] zMt_t;
	delete[] zMt_b;
	delete[] mt_val_t;
	delete[] mt_val_b;
	delete[] Sa_b;
	delete[] Sa_t;
	delete[] ru_t;
	delete[] ru_b;
	delete[] z50_t;
	delete[] z50_b;
	delete[] tult_t;
	delete[] tult_b;
	for(int i=0;i<NumMat;i++)
		delete[] MatType[i];
	delete[] MatType;
	delete[] tzType;
}

///////////////////////////////////////////////////////////////////////////////////////
// Function to call appropriate subroutines given input from main
void TzSimple1Gen::WriteTzSimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5)
{
	GetTzSimple1(file1, file2, file3, file4, file5);
}

void TzSimple1Gen::WriteTzSimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5, const char *file6)
{
	GetTzSimple1(file1, file2, file3, file4, file5);
	GetPattern(file6);
}

///////////////////////////////////////////////////////////////////////////////////////
// Function to write an output file containing tz materials
// given nodes, pile elements and tz elements
void TzSimple1Gen::GetTzSimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5)
{

	// Define local variables
	int i, j, k, tzelenum, TzIndex, tzmat, stype, NODENUM;
	double z, maxz, c, mt;
	double ztrib1, ztrib2, dzsub, zsub, depthsub, sublength, tult, z50, numtzshared;
	mt = 1.0;
	char *mattype;

	// Initialize output stream
	ofstream TzOut;
	TzOut.open(file5, ios::out);

	// Write headers for output file
	TzOut << "#######################################################################################" << endln;
	TzOut << "##" << endln;
	TzOut << "## This file contains TzSimple1 materials associated with pre-defined nodes, zeroLength" << endln;
	TzOut << "## elements and pile beam column elements.  The file was created using the program" << endln;
	TzOut << "## TzSimple1Gen.cpp written by Scott Brandenberg (sjbrandenberg@ucdavis.edu)" << endln;
	TzOut << "##" << endln;
	TzOut << "#######################################################################################" << endln << endln;
	TzOut << "########################################################################################" << endln;
	TzOut << "## Material Properties for tz Elements" << endln << endln;

	// Call functions to open input files
	
	GetSoilProperties(file1);
	GetNodes(file2);
	GetTzElements(file3);
	GetPileElements(file4);
	
	// Loop over nodes
	for(i=0;i<NumTzEle;i++)
	{
		// Initialize variables to zero.  Note that elements and nodes must be assigned numbers larger than zero
		tzelenum = 0;
		z = 0;
		TzIndex = -1;

		// Find number of tz element that shares the node
		for(j=0;j<NumNodes;j++)
		{
			if(NodeNum[j] == TzNode1[i]) // Only calculate values for TzNode that also attaches to pile node
			{
				for(k=0;k<NumPileEle;k++)
				{
					if(PileNode1[k] == TzNode1[i] || PileNode2[k] == TzNode1[i])
					{
                        tzelenum = TzEleNum[i];
						TzIndex = i;
						tzmat = TzMat[i];
						z = Nodey[j];
						NODENUM = NodeNum[j];
					}
				}
			}
			else if(NodeNum[j] == TzNode2[i])
			{
				for(k=0;k<NumPileEle;k++)
				{
					if(PileNode1[k] == TzNode2[i] || PileNode2[k] == TzNode2[i])
					{
                        tzelenum = TzEleNum[i];
						TzIndex = i;
						tzmat = TzMat[i];
						z = Nodey[j];
						NODENUM = NodeNum[j];
					}
				}
			}
			if(TzIndex == -1)
				continue;
		}

		// Find depth of node
		maxz = z_t[0];
		for (j=0;j<NumMat;j++)
		{
			if(z_t[j] > maxz)
				maxz = z_t[j];
		}

		depth = maxz - z;

		GetTributaryCoordsTz(NODENUM);
		ztrib1 = tribcoord[0];
		ztrib2 = tribcoord[1];

		// make sure coordinates of tributary length lie within soil layer for assigning tributary lengths
		// for the tz elements
		if(ztrib1 > maxz)
			ztrib1 = maxz;
		if(ztrib2 > maxz)
			ztrib2 = maxz;

		// Calculate tz material properties and write to file
		if(TzIndex != -1)
		{
			// subdivide tributary length into 10 sublayers, and integrate pult over tributary length
			dzsub = (ztrib2 - ztrib1)/10.0; // sublayer incremental depth change
			sublength = fabs(dzsub);  // thickness of sublayer
			tult = 0.0;
			for(k=0;k<10;k++)
			{			
				zsub = ztrib1 + dzsub/2.0 + k*dzsub; // z-coordinate at sublayer center
				depthsub = maxz - zsub;
			
				// Find properties at node location
				for(j=0;j<NumMat;j++)
				{
					if(zsub<=z_t[j] && zsub>=z_b[j])
					{
						mattype = MatType[j];
						// linearly interpolate parameters at z
						p = linterp(z_t[j], z_b[j], p_t[j], p_b[j], zsub);
						ca = linterp(z_t[j], z_b[j], ca_t[j], ca_b[j], zsub);
						delta = linterp(z_t[j], z_b[j], delta_t[j], delta_b[j], zsub);
						c = linterp(z_t[j], z_b[j], c_t[j], c_b[j], zsub);
						TULT = linterp(z_t[j], z_b[j], tult_t[j], tult_b[j], zsub);
						Z50 = linterp(z_t[j], z_b[j], z50_t[j], z50_b[j], zsub);
						ru = linterp(z_t[j], z_b[j], ru_t[j], ru_b[j], zsub);
						Sa = linterp(z_t[j], z_b[j], Sa_t[j], Sa_b[j], zsub);
						if(strcmp(mattype,"tz1")==0)
							stype = 1;
						else if(strcmp(mattype,"tz2")==0 || strcmp(mattype,"tz3")==0)
							stype = 2;
						else if(strcmp(mattype,"tz4")==0)
							stype = tzType[j];
						else					
						{
							opserr << "MatType must be tz1, tz2, tz3 or tz4.  " << mattype << " is not supported." << endln;
							exit(0);
						}
						break;
					}
				}
								
				for(j=0;j<NumMt;j++)
				{
					if(zsub<=zMt_t[j] && zsub>=zMt_b[j])
						mt = linterp(zMt_t[j], zMt_b[j], mt_val_t[j], mt_val_b[j], zsub);
					else mt = 1.0;
				}

				// calculate vertical effective stress and integrate over tributary length
				stress = GetVStress(zsub);
				tult = GetTult(mattype)*sublength*mt + tult;
			}

			z50 = GetZ50(mattype);

			// Calculate the number of t-z elements that share nodes with the current t-z element
			numtzshared = 1.0;
			for(j=0;j<NumTzEle;j++)
			{
				if(j!=i)
				{
					if(TzNode1[j] == TzNode1[i] || TzNode1[j] == TzNode2[i])
						numtzshared += 1.0;
				}
			}

			TzOut << "uniaxialMaterial TzSimple1 " << tzmat << " " << stype << " " << tult/numtzshared << " " << z50 << " " << c << endln;
		}
	}

	// Write footer for output file
	TzOut << endln << "## End Material Properties for tz Elements" << endln;
	TzOut << "########################################################################################" << endln;

	TzOut.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Function to get applied constraints
void TzSimple1Gen::GetPattern(const char *file6)
{
	double ztrib1, ztrib2, maxz, minz, dzsub, sublength, zsub, depthsub;
	int node, i, j, k;
	double patternvalue, z, load, sp;
	char patterntype[] = "trash";
	
	// Now open a stream to construct the constraints
	ofstream PatternOut;
	PatternOut.open(file6,ios::out);

	if(!PatternOut)
	{
		opserr << "Error opening " << file6 << " in TzSimple1Gen.cpp.  Must Exit." << endln;
		exit(-1);
	}

	patternvalue = 0.0;
	z = 0.0;

		// Write header for constraint file
	PatternOut << "#######################################################################################" << endln;
	PatternOut << "##" << endln;
	PatternOut << "## This file contains load patterns applied to pile nodes, and/or displacement" << endln;
	PatternOut << "## patterns applied to the free ends of tz elements.  The file was created using" << endln;
	PatternOut << "## TzSimple1Gen.cpp written by Scott Brandenberg (sjbrandenberg@ucdavis.edu)" << endln;
	PatternOut << "##" << endln;
	PatternOut << "#######################################################################################" << endln << endln;
	PatternOut << "#######################################################################################" << endln;
	PatternOut << "## Begin Pattern File" << endln << endln;
	
	// If loads are applied (i.e. PatternType = "load"), then the appropriate loads must be assigned to
	// the pile nodes.  The free ends of the tz elements below the nodal loads (i.e. in the soil
	// that is not spreading) must already be fixed when the mesh is generated in GiD.	
	// Write constraints file for pushover analyses
	
	for(i=0;i<NumNodes;i++)
	{
		z = Nodey[i];
		GetTributaryCoordsPile(NodeNum[i]);
		ztrib1 = tribcoord[0];
		ztrib2 = tribcoord[1];
	
		// Find depth of node
		maxz = z_t[0];  // initialize maxz to some value in the domain
		minz = z_b[0];
		for (j=0;j<NumMat;j++)
		{
			if(z_t[j] > maxz)
				maxz = z_t[j];
			if(z_b[j] < minz)
				minz = z_b[j];
		}

			// subdivide tributary length into 10 sublayers, and integrate distributed load over tributary length
		dzsub = (ztrib2 - ztrib1)/10.0; // sublayer incremental depth change
		sublength = fabs(dzsub);  // thickness of sublayer
		load = 0.0;
		for(k=0;k<10;k++)
		{	
			zsub = ztrib1 + dzsub/2.0 + k*dzsub; // z-coordinate at sublayer center
			depthsub = maxz - zsub;

			for(j=0;j<NumLoad;j++)
			{				
				if(zsub<=zLoad_t[j] && zsub>=zLoad_b[j])
				{
					load = linterp(zLoad_t[j], zLoad_b[j], load_val_t[j], load_val_b[j], zsub)*sublength + load;
					strcpy(patterntype,"load");
				}
				
			}
		}
		node = -1;
		if(strcmp(patterntype,"load")==0)
		{
			for(j=0;j<NumPileEle;j++)
			{					
				if(NodeNum[i] == PileNode1[j] || NodeNum[i] == PileNode2[j])
				{
					node = NodeNum[i];
				}
			}
			if(node!=-1)
			PatternOut << "load " << node << " 0.0 " << load << " 0.0" << endln;
		}
	
		for(j=0;j<NumSp;j++)
		{
			if(z<=zSp_t[j] && z>=zSp_b[j])
			{
				sp = linterp(zSp_t[j], zSp_b[j], sp_val_t[j], sp_val_b[j], z);
				strcpy(patterntype,"sp");
			}
		}		
	
		node = -1;
		if(strcmp(patterntype,"sp")==0)
		{
			for(k=0;k<NumTzEle;k++)
			{
				if(NodeNum[i] == TzNode1[k] || NodeNum[i] == TzNode2[k])
				{
					node = NodeNum[i];
					// Check if node is free or attached to pile
					for(j=0;j<NumPileEle;j++)
					{
						if(PileNode1[j] == NodeNum[i] || PileNode2[j] == NodeNum[i])
						{
							node = -1;
							break;
						}
					}
				}
			}
				
		// write to file
			if(node != -1)
				PatternOut << "sp " << node << " 2 " << sp << endln;
		}
		
	}

	PatternOut << endln << endln;
	PatternOut << "## End Tz Pattern File" << endln;
	PatternOut << "#######################################################################################" << endln;

	PatternOut.close();
}

/////////////////////////////////////////////////////////////////////////
// Function to get node numbers and coordinates
void TzSimple1Gen::GetNodes(const char *file)
{
	int i = 0;
	char *trash = new char[1000];
	char ch;
	
	ifstream in_file(file, ios::in);
	
	if(!in_file)
	{
		opserr << "File " << file << "does not exist.  Must exit." << endln;
		exit(-1);
	}
	
	NumNodes = NumRows(file,"node");
	NodeNum = new int[NumNodes];
	Nodex = new double[NumNodes];
	Nodey = new double[NumNodes];
	
	while(!in_file.eof())
	{
		if(in_file.peek()=='n')
		{
			in_file.get(trash,5);
			if(strcmp(trash,"node")==0)
			{
				in_file >> NodeNum[i] >> Nodex[i] >> Nodey[i];
				i+=1;
			}
		}
		while(in_file.get(ch))
		{
			if(ch=='\n')
				break;
		}
	}

	delete[] trash;
	in_file.close();
	return;
}

//////////////////////////////////////////////////////////////////////////////////
// Function to get tz element numbers, material numbers, and direction tags
void TzSimple1Gen::GetTzElements(const char *file)
{
	int i = 0;
	char *trash = new char[1000];
	char ch;

	ifstream in_file;
	in_file.open(file, ios::in);

	if(!in_file)
	{
		opserr << "File " << file << "does not exist.  Must exit." << endln;
		exit(-1);
	}

	NumTzEle = NumRows(file,"element");
	TzEleNum = new int[NumTzEle];
	TzNode1 = new int[NumTzEle];
	TzNode2 = new int[NumTzEle];
	TzMat = new int[NumTzEle];
	TzDir = new int[NumTzEle];

	while(!in_file.eof())
	{
		if(in_file.peek()=='e')
		{
			in_file.get(trash,8);
			if(strcmp(trash,"element")==0)
			{
				in_file >> trash >> TzEleNum[i] >> TzNode1[i] >> TzNode2[i] >> trash >> TzMat[i] >> trash >> TzDir[i];
				i+=1;
			}
			continue;
		}
		while(in_file.get(ch))
		{
			if(ch=='\n')
				break;
		}
	}

	delete[] trash;
	in_file.close();
	return;
}

//////////////////////////////////////////////////////////////////////////////////
// Function to get pile element numbers and node numbers
void TzSimple1Gen::GetPileElements(const char *file)
{
	int i = 0;
	char* trash = new char[1000];
	char ch;
	
	ifstream in_file;
	in_file.open(file, ios::in);

	if(!in_file)
	{
		opserr << "File " << file << "does not exist.  Must exit." << endln;
		exit(-1);
	}

	NumPileEle = NumRows(file,"element");
	PileEleNum = new int[NumPileEle];
	PileNode1 = new int[NumPileEle];
	PileNode2 = new int[NumPileEle];
	
	while(!in_file.eof())
	{
		if(in_file.peek()=='e')
		{
			in_file.get(trash,8);
			if(strcmp(trash,"element")==0)
			{
				in_file >> trash >> PileEleNum[i] >> PileNode1[i] >> PileNode2[i];
				i+=1;
			}
			continue;
		}
		while(in_file.get(ch))
		{
			if(ch=='\n')
				break;
		}
	}

	delete[] trash;
	in_file.close();
	return;
}

//////////////////////////////////////////////////////////////////////////////////
// Function to get soil properties
void TzSimple1Gen::GetSoilProperties(const char *file)
{
	int i = 0;
	int I = 0;
	int J = 0;
	int K = 0;
	char OptionalTag[10];
	
	ifstream in1;
	in1.open(file, ios::in);	
	
	if(!in1)
	{
		opserr << "File " << file << "does not exist.  Must exit." << endln;
		exit(0);
	}

	// Define number of rows containing properties to define TzSimple1 materials
	NumMat = NumRows(file, "tz1") + NumRows(file, "tz2") + NumRows(file,"tz3") + NumRows(file,"tz4");  // Number of tz materials defined in file
	NumMt = NumRows(file,"mt"); // Number of t-multiplier terms defined in file
	NumSp = NumRows(file,"sp");		    // number of applied displacements defined in file
	NumLoad = NumRows(file,"load");     // number of applied distributed loads defined in file
	NumMtLoadSp = NumMt + NumSp + NumLoad + NumRows(file,"Pattern") + NumRows(file,"pattern");		// total number of applied patterns defined in file

	// Dynamically allocate memory for arrays containing information for each soil layer.
	// Arguments general to all layers
	MatType = new char*[4];
	for(i=0;i<NumMat;i++)
		MatType[i] = new char[4];
	z_t = new double[NumMat];
	z_b = new double[NumMat];
	gamma_t = new double[NumMat];
	gamma_b = new double[NumMat];
	p_t = new double[NumMat];
	p_b = new double[NumMat];
	c_t = new double[NumMat];
	c_b = new double[NumMat];
	ca_t = new double[NumMat];
	ca_b = new double[NumMat];
	delta_t = new double[NumMat];
	delta_b = new double[NumMat];
	Sa_t = new double[NumMat];
	Sa_b = new double[NumMat];
	ru_t = new double[NumMat];
	ru_b = new double[NumMat];
	tzType = new int[NumMat];
	tult_t = new double[NumMat];
	tult_b = new double[NumMat];
	z50_t = new double[NumMat];
	z50_b = new double[NumMat];

	// Dynamically allocate memory for arrays containing information for p-multipliers
	zMt_t = new double[NumMt];
	zMt_b = new double[NumMt];
	mt_val_t = new double[NumMt];
	mt_val_b = new double[NumMt];
	
	// Dynamically allocate memory for arrays containing information for load pattern
	zLoad_t = new double[NumLoad];
	zLoad_b = new double[NumLoad];
	load_val_t = new double[NumLoad];
	load_val_b = new double[NumLoad];

	// Dynamically allocate memory for arrays containing information for displacement pattern
	zSp_t = new double[NumSp];
	zSp_b = new double[NumSp];
	sp_val_t = new double[NumSp];
	sp_val_b = new double[NumSp];

	for(i=0;i<NumMat;i++)
	{
		// initialize variables to zero, then redefine later
		c_t[i] = 0;
		c_b[i] = 0;
		ca_t[i] = 0;
		ca_b[i] = 0;
		delta_t[i] = 0;
		delta_b[i] = 0;

		// read in arguments that are common to all material types
		in1 >> MatType[i] >> z_t[i] >> z_b[i] >> gamma_t[i] >> gamma_b[i];
	
		// read in arguments that are specific to certain material types
		if(strcmp(MatType[i],"tz1")==0)
		{
			in1 >> p_t[i] >> p_b[i] >> ca_t[i] >> ca_b[i];
			if(in1.peek() != '\n')
				in1 >> c_t[i] >> c_b[i];
		}
		else if(strcmp(MatType[i],"tz2")==0)
		{
			in1 >> p_t[i] >> p_b[i] >> delta_t[i] >> delta_b[i];
			if(in1.peek() != '\n')
				in1 >> c_t[i] >> c_b[i];
		}		
		else if(strcmp(MatType[i],"tz3")==0)
		{
			in1 >> p_t[i] >> p_b[i] >> delta_t[i] >> delta_b[i] >> Sa_t[i] >> Sa_b[i] >> ru_t[i] >> ru_b[i];
			if(in1.peek() != '\n')	
				in1 >> c_t[i] >> c_b[i];
		}
		else if(strcmp(MatType[i],"tz4")==0)
		{
			in1 >> tzType[i] >> tult_t[i] >> tult_b[i] >> z50_t[i] >> z50_b[i];
			if(in1.peek() != '\n')
				in1 >> c_t[i] >> c_b[i];
		}
		else
		{
			opserr << "MatType " << MatType[i] << "Is not supported in TzSimple1Gen.cpp.";
			exit(0);
		}
		// read to next line or next character
		if(in1.peek()=='\n')
			in1.ignore(100000,'\n');
		if(in1.peek()==' ')
			in1.ignore(100000,' ');
	}
	
	// Read in values that define patterns (either loads applied directly to the pile nodes, or free-field
	// displacements applied to the backs of the tz elements).
	// Read in values that define patterns (either loads applied directly to the pile nodes, or free-field
	// displacements applied to the backs of the py elements).
	for(i=0;i<NumMtLoadSp;i++)
	{
		in1 >> OptionalTag;
		if(strcmp(OptionalTag,"load")==0)
		{
			in1 >> zLoad_t[I] >> zLoad_b[I] >> load_val_t[I] >> load_val_b[I];
			I+=1;
		}
		if(strcmp(OptionalTag,"sp")==0)
		{
			in1 >> zSp_t[J] >> zSp_b[J] >> sp_val_t[J] >> sp_val_b[J];
			J+=1;
		}
		if(strcmp(OptionalTag,"mt")==0)
		{
			in1 >> zMt_t[K] >> zMt_b[K] >> mt_val_t[K] >> mt_val_b[K];
			K+=1;
		}
		if(in1.peek()=='\n')
			in1.ignore(100000,'\n');
		if(in1.peek()==' ')
			in1.ignore(100000,' ');
	}

	in1.close();
}

	
///////////////////////////////////////////////////////////////////////////////////////
// Member function to calculate pult
double TzSimple1Gen::GetTult(const char *type)
{
	double tult_0, tult_1, tult_ru;

	// Calculate tult for clay
	if(strcmp(type,"tz1")==0)
	{
		return ca*p;
	}

	// Calculate tult for sand
	else if(strcmp(type,"tz2")==0)
	{
		if(depth == 0)
			return 0.00001;  // TzSimple1 does not support tult = 0;

		double Ko;
		double deg = 3.141592654/180.0;
		Ko = 0.4; // Use 0.4 like LPile
		return Ko*stress*tan(delta*deg)*p;
	}
	else if(strcmp(type,"tz3")==0)
	{
		double Ko;
		double deg = 3.141592654/180.0;
		// convert phi, alpha and beta from degrees to radians
		Ko = 0.4; // Use 0.4 like LPile
		tult_0 = Ko*stress*tan(delta*deg)*p;

		tult_1 = Sa*p*stress;

		tult_ru = linterp(0.0, 1.0, tult_0, tult_1, ru);

		return tult_ru;
	}
	else if(strcmp(type,"tz4")==0)
	{
		return TULT;
	}
	else
	{
		opserr << "TzType " << type << " is not supported in TzSimple1GenPushover::GetTult.  Setting tult = 0.00000001";
		return 0.00000001;
	}

}

///////////////////////////////////////////////////////////////////////////////////////
// Member function to return y50
double TzSimple1Gen::GetZ50(const char *type)
{
	if(strcmp(type,"tz4")==0)
		return Z50;

	// Set z50 such that zult = 0.5% of pile diameter.  Calculate pile diameter from perimeter for a circular section.
	else
        return 0.005*p/3.14159/8.0;
}

///////////////////////////////////////////////////////////////////////////////////////
// Member function to get the number of rows in a file that begin with a certain string
int TzSimple1Gen::NumRows(const char *file, const char *begin)
{
	if(!file)
	{
		opserr << "File " << file << "does not exist.  Must exit." << endln;
		exit(0);
	}

	ifstream in_file;
	in_file.open(file, ios::in);
	int i = 0;
	char *filein = new char[20];
	
	while(!in_file.eof())
	{
		// check for blank lines
		while(in_file.peek()=='\n')
			in_file.getline(filein,1,'\n');
		// Read first character string
		in_file.get(filein,19,' ');
		if(strcmp(filein, begin)==0)
			i = i+1;

		// Read remainder of line
		in_file.ignore(1000,'\n');
	}
	
	delete [] filein;

	in_file.close();
	return i;

}

/////////////////////////////////////////////////////////////////////////////////////////
// Member function to calculate vertical effective stress at a depth given the unit weight and depth arrays already read in.
double TzSimple1Gen::GetVStress(double z)
{
	double stress, maxz, minz, z_top, z_bot, gamma_top, gamma_bot, gamma_z;
	int i;
	stress = 0;
	maxz = z_t[0];
	minz = z_b[0];
	z_top = 0;
	z_bot = 0;
	gamma_top = 0;
	gamma_bot = 0;

	
	// Find maximum and minimum of depth range specified in z_t and z_b
	for (i=0;i<NumMat;i++)
	{
		if(z_t[i] >= maxz)
			maxz = z_t[i];
		if(z_b[i] <= minz)
			minz = z_b[i];
	}

	// Check that z lies within range of z_t and z_b
	if(z > maxz || z < minz)
	{
		opserr << "Depth lies out of range of specified depth vectors in function 'vstress' in PySimple1GenPushover. Setting stress = 0." << endln;
		return 0.0;
	}


	// Extract coordinates of top and bottom of layer
	for(i=0;i<NumMat;i++)
	{
		if(z >= z_b[i] && z <= z_t[i])
		{
			z_top = z_t[i];
			z_bot = z_b[i];
			gamma_top = gamma_t[i];
			gamma_bot = gamma_b[i];
		}
	}
	


	// Linearly interpolate unit weight at z
	gamma_z = linterp(z_top, z_bot, gamma_top, gamma_bot, z);

	// calculate stress
	for (i=0;i<NumMat;i++)
	{
		if(z <= z_b[i])
			stress = stress + 0.5*(gamma_t[i] + gamma_b[i])*(z_t[i] - z_b[i]);
		if(z > z_b[i] && z < z_t[i])
			stress = stress + 0.5*(gamma_t[i] + gamma_z)*(z_t[i] - z);
	}
	
	return stress;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Function to linearly interpolate
double TzSimple1Gen::linterp(double x1, double x2, double y1, double y2, double x3)
{
	return y1 + (x3-x1)*(y2-y1)/(x2-x1);
}

/////////////////////////////////////////////////////////////////////////////////
// Function that returns the coordinates of the ends of the tributary length
// based on t-z element locations.  Tributary length is based on 1/2 of the pile
// length above nodenum1 and 1/2 of the pile length below nodenum1 as long as
// the pile elements above and below nodenum1 both attach to t-z elements at both nodes.
void TzSimple1Gen::GetTributaryCoordsTz(int nodenum1)
{
	
	double coordnodenum1;
	int i, j, k, I, tzeletag;
	I = 0;

	// initialize tribcoord to the coordinate of nodenum1
	for(i=0; i<NumNodes; i++)
	{
		if(nodenum1 == NodeNum[i])
		{
			coordnodenum1 = Nodey[i];
			tribcoord[0] = Nodey[i];
			tribcoord[1] = Nodey[i];
		}
	}
	for(i=0; i<NumPileEle; i++)
	{
		if(PileNode1[i] == nodenum1)
		{
			tzeletag = 0;
			for(j=0; j<NumTzEle; j++)
			{
				if(TzNode1[j] == PileNode1[i] || TzNode2[j] == PileNode1[i])
				{
					for(k=0; k<NumTzEle; k++)
					{
						if(TzNode1[k] == PileNode2[i] || TzNode2[k] == PileNode2[i])
							tzeletag = 1;  // set pyeletag = 1 if PileNode1 is attached to a py element
					}
				}
			}
			if(tzeletag==1)
			{
				for(j=0; j<NumNodes; j++)
				{
					if(PileNode2[i] == NodeNum[j])
					{
						tribcoord[0] = coordnodenum1 + 0.5*(Nodey[j] - coordnodenum1);
					}
				}
			}
		}
		if(PileNode2[i] == nodenum1)
		{
			tzeletag = 0;
			for(j=0;j<NumTzEle;j++)
			{
				if(TzNode1[j] == PileNode2[i] || TzNode2[j] == PileNode2[i])
				{
					for(k=0; k<NumTzEle; k++)
					{
						if(TzNode1[k] == PileNode1[i] || TzNode2[k] == PileNode1[i])
                            tzeletag = 1;  // set pyeletag = 1 if PileNode2 is attached to a py element
					}
				}
			}
			if(tzeletag==1)
			{
				for(j=0; j<NumNodes; j++)
				{
					if(PileNode1[i] == NodeNum[j])
					{
						tribcoord[1] = coordnodenum1 + 0.5*(Nodey[j] - coordnodenum1);
					}
				}
			}
		}
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////
// Function that returns the coordinates of the ends of the tributary length
// based on pile element locations.  Tributary length is based on 1/2 of the pile
// length above nodenum1 and 1/2 of the pile length below nodenum1 even if
// the pile elements above and below nodenum1 do not both attach to p-y elements 
// at both nodes.
void TzSimple1Gen::GetTributaryCoordsPile(int nodenum1)
{
	
	double coordnodenum1;
	int i, j, I;
	I = 0;

	// initialize tribcoord to the coordinate of nodenum1
	for(i=0; i<NumNodes; i++)
	{
		if(nodenum1 == NodeNum[i])
		{
			coordnodenum1 = Nodey[i];
			tribcoord[0] = Nodey[i];
			tribcoord[1] = Nodey[i];
		}
	}
	for(i=0; i<NumPileEle; i++)
	{
		if(PileNode1[i] == nodenum1)
		{
			for(j=0; j<NumNodes; j++)
			{
				if(PileNode2[i] == NodeNum[j])
				{
					tribcoord[0] = coordnodenum1 + 0.5*(Nodey[j] - coordnodenum1);
				}
			}
		}
		if(PileNode2[i] == nodenum1)
		{
			for(j=0; j<NumNodes; j++)
			{
				if(PileNode1[i] == NodeNum[j])
				{
					tribcoord[1] = coordnodenum1 + 0.5*(Nodey[j] - coordnodenum1);
				}
			}
		}
	}

	return;
}
