//$Revision: 1.4 $
//$Date: 2004-06-30 00:27:40 $
//$Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PY/PySimple1Gen.h,v $

#include <fstream>
#include <cmath>
#include <iostream>
#include <ID.h>

class PySimple1Gen
{	
	char **pTest;
	// Variables used for reading input files:
	int NumNodes, NumPyEle, NumPileEle, NumLayer, NumMpLoadSp, NumLoad, NumSp, NumMp, NumMat;
	double pult, y50, b, maxz, minz, depth, cu, e50, stress, phi, sr, PULT, Y50, ru;
	int *NodeNum;								// Arrays for Nodes File
	double *Nodey, *Nodex;
	int *PyEleNum, *PyNode1, *PyNode2, *PyMat, *PyDir;	// Arrays for Py Elements File
	int *PileEleNum, *PileNode1, *PileNode2;			// Arrays for Pile Elements File
	int *pyType, stype;
	double *gamma_t, *gamma_b, *z_t, *z_b, *b_t, *b_b, *Cd_t, *Cd_b, *c_t, *c_b, // Arrays for Soil Properties File
		*cu_t, *cu_b, *e50_t, *e50_b, *phi_t, *phi_b, *Sr_t, *Sr_b, *pult_t, *pult_b,
		*y50_t, *y50_b, *zLoad_t, *zLoad_b, *load_val_t, *load_val_b, *zSp_t, *zSp_b, *sp_val_t,
		*sp_val_b, *zMp_t, *zMp_b, *mp_val_t, *mp_val_b, *ru_t, *ru_b, tribcoord[2];
	char **MatType, *PatternInfo;

	// Member functions for reading input files:
	void GetNodes(const char *file);
	void GetPyElements(const char *file);
	void GetPileElements(const char *file);
	void GetSoilProperties(const char *file);
	int NumRows(const char *file, const char *begin);


	// Member functions for generating output:
	void GetPySimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5);
	void GetPattern(const char *file6);

	// Member functions for calculating pult:
	double GetPult(const char *type);
	double GetY50(const char *type);
	double GetVStress(double z);
	double linterp(double x1, double x2, double y1, double y2, double x3);
	double GetMp(double *vx, double *vy, double x, int length);
	void GetTributaryCoordsPy(int nodenum1);
	void GetTributaryCoordsPile(int nodenum1);
	
public:

	// Public member functions accessed from TclModelBuilder.cpp
	void WritePySimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5);
	void WritePySimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5, const char *file6);

	PySimple1Gen();
	~PySimple1Gen();
};
