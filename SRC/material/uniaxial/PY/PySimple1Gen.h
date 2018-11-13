//$Revision: 1.4 $
//$Date: 2004-06-30 00:27:40 $
//$Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PY/PySimple1Gen.h,v $

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
