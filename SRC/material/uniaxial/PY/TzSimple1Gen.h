//$Revision: 1.4 $
//$Date: 2004-06-30 00:27:40 $
//$Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PY/TzSimple1Gen.h,v $

#include <iostream>
#include <fstream>
#include <cmath>
#include <ID.h> 

class TzSimple1Gen
{

  // Variables used for reading input files:
  int NumNodes, NumTzEle, NumPileEle, NPile, NumLayer, NumMtLoadSp, NumLoad, NumSp, NumMt, NumMat;
  double p, zground, TULT, Z50, ru, ca, depth, stress, delta, b, Sa;
  int *NodeNum;								// Arrays for Nodes File
  double *Nodey, *Nodex;
  int *TzEleNum, *TzNode1, *TzNode2, *TzMat, *TzDir;	// Arrays for Py Elements File
  int *PileEleNum, *PileNode1, *PileNode2;			// Arrays for Pile Elements File
  int *tzType;
  double *gamma_t, *gamma_b, *z_t, *z_b, *p_t, *p_b, *c_t, *c_b, *ca_t, *ca_b, *delta_t, *delta_b,
    *zLoad_t, *zLoad_b, *load_val_t, *load_val_b, *zSp_t, *zSp_b, *sp_val_t,
    *sp_val_b, *zMt_t, *zMt_b, *mt_val_t, *mt_val_b, tribcoord[2], *Sa_b, *Sa_t, *ru_t, *ru_b,
    *tult_t, *tult_b, *z50_t, *z50_b;
  char **MatType, *PatternInfo;
  
  
  // Member functions for reading input files:
  void GetNodes(const char *file);
  void GetTzElements(const char *file);
  void GetPileElements(const char *file);
  void GetSoilProperties(const char *file);
  double GetTult(const char *type);
  double GetZ50(const char *type);
  double GetMt(double *vx, double *vy, double x, int length);
  void GetTributaryCoordsTz(int nodenum1);
  void GetTributaryCoordsPile(int nodenum1);
  int NumRows(const char *file, const char *begin);
  
  // Member functions for generating output:
  void GetTzSimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5);
  void GetPattern(const char *file6);
  
  // Member functions for calculating tult:
  double GetVStress(double z);
  double linterp(double x1, double x2, double y1, double y2, double x3);
  
 public:
  
  void WriteTzSimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5);
  void WriteTzSimple1(const char *file1, const char *file2, const char *file3, const char *file4, const char *file5, const char *file6);
  
  TzSimple1Gen();
  ~TzSimple1Gen();
};
