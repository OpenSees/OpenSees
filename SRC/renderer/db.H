#ifndef DB_H
#define DB_H

#include "gMatrix.H"
#include "container.H"

class MYPOINT 
// MYpoints are members of FACEs. They have a color attribute
// and a position attribute.
{
public:
  MYPOINT(int _ID, float XX, float YY, float ZZ, 
        float _NX = 0, float _NY = 0, float _NZ = 0,
        float _VWX = 0, float _VWY = 0, float _VWZ = 0)
    : p(XX,YY,ZZ,1), nrm(_NX,_NY,_NZ,1), Viewp(_VWX, _VWY, _VWZ, 1)
    {
      copy = 0;
    }
  MYPOINT *MakeCopy()	
    // Makes a copy of the MYpoint the first time. Returns that copy
    // subsequent times.
    {
      if (copy == 0)
	{
	  copy = new MYPOINT(*this);
	}
      return copy;
    }
  void Transform(MATRIX *M);
    
  void ClearCopy()
    // Used to reset MakeCopy
    {
      if (copy != 0)
	delete copy;
      copy = 0;
    }
  float &X()
    {
      return p[0];
    }
  float &Y()
    {
      return p[1];
    }
  float &Z()
    {
      return p[2];
    }

  VECTOR p;			// The position of the point
  VECTOR nrm;                   // The normal field
  VECTOR Viewp;                 // The original coords in the view coord system

  float r,g,b;			// The color of the point

 private:
  MYPOINT *copy;
};


class FACE
// An orderered set of points along with a surface
{
public:
  FACE(){}

  void AddPoint(MYPOINT &_newPt) 
    // Adds a point to the end of the pointlist.
    {
      pointList.AddTail(_newPt);
    }

  VECTOR nrm;
  
  LIST<MYPOINT> pointList;	// List of points in face
};

#define TRANS_SCALE 1
#define TRANS_TRANSLATE 2
#define TRANS_ROTATE 3



class EDGE
{
public:
  EDGE(MYPOINT *P1, MYPOINT *P2 , int slope) {
       Y_Top = P1->Y(); X_Top = P1->X(); R_Top = P1->r; G_Top = P1->g; B_Top = P1->b; Z_Top = P1->Z();
       Y_Bot = P2->Y(); X_Bot = P2->X(); R_Bot = P2->r; G_Bot = P2->g; B_Bot = P2->b; Z_Bot = P2->Z();
       XV_Top = P1->Viewp[0]; YV_Top = P1->Viewp[1]; ZV_Top = P1->Viewp[2];
       XV_Bot = P2->Viewp[1]; YV_Bot = P2->Viewp[1]; ZV_Bot = P2->Viewp[2];
       down = slope;
     if ((Y_Top - Y_Bot) != 0.0) {
        dX_dY = (X_Top - X_Bot) / (Y_Bot - Y_Top);
        dR_dY = (R_Top - R_Bot) / (Y_Bot - Y_Top);
        dG_dY = (G_Top - G_Bot) / (Y_Bot - Y_Top);
        dB_dY = (B_Top - B_Bot) / (Y_Bot - Y_Top);
        dZ_dY = (Z_Top - Z_Bot) / (Y_Bot - Y_Top);
        dXV_dY = (XV_Top - XV_Bot) / (Y_Bot - Y_Top);
        dYV_dY = (YV_Top - YV_Bot) / (Y_Bot - Y_Top);
        dZV_dY = (ZV_Top - ZV_Bot) / (Y_Bot - Y_Top);

     }  /*  else we are just going to ignore anyway so do nothing */
    
     // set the initial current values for flr(Ytop) -1, i.e. for first scanline
     X_Cur = X_Top + (Y_Top -(int)(Y_Top-1)) * dX_dY;
     R_Cur = R_Top + (Y_Top -(int)(Y_Top-1)) * dR_dY;
     G_Cur = G_Top + (Y_Top -(int)(Y_Top-1)) * dG_dY;
     B_Cur = B_Top + (Y_Top -(int)(Y_Top-1)) * dB_dY;
     Z_Cur = Z_Top + (Y_Top -(int)(Y_Top-1)) * dZ_dY;
     XV_Cur = XV_Top + (Y_Top -(int)(Y_Top-1)) * dXV_dY;
     YV_Cur = YV_Top + (Y_Top -(int)(Y_Top-1)) * dYV_dY;
     ZV_Cur = ZV_Top + (Y_Top -(int)(Y_Top-1)) * dZV_dY;
  }

  // to increment the current edge attributes for going down to next scanline
  void INCREMENT(void) {
     X_Cur += dX_dY;
     R_Cur += dR_dY;
     G_Cur += dG_dY;
     B_Cur += dB_dY;
     Z_Cur += dZ_dY;
     XV_Cur += dXV_dY;
     YV_Cur += dYV_dY;
     ZV_Cur += dZV_dY;

  }

  int edge_YTop() { return (int)Y_Top; }

  int edge_YBot() { return (int)Y_Bot; }

  int edge_Xcur() { // return (int)X_Cur; }
	if (X_Top < X_Bot) 
	   if (X_Cur > X_Bot)
		return (int)X_Bot;
	   else
	        return (int)X_Cur;
     	else if (X_Cur < X_Bot)
		return (int)X_Bot;
	     else
		return (int)X_Cur;
  }     

  float edge_Rcur() { return R_Cur; }

  float edge_Gcur() { return G_Cur; }

  float edge_Bcur() { return B_Cur; }

  float edge_Zcur() { return Z_Cur; }
  float edge_XVcur() { return XV_Cur; }
  float edge_YVcur() { return YV_Cur; }
  float edge_ZVcur() { return ZV_Cur; }
  
  float Y_Top, X_Cur;  // needed for sorting in the list
 
  int edge_Slope() { return down; }

  void edge_SlopeSet(int i) { down = i; }

private:
  float Y_Bot, R_Cur, G_Cur, B_Cur, dX_dY, dR_dY, dG_dY, dB_dY;
  float R_Top, R_Bot, G_Top, G_Bot, B_Top, B_Bot, X_Top, X_Bot;
  float Z_Top, Z_Bot, Z_Cur, dZ_dY;
  float XV_Top, YV_Top, ZV_Top, XV_Cur, YV_Cur, ZV_Cur;
  float XV_Bot, YV_Bot, ZV_Bot, dXV_dY, dYV_dY, dZV_dY;

  int down;

};

#endif

