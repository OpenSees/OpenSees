/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */
//
// File: Patch.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef Patch_h 
#define Patch_h 

class Cell;
class OPS_Stream;

class Patch
{
  public:
    // edition functions

    virtual void setMaterialID (int materialID) = 0;
    
    // inquiring functions

    virtual int     getMaterialID() const = 0; 
    virtual int     getNumCells() const = 0;
    virtual Cell  **getCells() const = 0;
    virtual Patch  *getCopy() const = 0;

    virtual void Print(OPS_Stream &s, int flag =0) const =0;   
    friend OPS_Stream &operator<<(OPS_Stream &s, const Patch &patch);    

  protected:
    
  private:
};


#endif
