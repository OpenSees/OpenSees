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

class Patch {
public:

  // inquiring functions

  virtual int getMaterialID() const = 0;
  virtual int getNumCells() const   = 0;
  virtual Cell** getCells() const   = 0;

protected:
private:
};

#endif
