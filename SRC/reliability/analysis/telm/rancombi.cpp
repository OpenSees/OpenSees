/************************ RANCOMBI.CPP ****************** AgF 2001-10-18 *
*                                                                        *
*  This file defines a template class for combining two different        *
*  random number generators. A combination of two random number          *
*  generators is generally better than any of the two alone.             *
*  The two generators should preferably be of very different design.     *
*                                                                        *
*  Instructions:                                                         *
*  To make a combined random number generator, insert the class names    *
*  of any two random number generators, as shown in the example below.   *
*                                                                        *
*************************************************************************/

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/rancombi.cpp,v $


#include <randomc.h>

// This template class combines two different random number generators
// for improved randomness. R1 and R2 are any two different random number
// generator classes.
template <class RG1, class RG2>
class TRandomCombined : private RG1, private RG2 {
  public:
  TRandomCombined(int32 seed = 19) : RG1(seed), RG2(seed+1) {};

  void RandomInit(int32 seed) {        // re-seed
    RG1::RandomInit(seed);
    RG2::RandomInit(seed+1);}

  double Random() {
    long double r = RG1::Random() + RG2::Random();
    if (r >= 1.) r -= 1.;
    return r;}
    
  long IRandom(long min, long max){       // output random integer
    // get integer random number in desired interval
    int iinterval = max - min + 1;
    if (iinterval <= 0) return -1; // error
    int i = int(iinterval * Random()); // truncate
    if (i >= iinterval) i = iinterval-1;
    return min + i;}};

  
//////////////////////////////////////////////////////////////////////////
/* Example showing how to use the combined random number generator:
#include <stdio.h>
#include <time.h>
#include "randomc.h"
#include "mersenne.cpp"
#include "ranrotw.cpp"
#include "rancombi.cpp"

int main() {
  // Make an object of the template class. The names inside <> define the
  // class names of the two random number generators to combine.
  // Use time as seed.
  TRandomCombined<TRanrotWGenerator,TRandomMersenne> RG(time(0));

  for (int i=0; i<20; i++) {
    // generate 20 random floating point numbers and 20 random integers
    printf("\n%14.10f   %2i",  RG.Random(),  RG.IRandom(0,99));}

  return 0;}
*/
