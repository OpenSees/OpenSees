/************************* RANROTW.CPP ****************** AgF 1999-03-03 *
*  Random Number generator 'RANROT' type W                               *
*  This version is used when a resolution higher that 32 bits is desired.*
*                                                                        *
*  This is a lagged-Fibonacci type of random number generator with       *
*  rotation of bits.  The algorithm is:                                  *
*  Z[n] = (Y[n-j] + (Y[n-k] rotl r1)) modulo 2^(b/2)                     *
*  Y[n] = (Z[n-j] + (Z[n-k] rotl r2)) modulo 2^(b/2)                     *
*  X[n] = Y[n] + Z[n]*2^(b/2)                                            *
*                                                                        *
*  The last k values of Y and Z are stored in a circular buffer named    *
*  randbuffer.                                                           *
*  The code includes a self-test facility which will detect any          *
*  repetition of previous states.                                        *
*  The function uses a fast method for conversion to floating point.     *
*  This method relies on floating point numbers being stored in the      *
*  standard 64-bit IEEE format or the 80-bit long double format.         *
*                                                                        *
*  The theory of the RANROT type of generators and the reason for the    *
*  self-test are described at www.agner.org/random/theory                *
*                                                                        *
* ©2002 A. Fog. GNU General Public License www.gnu.org/copyleft/gpl.html *
*************************************************************************/

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/ranrotw.cpp,v $


#include "randomc.h"
#include <string.h> // some compilers require <mem.h> instead
#include <stdio.h>
#include <stdlib.h>

// If your system doesn't have a rotate function for 32 bits integers,
// then use the definition below. If your system has the _lrotl function 
// then remove this.

// constructor:
TRanrotWGenerator::TRanrotWGenerator(uint32 seed)
:UniformGenerator()
 {
  RandomInit(seed);
  // detect computer architecture
  randbits[2] = 0; randp1 = 1.0;
  if (randbits[2] == 0x3FFF) Architecture = EXTENDEDPRECISIONLITTLEENDIAN;
  else if (randbits[1] == 0x3FF00000) Architecture = LITTLEENDIAN;
  else if (randbits[0] == 0x3FF00000) Architecture = BIGENDIAN;
  else Architecture = NONIEEE;}


uint32 TRanrotWGenerator::BRandom() {
  // generate next random number
  uint32 y, z;
  // generate next number
  z = _lrotl(randbuffer[p1][0], R1) + randbuffer[p2][0];
  y = _lrotl(randbuffer[p1][1], R2) + randbuffer[p2][1];
  randbuffer[p1][0] = y; randbuffer[p1][1] = z;
  // rotate list pointers
  if (--p1 < 0) p1 = KK - 1;
  if (--p2 < 0) p2 = KK - 1;
  // perform self-test
  if (randbuffer[p1][0] == randbufcopy[0][0] &&
    memcmp(randbuffer, randbufcopy[KK-p1], 2*KK*sizeof(uint32)) == 0) {
      // self-test failed
      if ((p2 + KK - p1) % KK != JJ) {
        // note: the way of printing error messages depends on system
        // In Windows you may use FatalAppExit
        printf("Random number generator not initialized");}
      else {
        printf("Random number generator returned to initial state");}
      exit(1);}
  randbits[0] = y;
  randbits[1] = z;
  return y;}


double TRanrotWGenerator::Random() {
  // returns a random number between 0 and 1.
  uint32 z = BRandom();  // generate 64 random bits
  switch (Architecture) {
  case EXTENDEDPRECISIONLITTLEENDIAN:
    // 80 bits floats = 63 bits resolution
    randbits[1] = z | 0x80000000;  break;
  case LITTLEENDIAN:
    // 64 bits floats = 52 bits resolution
    randbits[1] = (z & 0x000FFFFF) | 0x3FF00000;  break;
  case BIGENDIAN:
    // 64 bits floats = 52 bits resolution
    randbits[0] = (randbits[0] & 0x000FFFFF) | 0x3FF00000;  break;
  case NONIEEE: default:
    // not a recognized floating point format. 32 bits resolution
    return (double)z * (1./((double)(uint32)(-1L)+1.));}
  return (double)(randp1 - 1.0);}


int TRanrotWGenerator::IRandom(int min, int max) {
  // get integer random number in desired interval
  int iinterval = max - min + 1;
  if (iinterval <= 0) return 0x80000000;  // error
  int i = int(iinterval * Random());      // truncate
  if (i >= iinterval) i = iinterval-1;
  return min + i;}


void TRanrotWGenerator::RandomInit (uint32 seed) {
  // this function initializes the random number generator.
  int i, j;

  // make random numbers and put them into the buffer
  for (i=0; i<KK; i++) {
    for (j=0; j<2; j++) {
      seed = seed * 2891336453UL + 1;
      randbuffer[i][j] = seed;}}
  // set exponent of randp1
  randbits[2] = 0; randp1 = 1.0;
  // initialize pointers to circular buffer
  p1 = 0;  p2 = JJ;
  // store state for self-test
  memcpy (randbufcopy, randbuffer, 2*KK*sizeof(uint32));
  memcpy (randbufcopy[KK], randbuffer, 2*KK*sizeof(uint32));
  // randomize some more
  for (i=0; i<31; i++) BRandom();
}

