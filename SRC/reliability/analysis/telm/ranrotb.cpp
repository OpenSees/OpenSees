/************************* RANROTB.CPP ****************** AgF 1999-03-03 *
*  Random Number generator 'RANROT' type B                               *
*                                                                        *
*  This is a lagged-Fibonacci type of random number generator with       *
*  rotation of bits.  The algorithm is:                                  *
*  X[n] = ((X[n-j] rotl r1) + (X[n-k] rotl r2)) modulo 2^b               *
*                                                                        *
*  The last k values of X are stored in a circular buffer named          *
*  randbuffer.                                                           *
*  The code includes a self-test facility which will detect any          *
*  repetition of previous states.                                        *
*                                                                        *
*  The theory of the RANROT type of generators and the reason for the    *
*  self-test are described at www.agner.org/random                       *
*                                                                        *
* © 2002 A. Fog. GNU General Public License www.gnu.org/copyleft/gpl.html*
*************************************************************************/

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/ranrotb.cpp,v $

#include "randomc.h"
#include <string.h>  // some compilers require <mem.h> instead
#include <stdio.h>
#include <stdlib.h>

// If your system doesn't have a rotate function for 32 bits integers,
// then use the definition below. If your system has the _lrotl function 
// then remove this.
// uint32 _lrotl (uint32 x, int r) {
//   return (x << r) | (x >> (sizeof(x)*8-r));}


// constructor:
TRanrotBGenerator::TRanrotBGenerator(uint32 seed) {
  RandomInit(seed);
  // detect computer architecture
  union {double f; uint32 i[2];} convert;
  convert.f = 1.0;
  if (convert.i[1] == 0x3FF00000) Architecture = LITTLEENDIAN;
  else if (convert.i[0] == 0x3FF00000) Architecture = BIGENDIAN;
  else Architecture = NONIEEE;}


// returns a random number between 0 and 1:
double TRanrotBGenerator::Random() {
  uint32 x;
  // generate next random number
  x = randbuffer[p1] = _lrotl(randbuffer[p2], R1) + _lrotl(randbuffer[p1], R2);
  // rotate list pointers
  if (--p1 < 0) p1 = KK - 1;
  if (--p2 < 0) p2 = KK - 1;
  // perform self-test
  if (randbuffer[p1] == randbufcopy[0] &&
    memcmp(randbuffer, randbufcopy+KK-p1, KK*sizeof(uint32)) == 0) {
      // self-test failed
      if ((p2 + KK - p1) % KK != JJ) {
        // note: the way of printing error messages depends on system
        // In Windows you may use FatalAppExit
        printf("Random number generator not initialized");}
      else {
        printf("Random number generator returned to initial state");}
      exit(1);}
  // conversion to float:
  union {double f; uint32 i[2];} convert;
  switch (Architecture) {
  case LITTLEENDIAN:
    convert.i[0] =  x << 20;
    convert.i[1] = (x >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case BIGENDIAN:
    convert.i[1] =  x << 20;
    convert.i[0] = (x >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case NONIEEE: default:
  ;} 
  // This somewhat slower method works for all architectures, including 
  // non-IEEE floating point representation:
  return (double)x * (1./((double)(uint32)(-1L)+1.));}


// returns integer random number in desired interval:
int TRanrotBGenerator::IRandom(int min, int max) {
  int iinterval = max - min + 1;
  if (iinterval <= 0) return -0x80000000; // error
  int i = iinterval * Random(); // truncate
  if (i >= iinterval) i = iinterval-1;
  return min + i;}
  

void TRanrotBGenerator::RandomInit (uint32 seed) {
  // this function initializes the random number generator.
  int i;
  uint32 s = seed;

  // make random numbers and put them into the buffer
  for (i=0; i<KK; i++) {
    s = s * 2891336453 + 1;
    randbuffer[i] = s;}

  // check that the right data formats are used by compiler:
  union {
    double randp1;
    uint32 randbits[2];};
  randp1 = 1.5;
  assert(randbits[1]==0x3FF80000); // check that IEEE double precision float format used

  // initialize pointers to circular buffer
  p1 = 0;  p2 = JJ;
  // store state for self-test
  memcpy (randbufcopy, randbuffer, KK*sizeof(uint32));
  memcpy (randbufcopy+KK, randbuffer, KK*sizeof(uint32));
  // randomize some more
  for (i=0; i<9; i++) Random();
}

