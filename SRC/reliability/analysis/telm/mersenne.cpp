/************************** MERSENNE.CPP ******************** AgF 2001-10-18 *
*  Random Number generator 'Mersenne Twister'                                *
*                                                                            *
*  This random number generator is described in the article by               *
*  M. Matsumoto & T. Nishimura, in:                                          *
*  ACM Transactions on Modeling and Computer Simulation,                     *
*  vol. 8, no. 1, 1998, pp. 3-30.                                            *
*  Details on the initialization scheme can be found at                      *
*  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html                  *
*                                                                            *
*  Experts consider this an excellent random number generator.               *
*                                                                            *
*  © 2001 - 2004 A. Fog.                                                     *
*  GNU General Public License www.gnu.org/copyleft/gpl.html                  *
*****************************************************************************/

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/mersenne.cpp,v $


#include <randomc.h>

TRandomMersenne::TRandomMersenne(uint32 seed)
:UniformGenerator()
{       // constructor
    RandomInit(seed);
}
void TRandomMersenne::RandomInit(uint32 seed)
 {
  // re-seed generator
  mt[0]= seed;
  for (mti=1; mti < MERS_N; mti++) {
    mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);}

  // detect computer architecture
  union {double f; uint32 i[2];} convert;
  convert.f = 1.0;
  // Note: Old versions of the Gnu g++ compiler may make an error here,
  // compile with the option  -fenum-int-equiv  to fix the problem
  if (convert.i[1] == 0x3FF00000) Architecture = LITTLEENDIAN;
  else if (convert.i[0] == 0x3FF00000) Architecture = BIGENDIAN;
  else Architecture = NONIEEE;}

  
void TRandomMersenne::RandomInitByArray(uint32 seeds[], int length) {
  // seed by more than 32 bits
  int i, j, k;
  RandomInit(19650218UL);
  if (length <= 0) return;
  i = 1;  j = 0;
  k = (MERS_N > length ? MERS_N : length);
  for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + seeds[j] + j;
    i++; j++;
    if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
    if (j >= length) j=0;}
  for (k = MERS_N-1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
    if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
  mt[0] = 0x80000000UL;} // MSB is 1; assuring non-zero initial array

  
uint32 TRandomMersenne::BRandom() {
  // generate 32 random bits
  uint32 y;

  if (mti >= MERS_N) {
    // generate MERS_N words at one time
    const uint32 LOWER_MASK = (1LU << MERS_R) - 1; // lower MERS_R bits
    const uint32 UPPER_MASK = -1L  << MERS_R;      // upper (32 - MERS_R) bits
    static const uint32 mag01[2] = {0, MERS_A};
    
    int kk;
    for (kk=0; kk < MERS_N-MERS_M; kk++) {    
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

    for (; kk < MERS_N-1; kk++) {    
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

    y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
    mti = 0;}

  y = mt[mti++];

  // Tempering (May be omitted):
  y ^=  y >> MERS_U;
  y ^= (y << MERS_S) & MERS_B;
  y ^= (y << MERS_T) & MERS_C;
  y ^=  y >> MERS_L;
  return y;}

  
double TRandomMersenne::Random() {
  // output random float number in the interval 0 <= x < 1
  union {double f; uint32 i[2];} convert;
  uint32 r = BRandom(); // get 32 random bits
  // The fastest way to convert random bits to floating point is as follows:
  // Set the binary exponent of a floating point number to 1+bias and set
  // the mantissa to random bits. This will give a random number in the 
  // interval [1,2). Then subtract 1.0 to get a random number in the interval
  // [0,1). This procedure requires that we know how floating point numbers
  // are stored. The storing method is tested in function RandomInit and saved 
  // in the variable Architecture. The following switch statement can be
  // omitted if the architecture is known. (A PC running Windows or Linux uses
  // LITTLEENDIAN architecture):
  switch (Architecture) {
  case LITTLEENDIAN:
    convert.i[0] =  r << 20;
    convert.i[1] = (r >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case BIGENDIAN:
    convert.i[1] =  r << 20;
    convert.i[0] = (r >> 12) | 0x3FF00000;
    return convert.f - 1.0;
  case NONIEEE: default:
  ;} 
  // This somewhat slower method works for all architectures, including 
  // non-IEEE floating point representation:
  return (double)r * (1./((double)(uint32)(-1L)+1.));}

  
int TRandomMersenne::IRandom(int min, int max) {
  // output random integer in the interval min <= x <= max
  int r; 
  r = int((max - min + 1) * Random()) + min; // multiply interval with random and truncate
  if (r > max) r = max;
  if (max < min) return 0x80000000;
  return r;}

