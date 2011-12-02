#ifndef MATRIX_H 
#define MATRIX_H

#include <math.h>
#include <iostream.h>

class MATRIX;
// Vectors need to know about matrices.

class VECTOR
// Implements 4D vectors
{
public:
  VECTOR(float a1=0, float a2=0, float a3=0, float a4=0)
    {
      v[0] = a1;
      v[1] = a2;
      v[2] = a3;
      v[3] = a4;
    }
  float& operator[](int x) 
    {
      return(v[x]);
    }
  int MakeNonHomogeneous();	// Divide all coordinates by w
  float Length()		
    // Returns the 3D length of the vector
    {
      return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
  int Normalize()
    // Treats this vector as if it were 3D and makes it unit
    // length.
    {
      float l = this->Length();
      if (l == 0)
	{
	  cerr << "ERROR: Trying to normalize a zero-vector" << endl;
	  return 0;
	}
      v[0] = v[0] / l;
      v[1] = v[1] / l;
      v[2] = v[2] / l;
      return 1;
    }
  float Dot(VECTOR &other)
   // 3D dot product
    { 
      return (v[0]*other.v[0] + v[1]*other.v[1] +v[2]*other.v[2]); 
    }

  void TimesMat(MATRIX &Mat);
      
  VECTOR operator%(VECTOR &other) 
    // 3D cross product of two vectors.
    {
      float x, y, z;
      x = (v[1] * other.v[2])-(v[2] * other.v[1]);
      y = (v[2] * other.v[0])-(v[0] * other.v[2]);
      z = (v[0] * other.v[1])-(v[1] * other.v[0]);
      return (*(new VECTOR(x,y,z)));
    }
  VECTOR operator*(MATRIX& _M); 
  // VECTOR * MATRIX
      
private:
    float v[4];
};

class MATRIX
{
public:
  MATRIX()
    // Constructor creates an identity matrix
    {
      int row;
      for (row = 0; row < 4; row++)
	{
	  for (int col = 0; col < 4; col++)
	    {
	      m[row][col] = 0;
	    }
	}
      for (row = 0; row < 4; row++)
	{
	  m[row][row] = 1;
	}
    }
  void Set(float a1, float a2, float a3, float a4,
	   float b1, float b2, float b3, float b4,
	   float c1, float c2, float c3, float c4,
	   float d1, float d2, float d3, float d4);
  // Sets the matrix to the given arguments

  MATRIX operator*(MATRIX &mat); // Matrix multiplications
  float m[4][4];

  MATRIX Transpose(void);

};

// Debugging output

extern ostream &operator<<(ostream &, VECTOR &);      
extern ostream &operator<<(ostream &, MATRIX &);      

#endif      



