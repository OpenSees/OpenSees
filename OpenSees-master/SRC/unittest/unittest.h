/**
 * A very simple unit testing class.
 *
 * @todo have each test() forked into a different
 *  process so that seg faults can be caught and the
 *  tests continued past the offending (and consequently
 *  marked) code.
 * 
 * Copyright (C) 2002, 2003 David M. Doolin
 *
 * This file is part of Geotechnica.
 *
 * Geotechnica is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * Geotechnica is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
 * for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with Geotechnica; see the file COPYING.  If not, write to the Free
 * Software Foundation, 59 Temple Place - Suite 330, Boston,
 * MA  02111-1307, USA.
 */


#ifndef __GEO_UNITTEST_H__
#define __GEO_UNITTEST_H__


/** Each unit test builds a table that can be 
 * traversed to exercise each function in the 
 * component.
 */
typedef struct _testfunc {

  bool       (*test)(void);
  const char * testname;
} TestFunc;



class UnitTest {

 public:

  void register_test_functions (TestFunc * testfunc);

 /**  Perform all of the unit testing specified in a 
  *  table.
  *
  *  @param TestFunc * an array of test functions that 
  *  are invoked sequentially from the calling function.
  *
  *  @return int TRUE if everything passes the unit test,
  *  false if any function fails its unit test.
  */
  bool test                    ();

  void print_header            (void * stream, 
                                const char * tag);

 private:

  TestFunc * testfunc;

};




/** 
 * @todo This probably needs to go into a utility
 * header for just dealing with various types of
 * valarrays.
 */
template <class T>
bool valarray_equals(std::valarray<T> & v1, std::valarray<T> & v2) {

  if (v1.size() != v2.size()) {
    std::cout << "Error: mismatched valarray sizes.\n";
    return false;
  }

  for (unsigned int i=0; i<v1.size(); i++) {
    if (v1[i] != v2[i]) {
      return false;
    }
  }

  return true;
}

template <class T>
void valarray_print(std::valarray<T> & v, const char * label) {

   for (unsigned int i=0; i<v.size(); i++) {
     std::cout << label << "[" << i << "]: " << v[i] << "\n";
   }
}

/**
 * @brief Test whether the arrays are equals within
 *  tolerance for the first <size> elements.
 *
 * @param double * a array of doubles.
 * @param double * b array of doubles.
 * @param int size the number of elements to traverse
 * @param double tol tolerance within which the arrays
 *  are assumed equal.
 *
 * @return bool true if equal within tolerance,
 *  false otherwise.
 *
 * @warning No error checking in this function.
 *  The user is responsible for ensuring that
 *  the allocated sizes are within the size
 *  parameter.
 */
bool   double_array_equals (const double * a,
			    const double * b,
			    const int size,
			    const double tol);


/** 
 * Individual test functions can be prototyped 
 * here for inclusion into other test program drivers.
 */


#endif  /* __GEO_UNITTEST_H__ */
