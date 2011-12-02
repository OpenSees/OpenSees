/**
 * A very simple unit testing class.
 *
 * @author David M. Doolin.
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


#include <stdio.h>

#include "unittest.h"


void 
UnitTest::print_header(void * stream, const char * tag) {

  char header[] = {"\n\n=================  %s  =================\n"};

  fprintf((FILE*)stream,header,tag);

}


void
UnitTest::register_test_functions(TestFunc * testfunc) {

  this->testfunc = testfunc;
}


/**
 * @todo Add an output stream here such that
 * stdout/cout can go to an arbitrary
 * file, socket, whatever.
 */
bool
UnitTest::test() {

    int i = 0;
    bool passed = true;

    if (testfunc == NULL) {
      return false;
    }

    while (testfunc[i].test != NULL) {

        print_header(stdout,testfunc[i].testname);

	if (testfunc[i].test()) {
	    fprintf(stdout, "Passed test_%s.\n", testfunc[i].testname);
	} else {
	    fprintf(stdout, "Failed test_%s.\n", testfunc[i].testname);
	    passed = false;
	}
	i++;
    }

    return passed;
}


