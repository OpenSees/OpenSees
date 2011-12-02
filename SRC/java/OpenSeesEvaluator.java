/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2007-04-06 03:43:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/java/OpenSeesEvaluator.java,v $
// 
// Written: fmk 
// Created: 07/04
//

import java.io.*;

public class OpenSeesEvaluator 
{

    // native functions written in c++
    public native int     openSeesInit();
    public native String  openSeesEval(String expression, int error);
    public native int     openSeesQuit();

    // main method .. 
    //   1. calls OpenSeesInit() to initialize OpenSees interpreter.
    //   2. loops until user invokes quit, evaluating OpenSees expressions
    //   3. calls OpenSeesQuit() to cleanup.
    public static void main (String[] args) throws IOException {

	// load the shared object library
	System.loadLibrary("OpenSeesEvaluator");
	//System.load("./OpenSeesEvaluator.so");

	// create an OpenSeesEvaluator
        OpenSeesEvaluator interpreter = new OpenSeesEvaluator();

	// call openSeesInit() to initialize OpenSees interpreter.
	interpreter.openSeesInit();

	// loop until user invokes quit, 
	//    read user i/p expression & evaluate using OpenSees interpreter
	//    by calling openSeesEval(i/p expression)
	BufferedReader stdin = new BufferedReader 
	    (new InputStreamReader(System.in));
	
	String expression = "";
	String result ="";
	int done = 0;
	
	while (done == 0) {
	    System.out.print ("OpenSees> ");
	    expression = stdin.readLine();
	    result =  interpreter.openSeesEval(expression, done);
	    System.out.println (result);
	    if (expression.equals("exit"))
		done = 1;
	} 

	// call openSeesQuit() to close OpenSees interpreter.
	interpreter.openSeesQuit();	
    }  
} 
