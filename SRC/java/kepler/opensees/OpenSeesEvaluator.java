
/* OpenSeesEvaluator actor for using OpenSees with Ptolemy II.

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
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2007-04-06 03:58:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/java/kepler/opensees/OpenSeesEvaluator.java,v $
                                                                        
// Written: fmk 
// Created: 04/07

/*

 All rights reserved.
 Permission is hereby granted, without written agreement and without
 license or royalty fees, to use, copy, modify, and distribute this
 software and its documentation for any purpose, provided that the above
 copyright notice and the following two paragraphs appear in all copies
 of this software.

 IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
 ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
 SUCH DAMAGE.

 THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
 PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
 CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
 ENHANCEMENTS, OR MODIFICATIONS.

                                        PT_COPYRIGHT_VERSION_2
                                        COPYRIGHTENDKEY
*/

// PACKAGE
package edu.opensees;

// IMPORTS
import ptolemy.actor.IOPort;
import ptolemy.actor.TypedIOPort;
import ptolemy.data.Token;
import ptolemy.data.StringToken;
import ptolemy.data.type.BaseType;
import ptolemy.data.type.Type;
import ptolemy.data.expr.Parameter;
import ptolemy.data.expr.FileParameter;

import ptolemy.actor.TypedAtomicActor;
import ptolemy.kernel.CompositeEntity;
import ptolemy.kernel.util.*;
import java.io.*;

//////////////////////////////////////////////////////////////////////////
//// OpenSees Evaluator
/**
 * OpenSeesEvaluator: 
"This is the implementation of a OpenSeesEvaluator actor using Ptolemy II. This actor has several input ports, has a file parameter and 
has a single output port. the ouput port outputs the result of running OpenSees and passing in the strings from the i/p port, followed
by a source of the file provided (if there is one)."
 *  @author Frank McKenna, UC Berkeley, March 2007
*/

public class OpenSeesEvaluator extends TypedAtomicActor {


   /** Construct a OpenSeesEvaluator source with the given container and name.
     *  @param OpenSeesEvaluator The name of this actor.
     *  @exception IllegalActionException If the entity cannot be contained
     *   by the proposed container.
     *  @exception NameDuplicationException If the container already has an
     *   actor with this name.
     */
    public OpenSeesEvaluator(CompositeEntity container, String name)
            throws NameDuplicationException, IllegalActionException  {
        super(container, name);

        fileName = new FileParameter(this, "fileName");

	output = new TypedIOPort(this, "output", false, true);
	output.setTypeEquals(BaseType.STRING);

	parameters  = new TypedIOPort(this, "parameters", true, false);
	parameters.setMultiport(true);
	parameters.setTypeEquals(BaseType.STRING);

	// Set the type constraint.
	output.setTypeEquals(BaseType.STRING);

        _attachText("_iconDescription", "<svg>\n" +
                "<rect x=\"0\" y=\"0\" "
                + "width=\"60\" height=\"20\" "
                + "style=\"fill:white\"/>\n" +
                "</svg>\n");
    }

    // native functions written in c++
    public native int     openSeesInit();
    public native String  openSeesEval(String expression, int error);
    public native int     openSeesQuit();

    ///////////////////////////////////////////////////////////////////
    ////                     ports and parameters                  ////
    /** The output port and the input ports.
     */
    public TypedIOPort output = null;
    public TypedIOPort parameters = null;

    /** Thefilename
     */
    public FileParameter fileName;

    ///////////////////////////////////////////////////////////////////
    ////                         public methods                    ////

    /** Send the token in the value  parameter to the output.
     *  @exception IllegalActionException If it is thrown by the
     *   send() method sending out the token.
     */

    public boolean prefire() throws IllegalActionException {

	// load the shared object library
        System.loadLibrary("OpenSeesEvaluator");
	int res = this.openSeesInit();
	if (res != 0)
	    return false;
	return super.prefire();
    }

    public void fire() throws IllegalActionException {
        super.fire();
	//	String userNameStr =  userName.getToken().toString();
	//      userNameStr = userNameStr.substring(1, userNameStr.length() - 1);
	//	output.send(0, new StringToken("Hello YES " + userNameStr + "!"));

	String tclCmds = "";
	String result = "";
	int error = 0;

        // If the parameterValue input port is connected and has data, then
        // get the parameter value  from there
        for (int i=0; i<parameters.getWidth(); i++) {
	    if (error == 0) {
		if (parameters.hasToken(i)) {
		    
		    String name = ((StringToken) parameters.get(i))
                        .stringValue();
		    
		    // Using setExpression() rather than setToken() allows
		    // the string to refer to variables defined in the
		    // scope of this actor.
		    tclCmds = tclCmds + name + ";";
		    result = this.openSeesEval(name + ";", error);
		}
	    }
        }

	String fileNameStr =  fileName.getToken().toString();
	fileNameStr = fileNameStr.substring(1, fileNameStr.length() - 1);	
	if (fileNameStr.length() != 0 && error == 0) {
	    tclCmds += "source " + fileNameStr + ";";
	    result = this.openSeesEval("source " + fileNameStr, error);
	}	
	
	output.send(0, new StringToken(result));
    }
}
