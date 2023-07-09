/* OpenSeesParameter actor for using OpenSees with Ptolemy II.

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
// $Source: /usr/local/cvs/OpenSees/SRC/java/kepler/opensees/OpenSeesParameter.java,v $
                                                                        
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

//////////////////////////////////////////////////////////////////////////
//// OpenSeesParameter
/**
 * OpenSeesParameter: 
"This is the implementation of a OpenSeesParameter actor using Ptolemy II. 
This actor has a parameter, an input port and an output port. The output
is the concatanation of the string parameter and the string in the input port.
 *  @author Frank McKenna, UC Berkeley, March 2007
*/

public class OpenSeesParameter extends TypedAtomicActor {

   /** Construct a OpenSeesParameter source with the given container and name.
     *  @param name The name of this actor.
     *  @exception IllegalActionException If the entity cannot be contained
     *   by the proposed container.
     *  @exception NameDuplicationException If the container already has an
     *   actor with this name.
     */
    public OpenSeesParameter(CompositeEntity container, String name)
            throws NameDuplicationException, IllegalActionException  {
        super(container, name);

        parameterName = new FileParameter(this, "parameterName");
        parameterValue = new FileParameter(this, "parameterValue");

	output = new TypedIOPort(this, "output", false, true);
	output.setTypeEquals(BaseType.STRING);

	parameterValuePort  = new TypedIOPort(this, "parameterValuePort", true, false);
	parameterValuePort.setTypeEquals(BaseType.STRING);

	// Set the type constraint.
	output.setTypeEquals(BaseType.STRING);

        _attachText("_iconDescription", "<svg>\n" +
                "<rect x=\"0\" y=\"0\" "
                + "width=\"60\" height=\"20\" "
                + "style=\"fill:white\"/>\n" +
                "</svg>\n");
    }

    ///////////////////////////////////////////////////////////////////
    ////                     ports and parameters                  ////
    /** The input and output port. The type of this port will be set to String.
     */
    public TypedIOPort output = null;
    public TypedIOPort parameterValuePort = null;

    /** The parameters
     */
    public FileParameter parameterName;
    public FileParameter parameterValue;

    
    ///////////////////////////////////////////////////////////////////
    ////                         public methods                    ////

    /** Send the token in the value  parameter to the output.
     *  @exception IllegalActionException If it is thrown by the
     *   send() method sending out the token.
     */

    public void fire() throws IllegalActionException {
        super.fire();

        // If the parameterValue input port is connected and has data, then
        // get the parameter value  from there
        if (parameterValuePort.getWidth() > 0) {
	    if (parameterValuePort.hasToken(0)) {
                String name = ((StringToken) parameterValuePort.get(0))
                        .stringValue();

                // Using setExpression() rather than setToken() allows
                // the string to refer to variables defined in the
                // scope of this actor.
                parameterValue.setExpression(name);
            }
        }

	String parameterNameStr =  parameterName.getToken().toString();
	parameterNameStr = parameterNameStr.substring(1, parameterNameStr.length() - 1);
	String parameterValueStr =  parameterValue.getToken().toString();
	parameterValueStr = parameterValueStr.substring(1, parameterValueStr.length() - 1);

	output.send(0, new StringToken(parameterNameStr + " " + parameterValueStr));
    }
}
