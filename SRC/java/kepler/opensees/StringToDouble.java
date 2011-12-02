/* Convert a string to a double

 Copyright (c) 2005 The Regents of the University of California.
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
import ptolemy.actor.TypedAtomicActor;
import ptolemy.actor.TypedIOPort;
import ptolemy.data.DoubleToken;
import ptolemy.data.StringToken;
import ptolemy.data.type.BaseType;
import ptolemy.kernel.CompositeEntity;
import ptolemy.kernel.util.IllegalActionException;
import ptolemy.kernel.util.NameDuplicationException;

///////////////////////////////////////////////////////////////
/// StringToDouble

/** <p>This actor converts a string to a double.
 *  </p><p>
 *  It is based on the Ptolemy II StringToInt actor.</p>
 *  @author Frank McKenna, UC Berkeley, March 2007
 */
public class StringToDouble extends TypedAtomicActor {

    /** Construct an actor with the given container and name.
     *  @param container The container.
     *  @param name The name of this actor.
     *  @exception IllegalActionException If the actor cannot be contained
     *   by the proposed container.
     *  @exception NameDuplicationException If the container already has an
     *   actor with this name.
     */
    public StringToDouble(CompositeEntity container, String name)
        throws NameDuplicationException, IllegalActionException {
        super(container, name);

        string = new TypedIOPort(this, "string", true, false);
        string.setTypeEquals(BaseType.STRING);

        doublePort = new TypedIOPort(this, "double", false, true);
        doublePort.setTypeEquals(BaseType.DOUBLE);
    }

    ///////////////////////////////////////////////////////////////////
    ////                         public variables                  ////

    /** The input port, which has type <i>string</i>. */
    public TypedIOPort string;
    /** The output port, which has type <i>int</i>. */
    public TypedIOPort doublePort;

    ///////////////////////////////////////////////////////////////////
    ////                         public methods                    ////

    /** Consume one string token on the input port and produce a new
     *  integer token on the output port.
     *  @exception IllegalActionException If there is no direc/restor.
     */
    public void fire() throws IllegalActionException {
        super.fire();
        String _input = ((StringToken)string.get(0)).stringValue();
        Double _output  = new Double(_input);
        doublePort.send(0, new DoubleToken(_output.doubleValue()));
    }

    ///////////////////////////////////////////////////////////////////
    ////                         protected members                 ////

    ///////////////////////////////////////////////////////////////////
    ////                         private methods                   ////

    ///////////////////////////////////////////////////////////////////
    ////                         private members                   ////
}
