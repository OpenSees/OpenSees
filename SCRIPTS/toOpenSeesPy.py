# Convert an OpenSees(Tcl) script to OpenSeesPy
#  Author: Michael H. Scott, michael.scott@oregonstate.edu
#  Date: June 2018
#
# Usage in a Python script
#   exec(open('toOpenSeesPy.py').read())
#   ...
#   outfile = open('model.py','w')
#   toOpenSeesPy('model.tcl',outfile)
#   toOpenSeesPy('anotherScript.tcl',outfile)
#   ...
#   outfile.close()
#
# - Assumes the OpenSees(.tcl) file defines the model line by line
#   without any loops, conditionals, variables, expressions, etc.
#   This is the format generated when you export a model from
#   OpenSees Navigator and perhaps from other front-ends to OpenSees.
#
# - The calling Python script should open and close the file stream for
#   for writing the converted .py file.  This allows you to call the
#   converter on multiple Tcl files in sequence, as shown above. 
#
# - If your OpenSees(.tcl) file uses any loops, conditionals, variables,
#   expressions, etc., you might be better off to port your OpenSees
#   model from Tcl to Python manually, or you can look in to Tkinter.
#
# - You may have some luck making your own "middleware" to convert your
#   OpenSees(.tcl) script to a model defined line by line by inserting
#   output statements in your loops and other constructs.  Even though
#   this won't get you to 100% conversion and you'll still have some
#   conversions to make here and there, it'll get you pretty far.
#
#   set output [open lineByLine.tcl w]
#   ...
#   for {set i 1} {$i <= $N} {incr i} {
#     element truss $i ...
#     puts $output "element truss $i ..."
#   }
#   ...
#   close $output
#
#   Then, in your Python script, call toOpenSeesPy with lineByLine.tcl as
#   the input file.
#
# - If you see any improvements to make to this toOpenSeesPy function,
#   please submit a pull request at OpenSees/OpenSees on github


# Helper function to deterine if a variable is a floating point value or not
#
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

# Function that does the conversion
#
def toOpenSeesPy(infile, outfile):
    outfile.write('\n\n')
    infile = open(infile,'r')
    for line in infile:
        info = line.split()
        N = len(info)
        
        # Ignore a close brace
        if N > 0 and info[0][0] == '}':
            continue
        # Echo a comment line
        if N < 2 or info[0][0] == '#':
            outfile.write(line)        
            continue
        
        # Needs to be a special case for now due to beam integration
        if info[1] == 'forceBeamColumn' or info[1] == 'dispBeamColumn':
            secTag = info[6]
            eleTag = info[2]
            Np = info[5]
            Np = 3
            outfile.write('ops.beamIntegration(\'Legendre\',%s,%s,%s)\n' % (eleTag,secTag,Np))
            outfile.write('ops.element(\'%s\',%s,%s,%s,%s,%s)\n' % (info[1],eleTag,info[3],info[4],info[7],eleTag))
            continue

        # Change print to printModel
        if info[0] == 'print':
            info[0] = 'printModel'
        
        # For everything else, have to do the first one before loop because of the commas
        if isfloat(info[1]):        
            outfile.write('ops.%s(%s' % (info[0],info[1]))            
        else:
            outfile.write('ops.%s(\'%s\'' % (info[0],info[1]))
        # Now loop through the rest with preceding commas
        writeClose = True
        for i in range (2,N):
            if info[i] == '{':
                writeClose = True
                break
            if info[i] == '}':
                writeClose = False                
                break
            if isfloat(info[i]):
                outfile.write(',%s' % info[i])
            else:
                outfile.write(',\'%s\'' % info[i])
        if writeClose:
            outfile.write(')\n')        
    infile.close()
