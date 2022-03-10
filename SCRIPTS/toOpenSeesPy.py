# Convert an OpenSees(Tcl) script to OpenSeesPy
#  Author: Michael H. Scott, michael.scott@oregonstate.edu
#  Date: June 2018
#
# EXAMPLE
#   exec(open('toOpenSeesPy.py').read())
#   ...
#   outfile = open('model.py','w')
#   outfile.write('import openseespy.opensees as ops\n\n')
#                  input      output alias
#   toOpenSeesPy('model.tcl',outfile,'ops')
#   toOpenSeesPy('anotherScript.tcl',outfile,'ops')
#   ...
#   outfile.close()
#
#
# - The third argument is the alias put before each OpenSeesPy command.
#   Default is no alias
#      element('dispBeamColumn',...)
#   If alias is used, e.g., 'ops'
#      ops.element('dispBeamColumn',...)
#   You write your own import statement as shown in the example above
#
# - The calling Python script should open and close the file stream
#   for writing the converted .py file.  This allows you to call the
#   converter on multiple Tcl files in sequence, as shown in the example above. 
#
# - Assumes the OpenSees(.tcl) files defines commands line by line
#   without any loops, conditionals, variables, expressions, etc.
#   This is the format generated when you export a model from
#   OpenSees Navigator and perhaps from other front-ends to OpenSees.
#
# - If your OpenSees(.tcl) files do use loops, conditionals, variables,
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
# - Or insert output statements that create Python commands
#
#   set output [open model.py w]
#   ...
#   for {set i 1} {$i <= $N} {incr i} {
#     element truss $i ...
#     puts $output "ops.element('truss',$i, ...)"
#   }
#   ...
#   close $output
#
# - If you see any improvements to make to this toOpenSeesPy function,
#   please submit a pull request at github.com/OpenSees/OpenSees


# Helper function to deterine if a variable is a floating point value
#
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

# Function that does the conversion
#
def toOpenSeesPy(infile, outfile, alias=''):
    # Add a dot if needed
    if len(alias) > 0 and alias[-1] != '.':
        alias = alias + '.'
        
    infile = open(infile,'r')
    for line in infile:
        info = line.split()
        N = len(info)

        # Skip a blank line
        if N < 1:
            outfile.write('\n')
            continue        
        # Ignore a close brace
        if info[0][0] == '}':
            continue
        # Echo a comment line
        if info[0][0] == '#':
            outfile.write(line)        
            continue

        # Change print to printModel
        if info[0] == 'print':
            info[0] = 'printModel'

        # A command with no arguments, e.g., wipe
        if N == 1:
            outfile.write(f"{alias}{info[0]}()\n")
            continue

        # Needs to be a special case due to beam integration
        if N >= 8 and info[1] in ['nonlinearBeamColumn','forceBeamColumn','dispBeamColumn']:
            eleTag = info[2]
            secTag = info[6]
            # The original element format
            # element beamColumn tag ndI ndJ Np secTag transfTag
            #    0        1       2   3   4   5    6       7
            if isfloat(secTag):
                Np = info[5]
                transfTag = info[7]
                if info[1] == 'dispBeamColumn':
                    outfile.write(f"{alias}beamIntegration('Legendre',{eleTag},{secTag},{Np})\n")
                else:
                    outfile.write(f"{alias}beamIntegration('Lobatto',{eleTag},{secTag},{Np})\n")
            # The newer element format
            # element beamColumn tag ndI ndJ transfTag integrType ...
            #    0        1       2   3   4      5         6
            else:
                transfTag = info[5]
                outfile.write(f"{alias}beamIntegration('{info[6]}',{eleTag}")
                for j in range(7,N):
                    outfile.write(f',{info[j]}')
                outfile.write(')\n')
            if info[1] == 'nonlinearBeamColumn':
                info[1] = 'forceBeamColumn'
            outfile.write(f"{alias}element('{info[1]}',{eleTag},{info[3]},{info[4]},{transfTag},{eleTag})\n")                
            continue

        # Have to do the first argument before loop because of the commas
        if isfloat(info[1]):        
            outfile.write(f'{alias}{info[0]}({info[1]}')
        else:
            outfile.write(f"{alias}{info[0]}('{info[1]}'")
        # Now loop through the remaining arguments with preceding commas
        writeClose = True
        for i in range (2,N):
            if info[i] == '{':
                writeClose = True
                break
            if info[i] == '}':
                writeClose = False                
                break
            if isfloat(info[i]):
                outfile.write(f',{info[i]}')
            else:
                outfile.write(f",'{info[i]}'")
        if writeClose:
            outfile.write(')\n')        
    infile.close()

    outfile.write('\n\n')
    
