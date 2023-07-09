# procedure to go through the directories and clean out stuff

proc cleanIt {startDir} {

    puts stderr $startDir

    set pwd [pwd]
    if [catch {cd $startDir} err] {
	puts stderr $err
	return
    }

    # remove unwanted stuff
    foreach file [glob -nocomplain *.o] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }
    foreach file [glob -nocomplain *.obj] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }
    foreach file [glob -nocomplain *~] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }
    foreach file [glob -nocomplain *.bak] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }
    foreach file [glob -nocomplain *.log] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }
    foreach file [glob -nocomplain *.LOG] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }
    foreach file [glob -nocomplain #*#] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }
    foreach file [glob -nocomplain .#*] {
	if [catch [exec rm -fr $file] error1] {set a 1}
    }

    if [catch [exec rm -fr core] error1] {set a 1}
    if [catch [exec rm -fr cxx_repository] error1] {set a 1}
    if [catch [exec rm -fr *.out] error1] {set a 1}
    if [catch [exec rm -fr ti_files] error1] {set a 1}

    # recursivily do the proc for any subdirectories
    # except those in the directories defined in DirNoModify
    foreach dirFile [glob -nocomplain  *] {
	if [file isdirectory $dirFile] {
	    cleanIt [file join $startDir $dirFile]
	}
    }
    cd $pwd
}
