set currentDir [pwd]
set openSeesDir "$currentDir/.."
set coreDir  "$currentDir/core"
set coreList "$currentDir/coreList";


# procedure to go through the directories and clean out stuff
proc copyIt {openSeesDir fileName coreDir} {
    set dirFile $openSeesDir
    set dirFile [file join $dirFile $fileName]
    file copy -force $dirFile $coreDir
}


proc makeCore {openSeesDir coreList coreDir} {

    file mkdir $coreDir

    set pwd [pwd]
    set pwd [file join $pwd "core"]

    set inFileID [open $coreList r]

    foreach line [split [read $inFileID] \n] {

	if {[llength $line] == 0} {
	    # Blank line --> do nothing
	    continue
	} else {
	    copyIt $openSeesDir $line $coreDir
	}
    }

    puts "DONE copying files from $coreList to $coreDir"
}

makeCore $openSeesDir $coreList $coreDir