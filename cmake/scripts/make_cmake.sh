#!/bin/bash
# Claudio Perez
# Loop over directories and create CMakeLists.txt in each.

# loop over directories
for i in */
do 
	# remove trailing slash
	i=${i%*/}
	#name=`echo $i | sed 's|-|_|g' | sed -r 's/\<./\U&/g'`;
	name="Renderer";
	if [ -f $i/CMakeLists.txt ]; then
		:
		#echo $i;
		#cat $i/CMakeLists.txt;
	else 
		#echo "OPS_$name"
		{
		echo "add_library(OPS_$name OBJECT)";
		echo "";
		cd "$i";
		echo "target_sources(OPS_$name";
		echo "    PRIVATE";
		ls -1 *.[cfF]* | sed 's|^\([A-z]\)|        \1|g' | grep -v Tcl;
		echo "    PUBLIC";
		ls -1 *.h | sed 's|^\([A-z]\)|        \1|g' | grep -v Tcl;
		echo ")";
		echo "target_include_directories(OPS_$name PUBLIC \${CMAKE_CURRENT_LIST_DIR})";
		echo "";
		for d in ./*/
		do
			echo "add_subdirectory($d)" | grep -v '/\*/';
		done
		cd ../;
		} 2>&1 #> $i/CMakeLists.txt
	fi;
done

