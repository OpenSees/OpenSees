#!/bin/bash
# Claudio Perez
# Loop over directories and create CMakeLists.txt in each.

# loop over directories
for i in */
do 
	# remove trailing slash
	i=${i%*/}
	name=`echo $i | sed 's|-|_|g'`;
	if [ -f $i/CMakeLists.txt ]; then
		:
		#echo $i;
		#cat $i/CMakeLists.txt;
	else 
		echo "OPS_Element_$name"
		{
		echo "add_library(OPS_Element_$name OBJECT)";
		echo "";
		cd "$i";
		echo "target_sources(OPS_Element_$name";
		echo "    PRIVATE";
		ls -1 *.[cfF]* | sed 's|^\([A-z]\)|        \1|g';
		echo "    PUBLIC";
		ls -1 *.h | sed 's|^\([A-z]\)|        \1|g';
		cd ../;
		echo ")";
		echo "target_include_directories(OPS_Element_$name PUBLIC \$(CMAKE_CURRENT_LIST_DIR))";
		echo "";
		} > $i/CMakeLists.txt 2>&1
	fi;
done

