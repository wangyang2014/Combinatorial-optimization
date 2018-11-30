#!/bin/bash
#Attention: first argument filters the processed instances
if make
then
	echo "Compilation OK"
else
    echo "Please solve the compilation problems first".
	exit
fi

day=`date | cut -c 5-19|sed -e "s/://g;s/ /-/g"`
echo "Today is:$day"
function thisFile {
#        day=""
    outputDir=./output/$day/
	mkdir -p $outputDir
	mkdir -p $outputDir/tex/
	mkdir -p $outputDir/src/
	cp ./src/*.* $outputDir/src/ 
	cp Makefile $outputDir/src/ 
	cp mmdp.cfg $outputDir/src/
	cp  *sh $outputDir/src/
    slashOccurences=`echo "$1"| tr -dc '/' |wc -c`
    slashOccurences=`expr $slashOccurences + 1`
    instanceName=`echo $1|cut -f $slashOccurences -d "/"`
    echo -n "Processing $instanceName:"
    ./mmdp "$1" $outputDir > "$outputDir/$instanceName.txt" 
    #> "$outputDir/$instanceName"
}
find instances/ -maxdepth 2 -name "*.txt"| grep "$1"|while read LINE ; do
  thisFile "$LINE"
done
