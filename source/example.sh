#!/bin/sh

TRY_WX=1

pathToImagesDir=".."
imagesDir="FITC-MCA0_N12_PI001_s1"
imStart=0
imEnd=61
outBase="../$imagesDir"
dx=0.412
dy=0.412
dz=0.492
unitName="microns"
intensityThresh=0.8

prog=minSurfTests

if [ ! -e ./$prog ]; then
	if [ $TRY_WX -eq 0 ]; then
		make noGUI > $outBase/make_out.txt
	else
		make wx > $outBase/make_out.txt
	fi
fi

if [ -e ./$prog ]; then
	nice ./$prog $pathToImagesDir/$imagesDir $imStart $imEnd $outBase $dx $dy $dz $unitName $intensityThresh > ${outBase}_log.txt
fi
