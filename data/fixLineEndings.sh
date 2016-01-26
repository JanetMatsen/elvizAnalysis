#!/bin/bash

TMP=`mktemp /tmp/temp.XXXX`
echo $TMP

for file in *.csv
do
	awk 'gsub("\r", "\n")' $file > $TMP
	\mv $TMP $file
done		
