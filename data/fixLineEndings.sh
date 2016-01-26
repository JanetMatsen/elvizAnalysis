#!/bin/bash

TMP=`mktemp`

for file in *.csv
do
	awk 'gsub("\r", "\n")' $file > $TMP
	\mv $TMP $file
done		
