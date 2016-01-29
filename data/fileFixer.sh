#!/bin/bash

tmp=`mktemp /tmp/fileFixer.sh.XXXXXXXXX`

# create a quick program to drop any character w/ byte val > 127
cat << EOF > fileFixer.c
#include <stdio.h>

int main(int argc, char *argv[]) {
	int c;

	while ((c = getchar()) != EOF) {
		if (c <= 127)
			printf("%c", c);
	}

	return 0;
}
EOF

gcc -o fileFixer fileFixer.c

for file in *.csv
do
	echo $file
	./fileFixer < $file > $tmp
	mv $tmp $file
done
