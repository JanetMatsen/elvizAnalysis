#include <stdio.h>

int main(int argc, char *argv[]) {
	int c;

	while ((c = getchar()) != EOF) {
		if (c <= 127)
			printf("%c", c);
	}

	return 0;
}
