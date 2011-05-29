#include <stdio.h>

#define HEADERSIZE 17

int main()
{
	int i, c, x;
	
	for (i = 0; i < HEADERSIZE; i++)
		getchar();

	x = 0;
	while ((c = getchar()) != EOF)
	{
		putchar(c);
		x++;
		if (x == 3656)
		{
			for (; x < 3672; x++)
				putchar(0);
			x = 0;
		}
	}
	
	if (x != 0)
		fprintf(stderr, "WARNING: ended on column %d!\n", x);
	return 0;
}
