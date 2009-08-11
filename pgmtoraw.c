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
		if ((x++) == 3656)
		{
			x = 0;
			for (; x < 3672; x++)
				putchar(128);
		}
	}
	return 0;
}
