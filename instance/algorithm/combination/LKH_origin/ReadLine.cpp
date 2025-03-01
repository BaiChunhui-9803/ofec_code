#include "./INCLUDE/LKH.h"
#include<cstring>
/*      
 * The ReadLine function reads the next input line from a file. The function
 * handles the problem that an input line may be terminated by a carriage
 * return, a newline, both, or EOF.
 */
namespace LKH {
	static thread_local char *Buffer;
	static thread_local int MaxBuffer;

	static int EndOfLine(FILE * InputFile, int c)
	{
		int EOL = (c == '\r' || c == '\n');
		if (c == '\r') {
			c = fgetc(InputFile);
			if (c != '\n' && c != EOF)
				ungetc(c, InputFile);
		}
		return EOL;
	}

	char *LKHAlg::ReadLine(FILE * InputFile)
	{
		int i, c;

		if (Buffer == 0)
			assert(Buffer = (char *)malloc(MaxBuffer = 80));
		for (i = 0; (c = fgetc(InputFile)) != EOF && !EndOfLine(InputFile, c);
			i++) {
			if (i >= MaxBuffer - 1) {
				MaxBuffer *= 2;
				assert(Buffer = (char *)realloc(Buffer, MaxBuffer));
			}
			Buffer[i] = (char)c;
		}
		Buffer[i] = '\0';
		if (!LastLine || (int)strlen(LastLine) < i) {
			free(LastLine);
			assert(LastLine = (char *)malloc((i + 1) * sizeof(char)));
		}
		strcpy(LastLine, Buffer);
		return c == EOF && i == 0 ? 0 : Buffer;
	}
}