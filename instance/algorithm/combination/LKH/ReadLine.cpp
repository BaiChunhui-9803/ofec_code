#include "./INCLUDE/LKH.h"
#include<cstring>
/*      
 * The ReadLine function reads the next input line from a file. The function
 * handles the problem that an input line may be terminated by a carriage
 * return, a newline, both, or EOF.
 */
namespace LKH {


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

	void LKHAlg::freeBuffer() {
		if (m_bufferInfo.Buffer) {
			free(m_bufferInfo.Buffer);
			m_bufferInfo.Buffer = 0;
		}
	}

	char *LKHAlg::ReadLine(FILE * InputFile)
	{
		int i, c;

		if (m_bufferInfo.Buffer == 0)
			assert(m_bufferInfo.Buffer = (char *)malloc(m_bufferInfo.MaxBuffer = 80));
		for (i = 0; (c = fgetc(InputFile)) != EOF && !EndOfLine(InputFile, c);
			i++) {
			if (i >= m_bufferInfo.MaxBuffer - 1) {
				m_bufferInfo.MaxBuffer *= 2;
				assert(m_bufferInfo.Buffer = (char *)realloc(m_bufferInfo.Buffer, m_bufferInfo.MaxBuffer));
			}
			m_bufferInfo.Buffer[i] = (char)c;
		}
		m_bufferInfo.Buffer[i] = '\0';
		if (!LastLine || (int)strlen(LastLine) < i) {
			free(LastLine);
			assert(LastLine = (char *)malloc((i + 1) * sizeof(char)));
		}
		strcpy(LastLine, m_bufferInfo.Buffer);
		return c == EOF && i == 0 ? 0 : m_bufferInfo.Buffer;
	}
}