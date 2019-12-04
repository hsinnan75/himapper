#include "structure.h"

int iChromsomeNum;
map<int64_t, int> ChrLocMap;
int64_t GenomeSize, TwoGenomeSize;
vector<Chromosome_t> ChromosomeVec;

bool CheckReadFormat(const char* filename)
{
	char buf[1];
	gzFile file = gzopen(filename, "rb");
	gzread(file, buf, 1); gzclose(file);

	if (buf[0] == '@') return true; // fastq
	else return false;
}

int IdentifyHeaderBegPos(char* str, int len)
{
	for (int i = 1; i < len; i++)
	{
		if (str[i] != '>' && str[i] != '@') return i;
	}
	return len - 1;
}

int IdentifyHeaderEndPos(char* str, int len)
{
	if (len > 100) len = 100;
	for (int i = 1; i < len; i++)
	{
		if (str[i] == ' ' || str[i] == '/' || isprint(str[i]) == 0) return i;
	}
	return len - 1;;
}

ReadItem_t GetNextEntry(FILE *file)
{
	int p1, p2;
	ssize_t len;
	size_t size = 0;
	ReadItem_t read;
	char *buffer = NULL;

	read.header = read.seq = NULL; read.rlen = 0;

	if ((len = getline(&buffer, &size, file)) != -1)
	{
		p1 = IdentifyHeaderBegPos(buffer, len); p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
		read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
		//len -= 1; read.header = new char[len]; strncpy(read.header, (buffer + 1), len - 1); read.header[len - 1] = '\0';
		if (FastQFormat)
		{
			if ((read.rlen = getline(&buffer, &size, file)) != -1)
			{
				read.rlen -= 1; read.seq = new char[read.rlen + 1]; strncpy(read.seq, buffer, read.rlen); read.seq[read.rlen] = '\0';
				getline(&buffer, &size, file); getline(&buffer, &size, file);
				
			}
			else read.rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read.rlen = (int)seq.length()) > 0)
			{
				read.seq = new char[read.rlen + 1];
				strcpy(read.seq, (char*)seq.c_str());
				read.seq[read.rlen] = '\0';
			}
		}
	}
	free(buffer);

	return read;
}

int GetRealPos(bool b, string str)
{
	string tmp;
	int len, p, pos;

	len = str.length(); 
	if (str[len - 1] == 'F' || str[len - 1] == 'R')
	{
		if (b == 0)
		{
			p = str.find("NC_000913:") + 10;
			tmp = str.substr(p, str.find_first_of('.', p) - p); pos = atoi(tmp.c_str());
		}
		else
		{
			p = str.find_last_of("..") + 2;
			tmp = str.substr(p, str.find_first_of(':', p) - p); pos = atoi(tmp.c_str());
		}
	}
	else
	{
		if (b == 0) str.resize(str.find_last_of(' '));
		p = str.find_last_of(':') + 1; tmp = str.substr(p); pos = atoi(tmp.c_str());
	}
	return pos;
}

int GetNextChunk(bool bSepLibrary, FILE *file, FILE *file2, ReadItem_t* ReadArr)
{
	int iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = GetNextEntry(file)).rlen == 0) break;
		//ReadArr[iCount].real_pos = GetRealPos(iCount % 2, ReadArr[iCount].header);
		
		iCount++;
		if (bSepLibrary) ReadArr[iCount] = GetNextEntry(file2);
		else ReadArr[iCount] = GetNextEntry(file);

		if (ReadArr[iCount].rlen == 0) break;

		//ReadArr[iCount].real_pos = GetRealPos(iCount % 2, ReadArr[iCount].header);
		if (++iCount == ReadChunkSize) break;
	}
	return iCount;
}

ReadItem_t gzGetNextEntry(gzFile file)
{
	int p1, p2;
	char* buffer;
	ReadItem_t read;
	int len, buf_size = 1000;
	
	buffer = new char[buf_size]; read.header = read.seq = NULL; read.rlen = 0;
	if (gzgets(file, buffer, buf_size) != NULL)
	{
		len = strlen(buffer); 
		if (len > 0 && (buffer[0] == '@' || buffer[0] == '>'))
		{
			p1 = IdentifyHeaderBegPos(buffer, len); p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
			read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
			//len -= 1; read.header = new char[len]; strncpy(read.header, buffer + 1, len - 1); read.header[len - 1] = '\0';
			//read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
			gzgets(file, buffer, buf_size); read.rlen = strlen(buffer) - 1; read.seq = new char[read.rlen + 1]; read.seq[read.rlen] = '\0';
			strncpy(read.seq, buffer, read.rlen);

			if (FastQFormat) gzgets(file, buffer, buf_size), gzgets(file, buffer, buf_size);
		}
	}
	delete[] buffer;

	return read;
}

int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr)
{
	int iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = gzGetNextEntry(file)).rlen == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		iCount++;

		if (bSepLibrary) ReadArr[iCount] = gzGetNextEntry(file2);
		else ReadArr[iCount] = gzGetNextEntry(file);

		if (ReadArr[iCount].rlen == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];

		iCount++;
		if (iCount == ReadChunkSize) break;
	}
	return iCount;
}


bool CheckBWAIndexFiles(string IndexPrefix)
{
	fstream file;
	string filename;
	bool bChecked = true;

	filename = IndexPrefix + ".ann"; file.open(filename.c_str(), ios_base::in);
	if (!file.is_open()) return false; else file.close();

	filename = IndexPrefix + ".amb"; file.open(filename.c_str(), ios_base::in);
	if (!file.is_open()) return false; else file.close();

	filename = IndexPrefix + ".pac"; file.open(filename.c_str(), ios_base::in);
	if (!file.is_open()) return false; else file.close();

	return bChecked;
}
