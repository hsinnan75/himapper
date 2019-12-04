#include "structure.h"

map<string, string> REmap;
vector<int64_t> CutSiteVec;

void InitializeREmap()
{
	REmap.insert(make_pair("mboi", "^GATC"));
	REmap.insert(make_pair("dpnii", "^GATC"));
	REmap.insert(make_pair("bglii", "A^GATCT"));
	REmap.insert(make_pair("hindiii", "A^AGCTT"));
}

void DigestGenome(string enzyme)
{

}
