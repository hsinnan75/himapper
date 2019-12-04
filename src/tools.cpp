#include "structure.h" 

static const char ReverseMap[255] =
{
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
	'\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
	'\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
	'\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', /*  90 -  99 */
	'\0', '\0', '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
	'N',  '\0', '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', /* 110 - 119 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
	'\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
	'\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
};

void GetComplementarySeq(int len, char* seq, char* rseq)
{
	int i, j;

	for (j = len - 1, i = 0; i<j; i++, j--)
	{
		rseq[i] = ReverseMap[(int)seq[j]];
		rseq[j] = ReverseMap[(int)seq[i]];
	}
	if (i == j) rseq[i] = ReverseMap[(int)seq[i]];
}

//void SelfComplementarySeq(int len, char* seq)
//{
//	int i, j;
//	char aa1, aa2;
//
//	for (j = len - 1, i = 0; i<j; i++, j--)
//	{
//		aa1 = seq[i]; aa2 = seq[j];
//		seq[i] = ReverseMap[(int)aa2];
//		seq[j] = ReverseMap[(int)aa1];
//	}
//	if (i == j) seq[i] = ReverseMap[(int)seq[i]];
//}

int CalFragPairIdenticalBases(int len, char* frag1, char* frag2)
{
	int i, c;

	for (c = 0, i = 0; i < len; i++) if (frag1[i] != frag2[i]) c++;

	return (len - c);
}

int CalFragPairMismatchBases(int len, char* frag1, char* frag2)
{
	int i, c;

	for (c = 0, i = 0; i < len; i++) if (frag1[i] != frag2[i]) c++;

	return c;
}

void ShowSeedInfo(vector<SeedPair_t>& SeedPairVec)
{
	for (vector<SeedPair_t>::const_iterator iter = SeedPairVec.begin(); iter != SeedPairVec.end(); iter++)
		printf("\t\tseed#%d: R[%d-%d] G[%lld-%lld]\n", (int)(iter - SeedPairVec.begin() + 1), iter->rPos, iter->rPos + iter->Len - 1, (long long)iter->gPos, (long long)(iter->gPos + iter->Len - 1));
	printf("\n\n"); fflush(stdout);
}

void ShowSeedLocationInfo(vector<SeedPair_t>& SeedPairVec)
{
	int64_t gPos1, gPos2;

	if (SeedPairVec.begin()->gPos < GenomeSize)
	{
		gPos1 = SeedPairVec.begin()->gPos;
		gPos2 = SeedPairVec.rbegin()->gPos + SeedPairVec.rbegin()->Len - 1;
	}
	else
	{
		gPos1 = TwoGenomeSize - (SeedPairVec.rbegin()->gPos + SeedPairVec.rbegin()->Len);
		gPos2 = TwoGenomeSize - SeedPairVec.begin()->gPos - 1;
	}
	printf("Loc=%lld-%lld\n", (long long)gPos1, (long long)gPos2);
}

bool CheckMappedPos(int real_pos1, int real_pos2, vector<AlignmentCandidate_t>& AlignmentVec)
{
	bool bCheck = false;
	int64_t gPos1, gPos2;
	vector<AlignmentCandidate_t>::iterator iter;

	for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
	{
		if (iter->Score == 0) continue;
		if (iter->SeedVec.begin()->gPos < GenomeSize)
		{
			gPos1 = iter->SeedVec.begin()->gPos;
			gPos2 = iter->SeedVec.rbegin()->gPos + iter->SeedVec.rbegin()->Len - 1;
		}
		else
		{
			gPos1 = TwoGenomeSize - (iter->SeedVec.rbegin()->gPos + iter->SeedVec.rbegin()->Len);
			gPos2 = TwoGenomeSize - iter->SeedVec.begin()->gPos - 1;
		}
		if (abs(gPos1 - real_pos1) / 1000 == 0 || abs(gPos2 - real_pos1) / 1000 == 0 || abs(gPos1 - real_pos2) / 1000 == 0 || abs(gPos2 - real_pos2) / 1000 == 0)
		{
			bCheck = true;
			break;
		}
	}
	return bCheck;
}

void ShowAlignmentCandidateInfo(ReadItem_t* read, vector<AlignmentCandidate_t>& AlignmentVec)
{
	vector<AlignmentCandidate_t>::iterator iter;

	printf("%s\n", string().assign(100, '-').c_str()); printf("%s\n", read->header);
	for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
	{
		printf("Score=%d: ", iter->Score); ShowSeedLocationInfo(iter->SeedVec);
		ShowSeedInfo(iter->SeedVec);
	}
	//printf("%s\n\n", string().assign(100, '-').c_str());
	fflush(stdout);
}

int64_t GetAlignmentBoundary(int64_t gPos)
{
	map<int64_t, int>::iterator iter = ChrLocMap.lower_bound(gPos);

	return iter->first;
}

bool CheckFragValidity(SeedPair_t SeedPair)
{
	map<int64_t, int>::iterator iter1, iter2;

	iter1 = ChrLocMap.lower_bound(SeedPair.gPos);
	iter2 = ChrLocMap.lower_bound(SeedPair.gPos + SeedPair.Len - 1);

	return (iter1 != ChrLocMap.end() && iter2 != ChrLocMap.end() && iter1->first == iter2->first);
}

bool CheckAlignmentValidity(vector<SeedPair_t>& SeedPairVec)
{
	map<int64_t, int>::iterator iter1, iter2;

	iter1 = ChrLocMap.lower_bound(SeedPairVec.begin()->gPos);
	iter2 = ChrLocMap.lower_bound(SeedPairVec.rbegin()->gPos + SeedPairVec.rbegin()->Len - 1);

	return (iter1 != ChrLocMap.end() && iter2 != ChrLocMap.end() && iter1->first == iter2->first);
}

void FreeReadArrMemory(int ReadNum, ReadItem_t* ReadArr)
{
	for (int i = 0; i != ReadNum; i++)
	{
		delete[] ReadArr[i].header;
		delete[] ReadArr[i].seq;
	}
}
