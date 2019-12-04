#include "structure.h"

bool CompByPosDiff(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.rPos < p2.rPos;
	else return p1.PosDiff < p2.PosDiff;
}

bool CompByGenomePos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.rPos < p2.rPos;
	else return p1.gPos < p2.gPos;
}

bool CompByReadPos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	return p1.rPos < p2.rPos;
}

vector<SeedPair_t> IdentifySeedPairs(int rlen, uint8_t* EncodeSeq)
{
	SeedPair_t SeedPair;
	int i, pos, end_pos;
	vector<SeedPair_t> SeedPairVec;
	bwtSearchResult_t bwtSearchResult;

	pos = 0, end_pos = rlen - MinSeedLength;
	while (pos < end_pos)
	{
		if (EncodeSeq[pos] > 3) pos++;
		else
		{
			bwtSearchResult = BWT_Search(EncodeSeq, pos, rlen);
			if (bwtSearchResult.freq > 0)
			{
				SeedPair.rPos = pos; SeedPair.Len =  bwtSearchResult.len;
				for (i = 0; i != bwtSearchResult.freq; i++)
				{
					SeedPair.PosDiff = (SeedPair.gPos = bwtSearchResult.LocArr[i]) - SeedPair.rPos;
					if (SeedPair.PosDiff >= 0) SeedPairVec.push_back(SeedPair);
				}
				delete[] bwtSearchResult.LocArr;
			}
			pos += (bwtSearchResult.len + 1);
		}
	}
	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByPosDiff);

	return SeedPairVec;
}

vector<AlignmentCandidate_t> GenerateAlignmentCandidates(int rlen, vector<SeedPair_t> SeedPairVec)
{
	int64_t gPos_end;
	int i, j, k, thr, num;
	AlignmentCandidate_t AlignmentCandidate;
	vector<AlignmentCandidate_t> AlignmentVec;

	thr = (rlen >> 1); num = (int)SeedPairVec.size();
	i = 0; while (i < num && SeedPairVec[i].PosDiff < 0) i++;
	for (i = 0; i < num;)
	{
		AlignmentCandidate.Score = SeedPairVec[i].Len; gPos_end = GetAlignmentBoundary(SeedPairVec[i].gPos); AlignmentCandidate.SeedVec.clear();
		for (j = i, k = i + 1; k < num; k++)
		{
			if ((SeedPairVec[k].PosDiff - SeedPairVec[j].PosDiff) > MaxGaps || SeedPairVec[k].gPos > gPos_end) break;
			else
			{
				AlignmentCandidate.Score += SeedPairVec[k].Len;
				j = k;
			}
		}
		if (AlignmentCandidate.Score > thr)
		{
			copy(SeedPairVec.begin() + i, SeedPairVec.begin() + k, back_inserter(AlignmentCandidate.SeedVec));
			sort(AlignmentCandidate.SeedVec.begin(), AlignmentCandidate.SeedVec.end(), CompByGenomePos);
			AlignmentVec.push_back(AlignmentCandidate);
		}
		i = k;
	}
	return AlignmentVec;
}

void GenMappingReport(ReadItem_t& read, vector<AlignmentCandidate_t>& AlignmentVec)
{
	//int i, j, num, s;
	//map<int, int>::iterator iter;
	//map<int64_t, int>::iterator ChrIter;

	//if (bDebugMode) printf("\n\n%s\nGenerate alignment for read %s (%d cans)\n", string().assign(100, '=').c_str(), read.header, (int)AlignmentVec.size()), fflush(stdout);
	
	//read.score = 0;
	//if ((read.CanNum = (int)AlignmentVec.size()) > 0)
	//{
	//	read.AlnReportArr = new AlignmentReport_t[read.CanNum];
	//	for (i = 0; i != read.CanNum; i++)
	//	{
	//		read.AlnReportArr[i].AlnScore = 0;
	//		read.AlnReportArr[i].PairedAlnCanIdx = AlignmentVec[i].PairedAlnCanIdx;

	//		if (AlignmentVec[i].Score == 0) continue;
	//		IdentifyNormalPairs(read.rlen, -1, AlignmentVec[i].SeedVec); // fill missing framgment pairs (normal pairs) between simple pairs
	//		//if (bDebugMode) printf("Process candidate#%d (Score = %d, SegmentPair#=%d): \n", i + 1, AlignmentVec[i].Score, (int)AlignmentVec[i].SeedVec.size()), ShowSeedInfo(AlignmentVec[i].SeedVec);
	//		if (CheckCoordinateValidity(AlignmentVec[i].SeedVec) == false) continue;

	//		cigar_vec.clear();
	//		for (num = (int)AlignmentVec[i].SeedVec.size(), j = 0; j != num; j++)
	//		{
	//			if (AlignmentVec[i].SeedVec[j].rLen == 0 && AlignmentVec[i].SeedVec[j].gLen == 0) continue;
	//			else if (AlignmentVec[i].SeedVec[j].bSimple)
	//			{
	//				//if (bDebugMode) ShowFragmentPair(AlignmentVec[i].SeedVec[j]);
	//				cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[j].rLen, 'M'));
	//				read.AlnReportArr[i].AlnScore += AlignmentVec[i].SeedVec[j].rLen;
	//			}
	//			else
	//			{
	//				//if (bDebugMode) printf("Check normal pair#%d: R[%d-%d]=%d G[%lld-%lld]=%d\n", j + 1, AlignmentVec[i].SeedVec[j].rPos, AlignmentVec[i].SeedVec[j].rPos + AlignmentVec[i].SeedVec[j].rLen - 1, AlignmentVec[i].SeedVec[j].rLen, AlignmentVec[i].SeedVec[j].gPos, AlignmentVec[i].SeedVec[j].gPos + AlignmentVec[i].SeedVec[j].gLen - 1, AlignmentVec[i].SeedVec[j].gLen);
	//				if (j == 0)
	//				{
	//					s = ProcessHeadSequencePair(read.seq, AlignmentVec[i].SeedVec[0], cigar_vec);
	//					read.AlnReportArr[i].AlnScore += s;
	//					if (s == 0)
	//					{
	//						AlignmentVec[i].SeedVec[0].gPos = AlignmentVec[i].SeedVec[1].gPos;
	//						AlignmentVec[i].SeedVec[0].rLen = AlignmentVec[i].SeedVec[0].gLen = 0;
	//					}
	//				}
	//				else if (j == num - 1)
	//				{
	//					s = ProcessTailSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
	//					read.AlnReportArr[i].AlnScore += s;
	//					if (s == 0)
	//					{
	//						AlignmentVec[i].SeedVec[j].gPos = AlignmentVec[i].SeedVec[j - 1].gPos + AlignmentVec[i].SeedVec[j - 1].gLen - 1;
	//						AlignmentVec[i].SeedVec[j].rLen = AlignmentVec[i].SeedVec[j].gLen = 0;
	//					}
	//				}
	//				else
	//				{
	//					read.AlnReportArr[i].AlnScore += ProcessNormalSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
	//				}
	//			}
	//		}
	//		//if (bDebugMode) printf("Alignment score = %d (rlen=%d) \n", read.AlnReportArr[i].AlnScore, read.rlen), fflush(stdout);
	//		if (cigar_vec.size() > 1)
	//		{
	//			read.AlnReportArr[i].AlnScore -= GapPenalty(cigar_vec);
	//			if (read.AlnReportArr[i].AlnScore <= 0)
	//			{
	//				read.AlnReportArr[i].AlnScore = 0;
	//				continue;
	//			}
	//		}
	//		if(cigar_vec.size() == 0 || (read.AlnReportArr[i].coor = GenCoordinateInfo(AlignmentVec[i].SeedVec[0].gPos, (AlignmentVec[i].SeedVec[num - 1].gPos + AlignmentVec[i].SeedVec[num - 1].gLen - 1), cigar_vec)).gPos <= 0) read.AlnReportArr[i].AlnScore = 0;

	//		if (read.AlnReportArr[i].AlnScore > read.score)
	//		{
	//			read.iBestAlnCanIdx = i;
	//			read.sub_score = read.score;
	//			read.score = read.AlnReportArr[i].AlnScore;
	//		}
	//		else if (read.AlnReportArr[i].AlnScore == read.score)
	//		{
	//			read.sub_score = read.score;
	//			if (!bMultiHit && ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].len > ChromosomeVec[read.AlnReportArr[read.iBestAlnCanIdx].coor.ChromosomeIdx].len) read.iBestAlnCanIdx = i;
	//		}
	//	}
	//}
	//else
	//{
	//	read.CanNum = 1; read.iBestAlnCanIdx = 0;
	//	read.AlnReportArr = new AlignmentReport_t[1];
	//	read.AlnReportArr[0].AlnScore = 0;
	//	read.AlnReportArr[0].PairedAlnCanIdx = -1;
	//}
}
