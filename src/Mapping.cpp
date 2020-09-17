#include "structure.h"

#define MAPQ_COEF 50

typedef struct
{
	string header;
	int64_t gPos1, gPos2;
	bool strand1, strand2;
	int chr1, chr2, pos1, pos2, mapq1, mapq2;
} AlnReport_t;

int thread_count=0;
int64_t iDistance = 0;
int TotalPair, CorPair;
time_t StartProcessTime;
bool bSepLibrary = false;

vector<AlnReport_t> AlnReportVec;
gzFile gzReadFileHandler1, gzReadFileHandler2;
static pthread_mutex_t LibraryLock, OutputLock;
FILE *ReadFileHandler1, *ReadFileHandler2, *AlnSummaryFileHandler;
int64_t iTotalReadNum = 0, iUniqueMapping = 0, iUnMapping = 0, iPaired = 0, iMultiHits = 0;

bool CompByChr(const AlnReport_t& p1, const AlnReport_t& p2)
{
	if (p1.chr1 == p2.chr1)
	{
		if (p1.chr2 == p2.chr2)
		{
			if (p1.pos1 == p2.pos1) return p1.pos2 < p2.pos2;
			else return p1.pos1 < p2.pos1;
		}
		else return p1.chr2 < p2.chr2;
	}
	else return p1.chr1 < p2.chr1;
}

void RemoveRedundantCandidates(vector<AlignmentCandidate_t>& AlignmentVec)
{
	vector<AlignmentCandidate_t> vec;
	vector<AlignmentCandidate_t>::iterator iter;

	if (AlignmentVec.size() <= 1) return;
	else
	{
		int score = 0;
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++) if (iter->Score > score) score = iter->Score;
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++) if (iter->Score == score) vec.push_back(*iter);
	}
	if (vec.size() < AlignmentVec.size()) AlignmentVec.swap(vec);
}

int Cal_mapq(int rlen, int score)
{
	float similarity = 1.0 * score / rlen;
	return MAPQ_COEF* similarity;
}

AlnReport_t SwapAlnPair(AlnReport_t AlnReport)
{
	AlnReport_t retAlnRep;

	retAlnRep.chr1 = AlnReport.chr2; retAlnRep.chr2 = AlnReport.chr1;
	retAlnRep.pos1 = AlnReport.pos2; retAlnRep.pos2 = AlnReport.pos1;
	retAlnRep.mapq1 = AlnReport.mapq2; retAlnRep.mapq2 = AlnReport.mapq1;
	retAlnRep.strand1 = AlnReport.strand2; retAlnRep.strand2 = AlnReport.strand1;

	return retAlnRep;
}

AlnReport_t GenAlnReport(ReadItem_t& read, AlignmentCandidate_t& Aln1, AlignmentCandidate_t& Aln2, char* aln_summary)
{
	int64_t gPos1, gPos2;
	AlnReport_t AlnReport;
	map<int64_t, int>::iterator iter1, iter2;

	if (Aln1.SeedVec.begin()->gPos < GenomeSize) AlnReport.strand1 = true, gPos1 = Aln1.SeedVec.begin()->gPos;
	else AlnReport.strand1 = false, gPos1 = TwoGenomeSize - (Aln1.SeedVec.rbegin()->gPos + Aln1.SeedVec.rbegin()->Len);
	AlnReport.mapq1 = Cal_mapq(read.rlen, Aln1.Score); AlnReport.mapq2 = Cal_mapq(read.rlen, Aln2.Score);

	if (Aln2.SeedVec.begin()->gPos < GenomeSize) AlnReport.strand2 = true, gPos2 = Aln2.SeedVec.begin()->gPos;
	else AlnReport.strand2 = false, gPos2 = TwoGenomeSize - (Aln2.SeedVec.rbegin()->gPos + Aln2.SeedVec.rbegin()->Len);

	if ((iter1 = ChrLocMap.lower_bound(gPos1)) != ChrLocMap.end())
	{
		AlnReport.chr1 = iter1->second;
		AlnReport.pos1 = (int)(gPos1 - ChromosomeVec[AlnReport.chr1].FowardLocation);
	}
	else AlnReport.chr1 = -1, AlnReport.pos1 = 0;
	if ((iter2 = ChrLocMap.lower_bound(gPos2)) != ChrLocMap.end())
	{
		AlnReport.chr2 = iter2->second;
		AlnReport.pos2 = (int)(gPos2 - ChromosomeVec[AlnReport.chr2].FowardLocation);
	}
	else AlnReport.chr2 = -1, AlnReport.pos2 = 0;

	if (AlnReport.chr1 > AlnReport.chr2 || (AlnReport.chr1 == AlnReport.chr2 && AlnReport.pos1 > AlnReport.pos2)) AlnReport = SwapAlnPair(AlnReport);

	return AlnReport;
}

void *ReadMapping(void *arg)
{
	uint8_t* EncodeSeq;
	//int64_t readID_base;
	AlnReport_t AlnReport;
	char aln_summary[1024];
	ReadItem_t* ReadArr = NULL;
	pair<int, int> contact_pair;
	vector<AlnReport_t> myAlnReportVec;
	vector<SeedPair_t> SeedPairVec1, SeedPairVec2;
	vector<AlignmentCandidate_t> AlignmentVec1, AlignmentVec2;
	int i, j, k, aln_num1, aln_num2, ReadNum, myPairedMapping, myUnMapping, myMultiHits;

	myPairedMapping = myMultiHits = myUnMapping = 0; ReadArr = new ReadItem_t[ReadChunkSize];
	while (true)
	{
		pthread_mutex_lock(&LibraryLock);
		if(gzCompressed) ReadNum = gzGetNextChunk(bSepLibrary, gzReadFileHandler1, gzReadFileHandler2, ReadArr);
		else ReadNum = GetNextChunk(bSepLibrary, ReadFileHandler1, ReadFileHandler2, ReadArr);
		if (!bSilent) fprintf(stderr, "\r%lld reads have been processed in %ld seconds...", (long long)iTotalReadNum, (long)(time(NULL) - StartProcessTime)); fflush(stdout);
		iTotalReadNum += ReadNum;
		pthread_mutex_unlock(&LibraryLock);
		
		if (ReadNum == 0) break;
		if (ReadNum % 2 == 0)
		{
			//if (bDebugMode) printf("iDistance = %ld, iPaired=%d, EstiDistance=%d\n", iDistance, iPaired, EstDistance);
			for (i = 0, j = 1; i != ReadNum; i += 2, j += 2)
			{
				//if (bDebugMode) printf("Mapping paired reads#%d %s (len=%d) and %s (len=%d):\n", i + 1, ReadArr[i].header, ReadArr[i].rlen, ReadArr[j].header, ReadArr[j].rlen);
				EncodeSeq = new uint8_t[ReadArr[i].rlen]; for (k = 0; k < ReadArr[i].rlen;k++) EncodeSeq[k] = nst_nt4_table[(int)ReadArr[i].seq[k]];
				SeedPairVec1 = IdentifySeedPairs(ReadArr[i].rlen, EncodeSeq); AlignmentVec1 = GenerateAlignmentCandidates(ReadArr[i].rlen, SeedPairVec1);
				delete[] EncodeSeq;

				EncodeSeq = new uint8_t[ReadArr[j].rlen]; for (k = 0; k < ReadArr[j].rlen; k++) EncodeSeq[k] = nst_nt4_table[(int)ReadArr[j].seq[k]];
				SeedPairVec2 = IdentifySeedPairs(ReadArr[j].rlen, EncodeSeq); AlignmentVec2 = GenerateAlignmentCandidates(ReadArr[j].rlen, SeedPairVec2);
				delete[] EncodeSeq;

				RemoveRedundantCandidates(AlignmentVec1); RemoveRedundantCandidates(AlignmentVec2);

				aln_num1 = (int)AlignmentVec1.size(); aln_num2 = (int)AlignmentVec2.size();
				if (aln_num1 == 1 && aln_num2 == 1) // unique & paired-mapping
				{
					myPairedMapping += 2;
					AlnReport = GenAlnReport(ReadArr[i], AlignmentVec1[0], AlignmentVec2[0], aln_summary);
					if (AlnReport.chr1 != -1 && AlnReport.chr2 != -1)
					{
						AlnReport.header = ReadArr[i].header;
						myAlnReportVec.push_back(AlnReport);
					}
				}
				else if (aln_num1 > 1 && aln_num2 > 1) myMultiHits += 2;
				else
				{
					if (aln_num1 == 0) myUnMapping++; if (aln_num2 == 0) myUnMapping++;
				}
			}
		}
		FreeReadArrMemory(ReadNum, ReadArr);
		//pthread_mutex_lock(&DataLock);
		//pthread_mutex_unlock(&DataLock);
		//if (iTotalReadNum > 100000) break;
	}
	delete[] ReadArr;

	sort(myAlnReportVec.begin(), myAlnReportVec.end(), CompByChr);

	pthread_mutex_lock(&OutputLock);

	if (thread_count == 0) AlnReportVec.reserve((iTotalReadNum / 2));
	thread_count++;
	if (!bSilent) fprintf(stderr, "\33[2K\rMerge alignment results (%d / %d)...", thread_count, iThreadNum); fflush(stdout);
	iPaired += myPairedMapping; iMultiHits += myMultiHits; iUnMapping += myUnMapping;
	ReadNum = (int)myAlnReportVec.size();
	copy(myAlnReportVec.begin(), myAlnReportVec.end(), back_inserter(AlnReportVec));
	inplace_merge(AlnReportVec.begin(), AlnReportVec.end() - ReadNum, AlnReportVec.end(), CompByChr);
	pthread_mutex_unlock(&OutputLock);

	myAlnReportVec.clear();

	return (void*)(1);
}

void RemoveDuplications()
{
	int i, j, n, count;
	map<pair<int, int64_t>, bool> DupMap;

	n = (int)AlnReportVec.size();
	for (i = 0; i < n; i++)
	{
		for (count = 0, j = i + 1; j < n; j++)
		{
			if (AlnReportVec[i].pos1 == AlnReportVec[j].pos1) count++;
			else break;
		}
		if (count > MaxDuplicates) DupMap.insert(make_pair(make_pair(AlnReportVec[i].chr1, AlnReportVec[i].pos1), true));
		i += count;
	}
	fprintf(stderr, "Duplicated coordinates: %d\n", (int)DupMap.size());
	//for (map<pair<int, int64_t>, bool>::iterator iter = DupMap.begin(); iter != DupMap.end(); iter++) fprintf(stderr, "%s:%lld\n", (char*)ChromosomeVec[iter->first.first].name, (long long)iter->first.second);
	for (i = 0; i < n; i++)
	{
		if (DupMap.find(make_pair(AlnReportVec[i].chr1, AlnReportVec[i].pos1)) != DupMap.end()) AlnReportVec[i].mapq1 = 0;
		else if(DupMap.find(make_pair(AlnReportVec[i].chr2, AlnReportVec[i].pos2)) != DupMap.end()) AlnReportVec[i].mapq2 = 0;
	}
}

void Mapping()
{
	int i;
	vector<int> vec(iThreadNum); for (i = 0; i < iThreadNum; i++) vec[i] = i;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	for (MinSeedLength = 13; MinSeedLength < 16; MinSeedLength++) if (TwoGenomeSize < pow(4, MinSeedLength)) break;

	TotalPair = CorPair = 0;
	if (bDebugMode) iThreadNum = 1;
	if (bSilent) fprintf(stderr, "Start read mapping...\n");
	StartProcessTime = time(NULL);

	for (int LibraryID = 0; LibraryID < (int)ReadFileNameVec1.size(); LibraryID++)
	{
		gzReadFileHandler1 = gzReadFileHandler2 = NULL; ReadFileHandler1 = ReadFileHandler2 = NULL;

		if (ReadFileNameVec1[LibraryID].substr(ReadFileNameVec1[LibraryID].find_last_of('.') + 1) == "gz") gzCompressed = true;
		else gzCompressed = false;

		FastQFormat = CheckReadFormat(ReadFileNameVec1[LibraryID].c_str());
		//fprintf(stdout, "gz=%s, format=%s\n", gzCompressed ? "Yes" : "No", FastQFormat ? "Fastq" : "Fasta");

		if (gzCompressed) gzReadFileHandler1 = gzopen(ReadFileNameVec1[LibraryID].c_str(), "rb");
		else ReadFileHandler1 = fopen(ReadFileNameVec1[LibraryID].c_str(), "r");

		if (ReadFileNameVec1.size() == ReadFileNameVec2.size())
		{
			bSepLibrary = true;
			if (FastQFormat == CheckReadFormat(ReadFileNameVec2[LibraryID].c_str()))
			{
				if (gzCompressed) gzReadFileHandler2 = gzopen(ReadFileNameVec2[LibraryID].c_str(), "rb");
				else ReadFileHandler2 = fopen(ReadFileNameVec2[LibraryID].c_str(), "r");
			}
			else
			{
				fprintf(stderr, "Error! %s and %s are with different format...\n", (char*)ReadFileNameVec1[LibraryID].c_str(), (char*)ReadFileNameVec2[LibraryID].c_str());
				continue;
			}
		}
		else bSepLibrary = false;

		if (ReadFileHandler1 == NULL && gzReadFileHandler1 == NULL) continue;
		if (bSepLibrary && ReadFileHandler2 == NULL && gzReadFileHandler2 == NULL) continue;

		for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, ReadMapping, &vec[i]);
		for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

		if (gzCompressed)
		{
			if (gzReadFileHandler1 != NULL) gzclose(gzReadFileHandler1);
			if (gzReadFileHandler2 != NULL) gzclose(gzReadFileHandler2);
		}
		else
		{
			if (ReadFileHandler1 != NULL) fclose(ReadFileHandler1);
			if (ReadFileHandler2 != NULL) fclose(ReadFileHandler2);
		}
	}
	fprintf(stderr, "\rAll the %lld reads have been processed in %lld seconds.\n", (long long)iTotalReadNum, (long long)(time(NULL) - StartProcessTime));
	delete[] ThreadArr;

	if(iTotalReadNum > 0)
	{
		fprintf(stderr, "\t# of total mapped sequences = %lld (sensitivity = %.2f%%)\n", (long long)(iTotalReadNum - iUnMapping), (int)(10000 * (1.0*(iTotalReadNum - iUnMapping) / iTotalReadNum) + 0.5) / 100.0);
		fprintf(stderr, "\t# of unique paired alignments = %lld (%.2f%%)\n", (long long)iPaired, (int)(10000 * (1.0*iPaired / (iPaired + iMultiHits)) + 0.5) / 100.0);
		fprintf(stderr, "\t# of ambiguous paired alignments = %lld (%.2f%%)\n", (long long)iMultiHits, (int)(10000 * (1.0*iMultiHits / (iPaired + iMultiHits)) + 0.5) / 100.0);

		fprintf(stderr, "Remove duplicated reads...\n"); RemoveDuplications();
	}
	if (AlnReportVec.size() > 0)
	{
		fprintf(stderr, "Writing alignment summaries to [%s]\n", OutputFileName);
		AlnSummaryFileHandler = fopen(OutputFileName, "w");
		for (vector<AlnReport_t>::iterator iter = AlnReportVec.begin(); iter != AlnReportVec.end(); iter++)
		{
			if (iter->mapq1 == 0 || iter->mapq2 == 0) continue;
			fprintf(AlnSummaryFileHandler, "%s %d %s %d 1 %d %s %d 2 %d %d\n", iter->header.c_str(), (iter->strand1 ? 0 : 16), ChromosomeVec[iter->chr1].name, iter->pos1, (iter->strand2 ? 0 : 16), ChromosomeVec[iter->chr2].name, iter->pos2, iter->mapq1, iter->mapq2);
		}
		fclose(AlnSummaryFileHandler);
	}
}
