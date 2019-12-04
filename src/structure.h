#include <iostream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <ctype.h>
#include <zlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <inttypes.h>

#define ReadChunkSize 4000

using namespace std;

typedef uint64_t bwtint_t;
typedef unsigned char ubyte_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

typedef struct
{
	bwtint_t x[3];
} bwtintv_t;

typedef struct
{
	int len;
	int freq;
	bwtint_t* LocArr;
} bwtSearchResult_t;

typedef struct
{
	char* name; // chromosome name
	int64_t FowardLocation;
	int64_t ReverseLocation;
	int64_t len; // chromosome length
} Chromosome_t;

typedef struct
{
	int rPos; // read position
	int64_t gPos; // genome position
	int Len; // read block size
	int64_t PosDiff; // gPos-rPos
} SeedPair_t;

typedef struct
{
	int Score; // seed score
	vector<SeedPair_t> SeedVec;
} AlignmentCandidate_t;

typedef struct
{
	int rlen;
	char* seq;
	char* header;
} ReadItem_t;

// Global variables
extern bwt_t *Refbwt;
extern bwaidx_t *RefIdx;
extern vector<string> ReadVec;
extern map<string, string> REmap;
extern vector<int64_t> CutSiteVec;
extern unsigned char nst_nt4_table[256];
extern int64_t GenomeSize, TwoGenomeSize;
extern vector<Chromosome_t> ChromosomeVec;

extern const char* VersionStr;
extern map<int64_t, int> ChrLocMap;
extern char *IndexFileName, *OutputFileName;
extern vector<string> ReadFileNameVec1, ReadFileNameVec2;
extern bool bDebugMode, gzCompressed, FastQFormat, bSilent;
extern int MaxInsertSize, iThreadNum, iChromsomeNum, WholeChromosomeNum, MaxGaps, MinSeedLength;

// bwt_index.cpp
extern void RestoreReferenceInfo();
extern void bwa_idx_destroy(bwaidx_t *idx);
extern bwaidx_t *bwa_idx_load(const char *hint);

// bwt_search.cpp
extern bwtSearchResult_t BWT_Search(uint8_t* seq, int start, int stop);

// Mapping.cpp
extern void Mapping();

// Digest.cpp
extern void InitializeREmap();

// AlignmentCandidates.cpp
extern void RemoveShortSeeds(vector<SeedPair_t>& SeedVec, int thr);
extern bool CompByPosDiff(const SeedPair_t& p1, const SeedPair_t& p2);
extern bool CompByGenomePos(const SeedPair_t& p1, const SeedPair_t& p2);
extern vector<SeedPair_t> IdentifySeedPairs(int rlen, uint8_t* EncodeSeq);
extern void GenMappingReport(ReadItem_t& read, vector<AlignmentCandidate_t>& AlignmentVec);
extern vector<AlignmentCandidate_t> GenerateAlignmentCandidates(int rlen, vector<SeedPair_t> SeedPairVec);

// GetData.cpp
extern bool CheckReadFormat(const char* filename);
extern bool CheckBWAIndexFiles(string IndexPrefix);
extern int GetNextChunk(bool bSepLibrary, FILE *file, FILE *file2, ReadItem_t* ReadArr);
extern int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr);

// tools.cpp
extern void ShowSeedLocationInfo(int64_t MyPos);
extern int64_t GetAlignmentBoundary(int64_t gPos);
extern void ShowSeedInfo(vector<SeedPair_t>& SeedPairVec);
extern void FreeReadArrMemory(int ReadNum, ReadItem_t* ReadArr);
extern void GetComplementarySeq(int len, char* seq, char* rseq);
extern bool CheckAlignmentValidity(vector<SeedPair_t>& SeedPairVec);
extern bool CheckMappedPos(int real_pos1, int real_pos2, vector<AlignmentCandidate_t>& AlignmentVec);
extern void ShowAlignmentCandidateInfo(ReadItem_t* read, vector<AlignmentCandidate_t>& AlignmentVec);
//extern int ProcessSimpleSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);
