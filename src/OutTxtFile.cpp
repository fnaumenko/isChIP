/**********************************************************
OutTxtFile.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 03.04.2019
-------------------------
Provides output text files functionality
***********************************************************/

#include "OutTxtFile.h"

/************************ class Random ************************/
int Random::Seed;

// Sets  and returns real seed
//	@random: if true, random seed
int Random::SetSeed(UINT seed)
{
	if(seed)
		Seed = 12345678 + (seed<<24);	// any number with capacity 6-8 (needed RAND_XORSHIFT constructor)
	else {
		time_t tm;
		time(&tm);	// get current time; same as: timer = time(NULL)
		struct tm y2k = {0};
		y2k.tm_year = 117; y2k.tm_mday = 1;
		Seed = int(difftime(tm, mktime(&y2k)));	// seconds since January 1, 2017
	}
	return Seed;
}

Random::Random()
{
#ifdef RAND_STD
	srand( (unsigned)Seed );
	_seed = Seed;
#elif defined RAND_MT
	mt[0]= Seed;
	for (mti=1; mti < MERS_N; mti++)
		mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
	for (int i = 0; i < 37; i++) rand();		// Randomize some more
#elif defined RAND_XORSHIFT
	x = Seed;
	// initialize to fix random generator. Any initialization of y, w, z in fact
	y = x >> 1;	 w = y + 1000;  z = w >> 1;
#endif
}

#ifdef RAND_MT
// Generates 32 random bits
uint32_t Random::rand()
{
	if (mti >= MERS_N) {
		// Generate MERS_N words at one time
		const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
		const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
		static const uint32_t mag01[2] = {0, MERS_A};
		int kk;
		for (kk=0; kk < MERS_N-MERS_M; kk++) {    
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

		for (; kk < MERS_N-1; kk++) {    
			y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
			mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

		y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
		mti = 0;
	}
	y = mt[mti++];

	// Tempering (May be omitted):
	y ^=  y >> MERS_U;
	y ^= (y << MERS_S) & MERS_B;
	y ^= (y << MERS_T) & MERS_C;
	y ^=  y >> MERS_L;

	return y;
}
#elif defined RAND_XORSHIFT
// Generates 32 random bits
uint32_t Random::rand()
{
    uint32_t t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = w ^ (w >> 19) ^ t ^ (t >> 8);
}
#endif	

// Returns random integer within interval [1, max]
int Random::Range(int max) {
#ifdef RAND_MT
	if( max == 1 )	return 1;
	return min(int(DRand() * (max-1) + 1), max);
#else
	return int(DRand() * --max + 1);
#endif
}

// Returns true with given likelihood
//	@sample: probability of returning true; from 0.0. to 1.0
bool Random::Sample(float sample) 
{
	return sample >= 1.0 ? true : (sample == 0.0 ? false : (DRand() <= sample));
}

/************************ end of class Random ************************/

const char* GM::title[] = {"test","control"};	// title: printed member's name

float DistrParams::lnMean;	// mean of initial lognormal distribution
float DistrParams::lnSigma;	// sigma of initial lognormal distribution
float DistrParams::ssMean;	// mean of size selection normal distribution
float DistrParams::ssSigma;	// sigma of size selection normal distribution

/************************ class ReadQualPattern ************************/

// Creates Read quality pattern buffer and fills it by first valid line from file.
//	@rqPattFName: name of valid file with a quality line
ReadQualPattern::ReadQualPattern(const char* rqPattFName) : _rqPatt(NULL), _rqTempl(NULL)
{
	if(!rqPattFName) return;

	TabFile file(rqPattFName);
	const char* line = file.GetLine();	// first line in file
	if(!line)	Err("no lines", rqPattFName).Throw();
	_rqPatt = new char[Read::Len];		// don't check allocation: small size
	readlen lineLen = (readlen)file.LineLength();
	if(lineLen < Read::Len)
		Read::FillBySeqQual(_rqPatt);	// keep the rest if buffer filled by default quality
	else	lineLen = Read::Len;	// cut off LineLen in any case
	memcpy(_rqPatt, line, lineLen);
}

// Fills external FQ|SAM template by Read quality pattern and remembers it
//	@templ: pointer to external FQ|SAM Read quality template
void ReadQualPattern::Fill(char* templ)
{ 
	if(_rqPatt)	memcpy(_rqTempl = templ, _rqPatt, Read::Len);
	else		Read::FillBySeqQual(templ);
}

/************************ end of class ReadQualPattern ************************/

/************************ struct ReadName ************************/

BYTE ReadName::MaxLen = 0;	// Initialized in ExtOutFile::Init();
							// if initialize by declaration, seg fault by compiling with gcc 4.1.2

void ReadName::Init()
{
	MaxLen =
		Product::Title.length() + 1 +	// + 1 delimiter
		Chrom::MaxAbbrNameLength + 1 +	// length of chrom's name + 1 delimiter
		2*CHRLEN_CAPAC + 1 +			// length of PE Read name + 1 delimiter
		ExtOutFile::MateLen;			// length of Mate suffix
}

// Creates instance and filles it by constant part of Read name
ReadName::ReadName() : _headLen(Product::Title.length()), _len(0)
{
	_name = new char[Read::IsNameEmpty() ? _headLen : MaxLen];
	memcpy(_name, Product::Title.c_str(), _headLen);
	if(!Read::IsNameEmpty())
		_headLen += sprintf(_name + _headLen, "%c%s", Read::NmDelimiter, Chrom::Abbr);
	_len = _headLen;
}

// Copy constructor
ReadName::ReadName(const ReadName& hrName) : _headLen(hrName._headLen)
{
	_name = new char[Read::IsNameEmpty() ? Product::Title.length() : MaxLen];
	memcpy(_name, hrName._name, _headLen);
}

// Sets current chrom's mark
//	@cID: chrom's ID
void ReadName::SetChrom(chrid cID)
{
	if(!Read::IsNameEmpty())
		_chrLen = _headLen + 
			sprintf(_name + _headLen, "%s%c", Chrom::Mark(cID).c_str(), Read::GetDelimiter());
}

/************************ struct ReadName: end ************************/


/************************ class ExtOutFile ************************/

const char*	ExtOutFile::Mate[] = {"/1", "/2"};					// Array of Mate suffixes
const BYTE	ExtOutFile::MateLen = strlen(ExtOutFile::Mate[0]);	// The length of Mate suffix
bool ExtOutFile::Zipped;										// true if filed should be zippped

// Adds byte to the current position in the line write buffer,
//	adds delimiter after byte and increases current position.
//	@ch: value to be set
//	@addDelim: if true then adds delimiter and increases current position
void ExtOutFile::LineAddChar(char ch, bool addDelim)
{ 
	_lineBuff[_lineBuffOffset++] = ch;
	if(addDelim)	_lineBuff[_lineBuffOffset++] = _delim;
}

// Copies block of chars before the current position in the line write buffer.
//	@src:  pointer to the block of chars
//	@num: number of chars
//	@addDelim: if true then adds delimiter before string and decreases current position
void ExtOutFile::LineAddCharsBack(const char* src, size_t num, bool addDelim)
{
	if(addDelim)	_lineBuff[--_lineBuffOffset] = _delim;
	_lineBuffOffset -= num;
	memcpy(_lineBuff + _lineBuffOffset, src, num);
}

// Copies default Read name and pos extention after current position in the line write buffer,
//	and decreases current position.
//	@mate: mate number for PE Reads, or 0 for SE Read
void ExtOutFile::LineAddReadName(BYTE mate)
{
	LineAddStr(_rName.Name(), _rName.Length(), !mate);
	if(mate)	LineAddChars(Mate[mate-1], MateLen, true);
}

// Copies default Read name and pos extention  before current position in the line write buffer,
//	and decreases current position.
//	@mate: mate number for PE Reads, or 0 for SE Read
//	@addDelim: if true then adds delimiter after Read Name and decreases current position
void ExtOutFile::LineAddReadNameBack(BYTE mate, bool addDelim )
{
	if(mate)
		LineAddCharsBack(Mate[mate-1], MateLen, addDelim), addDelim = false;
	LineAddStrBack(_rName.Name(), _rName.Length(), addDelim);
}

/************************ class ExtOutFile: end ************************/


Seq::sMode	Seq::mode;		// sequence mode
ULLONG	Seq::maxFragCnt;	// up limit of saved fragments


/************************ class BedROutFile ************************/

// Creates new instance for writing and initializes line write buffer.
//	@fName: file name without extention
//	@rName: Read's name
//	@commLine: command line to add as a comment in the first line
BedROutFile::BedROutFile(const string& fName, const ReadName& rName, const string& commLine)
	: ExtOutFile(fName + FT::RealExt(FT::BED, ExtOutFile::Zipped), rName)
{
	CommLineToIOBuff(commLine);
	if(OutFiles::MultiThread)		Write();
	SetLineBuff(
		Chrom::MaxAbbrNameLength +		// length of chrom name
		ReadName::MaxLength() + 		// length of Read name
		2 * CHRLEN_CAPAC +				// start + stop positions
		OutFiles::MapQual.length() +	// score
		1 + 2 + 6);						// strand + HASH + BLANK + 5 TABs + EOL
}

// Sets treated chrom's name to line write buffer
void BedROutFile::SetChrom(chrid cID)
{
	LineSetOffset(0);
	_offset = LineAddStr(Chrom::AbbrName(cID));
}

// Adds Read to the line's write buffer.
//	@pos: valid Read's start position
//	@reverse: if true then set reverse strand, otherwise set forward
//	@mate: mate number for PE Reads, or 0 for SE Read
void BedROutFile::AddRead(chrlen pos, bool reverse, BYTE mate)
{
	LineAddInts(pos, pos + Read::Len);			// start, end
	LineAddReadName(mate);						// Read name
	LineAddStr(OutFiles::MapQual);				// score
	LineAddChar(Read::Strands[int(reverse)]);	// strand
	LineToIOBuff(_offset);
}

/************************ class BedROutFile: end ************************/

/************************ class FqOutFile ************************/

rowlen FqOutFile::ReadStartPos = 0;	// Read field constant start position 

// Creates new instance for writing
//	@fName: file name without extention
//	@rName: Read's name
//	@qPatt: external Read quality pattern
FqOutFile::FqOutFile(const string& fName, const ReadName& rName, ReadQualPattern& qPatt)
	: ExtOutFile(fName + FT::RealExt(FT::FQ, ExtOutFile::Zipped), rName, EOL) 
{
	if(!ReadStartPos)								// if not initialized yet
		ReadStartPos = ReadName::MaxLength() + 2;	// Read name + AT + EOL
	//== set line write buffer and writing start position
	rowlen pos = ReadStartPos + Read::Len + 3;
	SetLineBuff(pos + Read::Len);
	LineSetOffset(pos - 2);
	LineAddChar(PLUS, true);
	LineFillReadPatt(pos, qPatt);
}

// Forms Read from fragment and adds it to the file.
//	@read: valid read
//	@reverse: if true then complement added read 
//	@mate: mate number for PE Reads, or 0 for SE Read
void FqOutFile::AddRead(const char* read, bool reverse, BYTE mate)
{
	LineSetOffset(ReadStartPos);
	if(reverse)		Read::CopyComplement(LineCurrPosBuf(), read);
	else			LineCopyChars(read, Read::Len);
	LineAddReadNameBack(mate);
	LineAddCharBack(AT);
	LineBackToBuffer();
}

// Prints output file names separated by comma
//void FqOutFile::PrintNames() const
//{
//	//cout << _file.FileName();
//	//if(_filePos)
//	//	cout << SepCm << _filePos->FileName() << SepCm << _fileNeg->FileName();
//}

/************************ class FqOutFile: end ************************/

/************************ class SamOutFile ************************/

#define	TAG_LEN	2
/*
SAM format:
Col	Field	Type	Regexp/Range Brief		description
----------------------------------------------------------------------------
1	QNAME	String	[!-?A-~]f1,255g			Query template NAME
2	FLAG	Int		[0,216-1]				bitwise FLAG
3	RNAME	String	\*|[!-()+-<>-~][!-~]*	Reference sequence NAME
4	POS		Int		[0,231-1]				1-based leftmost mapping POSition
5	MAPQ	Int		[0,28-1]				MAPping Quality
6	CIGAR	String	\*|([0-9]+[MIDNSHPX=])+	CIGAR string
7	RNEXT	String	\*|=|[!-()+-<>-~][!-~]*	Ref. name of the mate/next read
8	PNEXT	Int		[0,231-1]				Position of the mate/next read
9	TLEN	Int		[-231+1,231-1]			observed Template LENgth
10	SEQ		String	\*|[A-Za-z=.]+			segment SEQuence
11	QUAL	String	[!-~]+					ASCII of Phred-scaled base QUALity+33
----------------------------------------------------------------------------

bitwise FLAG
Bit		Description
--------------------------------------------------------
1 0x1	template having multiple segments in sequencing
2 0x2	each segment properly aligned according to the aligner
3 0x4	segment unmapped
4 0x8	next segment in the template unmapped
5 0x10	SEQ being reverse complemented
6 0x20	SEQ of the next segment in the template being reversed
7 0x40	the rst segment in the template
8 0x80	the last segment in the template
9 0x100	secondary alignment
0 0x200	not passing quality controls
1 0x400	PCR or optical duplicate
2 0x800	supplementary alignment
--------------------------------------------------------

 0	   00000	SE +
16	   10000	SE -
99	01100011	PE +
147	10010011	PE -

*/

rowlen SamOutFile::ReadStartPos = 0;
string SamOutFile::Fields5_6;		// combined value from 5 to 6 field: defined in constructor
const string SamOutFile::Fields7_9 = "*\t0\t0";	// combined value from 7 to 9 field for SE mode
const char SamOutFile::RNEXT = '=';				// Ref. name of the mate/next read
const char* SamOutFile::FLAG[2];

// Creates new instance for writing, initializes line write buffer writes header.
//	@fName: file name without extention
//	@rName: Read's name
//	@cmLine: command line
//	@cSizes: chrom sizes
//	@qPatt: external Read quality pattern
SamOutFile::SamOutFile(const string& fName, const ReadName& rName,
	const string& cmLine, const ChromSizes& cSizes, ReadQualPattern& qPatt)
	: ExtOutFile(fName + FT::RealExt(FT::SAM, ExtOutFile::Zipped), rName)
{
	if(Seq::IsPE())	{ FLAG[0] = "99";	FLAG[1] = "147"; }
	else			{ FLAG[0] = "0";	FLAG[1] = "16";	 }
	
	//== preset Fields5_6
	ostringstream ss;
	ss << OutFiles::MapQual << TAB << Read::Len << 'M';	// MAPping Quality + CIGAR: Read length
	Fields5_6 = ss.str();

	//== write header
	StrToIOBuff("@HD\tVN:1.0\tSO:unsorted");
	for(ChromSizes::cIter it=cSizes.cBegin(); it!=cSizes.cEnd(); it++) {
		ss.str(strEmpty);	// clear stream
		ss << it->second.Size();
		StrToIOBuff("@SQ\tSN:" + Chrom::AbbrName(CID(it)) + "\tLN:" + ss.str());
	}
	StrToIOBuff("@PG\tID:" + Product::Title + "\tPN:" + Product::Title + 
		"\tVN:" + Product::Version + "\tCL:" + cmLine);
	if(OutFiles::MultiThread)		Write();

	//== maximal length of write line buffer without Read & Quality fields, with delimiters
	ReadStartPos =				
		ReadName::MaxLength() +		// QNAME: Read name
		3 +							// FLAG: bitwise FLAG
		Chrom::MaxAbbrNameLength +	// RNAME: AbbrChromName
		CHRLEN_CAPAC +				// POS: 1-based start pos
		2 +							// MAPQ: MAPping Quality
		4 +							// CIGAR: Read length in 3 digits + letter 'M' or '='
		1 +							// RNEXT: Ref. name of the mate/next read: SE '*', PE '='
		CHRLEN_CAPAC +				// PNEXT: Position of the mate/next read
		3 + 						// TLEN: observed Template LENgth
		10 + 1 + 1 +				// number of TABs + EOL  + one for safety
		Read::Len;					// only at the time of line buffer setting
	
	//== set line buffer
	SetLineBuff		(ReadStartPos + Read::Len);
	LineFillReadPatt(ReadStartPos, qPatt);
	LineSetChar		(--ReadStartPos, TAB);
	ReadStartPos -= Read::Len;		// restore true value
	
	//=== preset SE pattern
	if( !Seq::IsPE() ) {
		LineSetOffset(ReadStartPos - Fields7_9.length() - Fields5_6.length() - 2);
		LineAddStr(Fields5_6);
		LineAddStr(Fields7_9);
	}
}

// Adds Read with prepared to the line's write buffer.
//	@mate: suffix of PE Read's name; NULL for SE Read
//	@fields7_9: prepared 7-9 fields (RNEXT,PNEXT,TLEN)
//	@len7_9: lengths of 7-9 fields string
//	@read: valid Read
//	@flag: FLAG field value
//	@mate: mate number for PE Reads, or 0 for SE Read
void SamOutFile::AddRead(BYTE mate, const char* fields7_9,
	int len7_9, const char* read, const char* flag, chrlen pos)
{
	LineSetOffset(ReadStartPos);
	LineCopyChars(read, Read::Len);		// SEQ: Read
	LineAddStrBack(fields7_9, len7_9);	// RNEXT + PNEXT + TLEN
	LineAddStrBack(Fields5_6);			// MAPQ + CIGAR
	LineAddStrBack(_POS, sprintf(_POS, "%s\t%s\t%d", flag, _cName.c_str(), pos));	// FLAG,RNAME,POS
	LineAddReadNameBack(mate);			// QNAME: Read name

	LineBackToBuffer();
}

// Adds two mate Reads to the line's write buffer.
//	@read1: valid first mate Read
//	@read2: valid second mate Read
//	@pos1: valid first mate Read's start position
//	@pos2: valid second mate Read's start position
//	@fLen: fragment's length
void SamOutFile::AddTwoReads(const char* read1, const char* read2, chrlen pos1, chrlen pos2, int fLen)
{
	AddRead(1, _POS, AddFields7_9(++pos2, fLen), read1, FLAG[0], ++pos1);
	AddRead(2, _POS, AddFields7_9(pos1, -fLen),  read2, FLAG[1], pos2);
}

/************************ class SamOutFile: end ************************/

/************************ class Coverage ************************/

// Adds fragment to coverage
//	@pos: frag's position
//	@len: frag's length
void Coverage::AddFrag(chrlen pos, fraglen len)
{
	iter it2, it1 = find(pos);	// 'start', 'end' points iterators
	fraglen val = 0;			// 'end' point value
	bool newEnd;				// true if 'end' point is new 

	// *** set 'start' point
	if(it1 == end()) {								// 'start' point is new
		it2 = it1 = insert(pair<chrlen,fraglen>(pos, 1)).first;
		if(it1 != begin())							// 'start' point is not the first one at all
			it1->second += (--it2)->second;			// correct val by prev point; keep it1 unchanged
	}
	else if((--(it2=it1))->second == ++it1->second)	// 'start' point exists; --it2 may be end(), never mind
		erase(it1), it1 = it2;						// remove 'start' point as duplicate
			
	// *** set 'end' point
	it2 = find(pos+len);
	if(newEnd = it2 == end()) {			// 'end' point is new
		iter it = it2 = insert(pair<chrlen,fraglen>(pos+len, 0)).first;
		val = ++it == end() ? 1 :		// 'end' point is the last one at all
			(--(it = it2))->second;		// grab val by prev point; keep it2 unchanged
	}
		
	// *** correct range between 'start' and 'end', set 'end' point value
	for(it1++; it1!=it2; it1++)			// correct values within range
		val = ++it1->second;
	if((--it1)->second == it2->second)	// is the last added point a duplicate?
		erase(it2);						// remove duplicate point
	else if(newEnd)
		it2->second = --val;			// set new 'end' point value
}

/************************ class Coverage: end ************************/

/************************ class Coverages ************************/

// Adds SE fragment to coverage
//	@pos: frag's position
//	@len: frag's length
//	@reverse: true if read is reversed (neg strand)
void Coverages::AddFrag(chrlen pos, fraglen len, bool reverse)
{
	_covers[Count-1]->AddFrag(pos, len);					// common coverage
	if(_covers[0])	_covers[reverse]->AddFrag(pos, len);	// strand coverage
}

/************************ Coverages: end ************************/

/************************ class WigOutFile ************************/

// Initializes line write buffer, adds command and definition lines
//	@cl: command line to add as a comment in the first line
//	@strand: string denoted strand, or empty string
void WigOutFiles::WigOutFile::Init(const string& fName, const string& cl, const string& strand)
{
	SetLineBuff(Chrom::MaxAbbrNameLength + 3 * CHRLEN_CAPAC);
	CommLineToIOBuff(cl);
	StrToIOBuff("track type=bedGraph name=\"" + 
		fName + FT::RealExt(FT::WIG, ExtOutFile::Zipped) + 
		"\" description=\"" + Product::Title + " actual coverage" +
		(strand.length() ? ": " + strand + " strand": "") +
		"\" color=50,130,190 autoScale=on");
}

// Creates new instance for writing and initializes line write buffer.
//	@fName: file name without extention
//	@cmLine: command line to add as a comment in the first file line
//	@cSizes: chrom sizes
WigOutFiles::WigOutFile::WigOutFile(
	const string& fName, const string& cmLine, const ChromSizesExt& cSizes)
	: TxtOutFile(fName + FT::RealExt(FT::WIG, ExtOutFile::Zipped), TAB)
{
	Init(fName, cmLine);
	for(ChromSizes::cIter it=cSizes.cBegin(); it!=cSizes.cEnd(); it++)
		if( cSizes.IsTreated(it) )
			AddElem(CID(it), CoverProxy());
}

// Creates new strand-separated instance for writing and initializes line write buffer.
//	@fName: file name without extention
//	@cmLine: command line to add as a comment in the first line
//	@strand: string denoted strand, or empty string
//	@wig: pattern replacing basic Chrom container
WigOutFiles::WigOutFile::WigOutFile(
	const string& fName, const string& cmLine, const string& strand, const WigOutFile& wig)
	: TxtOutFile(fName + FT::RealExt(FT::WIG, ExtOutFile::Zipped), TAB)
{
	Init(fName, cmLine, strand);
	Assign(wig);
}

WigOutFiles::WigOutFile::~WigOutFile()
{
	// save unsaved coverages
	for(Iter it=Begin(); it!=End(); it++)
		if(it->second.IsUnsaved())
			WriteCover(CID(it), it->second);
}

// Writes coverage for given chrom to file
void WigOutFiles::WigOutFile::WriteCover(chrid cID, CoverProxy& proxy)
{
	if(!proxy.Cover->empty()) {
		LineSetOffset(0);
		const rowlen offset = LineAddStr(Chrom::AbbrName(cID));
		Coverage::citer it = proxy.Cover->begin();
		pair<chrlen,covr> pt = *it;

		for(++it; it!=proxy.Cover->end(); pt=*it++)
			if(pt.second) {
				LineAddInts(pt.first, it->first, pt.second, false);	// start, end, coverage
				LineToIOBuff(offset);
			}
		Write();
	}
	proxy.RemoveCover();
	proxy.Saved = true;
}

// Close coverage for given chrom for increment.
//	@cID: chrom ID
//	@i:	index of wig file
void WigOutFiles::WigOutFile::CloseCover(chrid cID, BYTE i)
{
	// Coverages are created and filled in different threads independently.
	// To save them in chrom sorted order the total pool is examined 
	// each time the next coverage is completed.
	// If previous coveraged are filled without gups, they are recorded and removed from the pool.
	Iter itThis = GetIter(cID);
	CoverProxy& proxy = itThis->second;
	Mutex::Lock(Mutex::eType(Mutex::OTHER1+i));	// diferent mutexes for different WIG files
	proxy.Closed = true;
	bool save = true;

	for(Iter it=Begin(); it!=itThis; it++) {
		CoverProxy& proxy0 = it->second;
		if(!proxy0.Saved)
			if(save = proxy0.Closed)
				WriteCover(CID(it), proxy0);
			else	break;
	}
	if(save)
		WriteCover(cID, proxy);
	Mutex::Unlock(Mutex::eType(Mutex::OTHER1+i));
}

/************************ class WigOutFile: end ************************/

/************************ class WigOutFiles ************************/

bool WigOutFiles::IsStrands;		// true if wigs with different strands should be generated

// Creates new instance for writing and initializes line write buffer.
//	@fName: file name without extention
//	@cmLine: command line to add as a comment in the first line
//	@cSizes: chrom sizes
WigOutFiles::WigOutFiles(const string& fName, const string& cmLine, const ChromSizesExt& cSizes)
{
	memset(_files, 0, Count * sizeof(WigOutFile*));

	_files[Count-1] = new WigOutFile(fName, cmLine, cSizes);
	if(IsStrands) {
		_files[0] = new WigOutFile(fName + "_pos", cmLine, "positive", *_files[Count-1]);
		_files[1] = new WigOutFile(fName + "_neg", cmLine, "negative", *_files[Count-1]);
	}
}

WigOutFiles::~WigOutFiles()
{
	for(BYTE i=!IsStrands*(Count-1); i<Count; i++)	// last item or all ones
		delete _files[i];
}

// Starts accumalating coverage for given chrom
//	@cID: chrom
//	@covrs: extern Coverage pointers
void WigOutFiles::StartChrom(chrid cID, Coverages& covrs)
{
	for(BYTE i=!IsStrands*(Count-1); i<Count; i++)	// last item or all ones
		covrs.SetCover(_files[i]->At(cID).CreateCover(), i);
}

// Stops accumalating coverage for given chrom
// record this one if all previous ones are already recorded
void WigOutFiles::StopChrom(chrid cID)
{
	for(BYTE i=!IsStrands*(Count-1); i<Count; i++)	// last item or all ones
		_files[i]->CloseCover(cID, i);
}

// Prints output file names separated by comma
void WigOutFiles::PrintNames() const
{
	cout << _files[Count-1]->FileName();
	if(IsStrands)
		for(int i=0; i<Count-1; i++)
			cout << SepCm << _files[i]->FileName();
}


/************************ class WigOutFiles: end ************************/

/************************ class OutFile ************************/

OutFiles::OutFile::tpAddRead	OutFiles::OutFile::pAddRead = &OutFiles::OutFile::AddReadSE;
OutFiles::OutFile::tpAddRInfo	OutFiles::OutFile::pAddReadInfo;

ULLONG	OutFiles::OutFile::rCnt = 0;		// total Read counter within instance
float	OutFiles::OutFile::StrandErrProb;	// the probability of strand error

// Creates and initializes new instance for writing.
//	@fName: common file name without extention
//	@cSizes: chrom sizes
//	@rqPatt: external Read quality pattern
//	@cmLine: command line
OutFiles::OutFile::OutFile(const string& fName, const ChromSizesExt& cSizes,
	ReadQualPattern& rqPatt, const string& cmLine) :
	_bedFile(NULL), _samFile(NULL), _wigFile(NULL), _primer(true)
{
	_fqFile1 = _fqFile2 = NULL;

	if(HasFormat(ofFQ)) {
		_fqFile1 = new FqOutFile(Seq::IsPE() ? fName + "_1" : fName, _rName, rqPatt);
		if(Seq::IsPE())
			_fqFile2 = new FqOutFile(fName + "_2", _rName, rqPatt);
	}
	if(HasFormat(ofBED))
		_bedFile = new BedROutFile(fName, _rName, cmLine);
	if(HasFormat(ofSAM))
		_samFile = new SamOutFile (fName, _rName, cmLine, cSizes, rqPatt);
	if(HasFormat(ofWIG))
		_wigFile = new WigOutFiles(fName, cmLine, cSizes);
}

// Clone constructor for multithreading
//	@file: original instance
//	@threadNumb: number of thread
OutFiles::OutFile::OutFile(const OutFile& file) : _primer(false)
{
	_fqFile1 = file._fqFile1 ? (FqOutFile*)new ExtOutFile(*file._fqFile1) : NULL;	// non own members
	_fqFile2 = file._fqFile2 ? (FqOutFile*)new ExtOutFile(*file._fqFile2) : NULL;	// non own members
	_bedFile = file._bedFile ? new BedROutFile(*file._bedFile, _rName) : NULL;	// own copy constructor
	_samFile = file._samFile ? new SamOutFile (*file._samFile, _rName) : NULL;	// own copy constructor
	_wigFile = file._wigFile;		// wigFile is common for all clones
}

OutFiles::OutFile::~OutFile()
{
	if(_fqFile1)	delete _fqFile1;
	if(_fqFile2)	delete _fqFile2;
	if(_bedFile)	delete _bedFile;
	if(_samFile)	delete _samFile;
	if(_wigFile
	&& _primer)		delete _wigFile;
}

// Start recording chrom
void OutFiles::OutFile::BeginWriteChrom(chrid cID)
{
	_rName.SetChrom(cID);
	if(_bedFile) _bedFile->SetChrom(cID);			// set chrom's name for writing.
	if(_samFile) _samFile->SetChrom(cID);			// set chrom's name for writing.
	if(_wigFile) _wigFile->StartChrom(cID, _covers);// set chrom's coverage as current
}

// Adds one SE Read
//	@seq: cutted reference chromosome
//	@pos: current fragment's position
//	@len: length of current fragment
//	@g: FG or BG; needs for strand error imitation
//	return:	1: fragment is out of range (end of chrom)
//			0: Read is added successfully
//			-1: N limit is exceeded
int OutFiles::OutFile::AddReadSE(const RefSeq& seq, chrlen pos, fraglen len, Gr::Type g)
{
	//if(RandomReverse && !primer)	// for amplified frag _reverse is the same
	//	_reverse = _rng.Boolean();	// true if read is reversed (neg strand)

	bool reverse = RandomReverse ? _rng.Boolean() : false;	// true if read is reversed (neg strand)
	if(_wigFile) {
		_covers.AddFrag(pos, len, reverse);
		if(HasWigOnly())	return 0;
	}
	chrlen rPos = reverse ? pos + len - Read::Len : pos;	// Read's position
	const char* read = seq.Read(rPos);
	int ret = Read::CheckNLimit(read);
	if(ret)		return ret;			// 1: NULL read, -1: N limit is exceeded

	(this->*pAddReadInfo)(rPos, 0);	// add additional info to the Read name

	if(RandomReverse && g==Gr::FG && _rng.Sample(OutFile::StrandErrProb) ) {
		//short diff = ftr->Centre() - rPos - (Read::Len>>1);
		//short halfDiffPeak = (200 - ftr->Length())>>1;
		//if(reverse)		diff += halfDiffPeak;
		//else			diff -= halfDiffPeak;

		//short diff = short( reverse ? (rPos + Read::Len - ftr->Start) : ftr->End - rPos);
		//if(!_rng.Sample(float(diff) / fLen)) {
			//if(diff/10 < _freq.Length()) _freq[diff/10]++;
			reverse = !reverse;			// imitate strand error
		//}
	}
	if(_fqFile1)	_fqFile1->AddRead(read, reverse);
	if(_bedFile)	_bedFile->AddRead(rPos, reverse);
	if(_samFile)	_samFile->AddRead(read, rPos, reverse);
	return 0;
}

// Adds two PE Reads
//	@seq: cutted reference chromosome
//	@pos: current fragment's position
//	@len: length of current fragment
//	@g: FG or BG; needs for strand error imitation; not used
//	return:	1: fragment is out of range (end of chrom)
//			0: Reads are added successfully
//			-1: N limit is exceeded
int OutFiles::OutFile::AddReadPE(const RefSeq& seq, chrlen pos, fraglen len, Gr::Type g)
{
	if(_wigFile) {
		_covers.AddFrag(pos, len);
		if(HasWigOnly())	return 0;
	}
	int ret;
	const char* read1 = seq.Read(pos);
	if(ret = Read::CheckNLimit(read1))	return ret;	// 1: NULL read, -1: N limit is exceeded
	chrlen pos2 = pos + len - Read::Len;
	const char* read2 = seq.Read(pos2);
	if(ret = Read::CheckNLimit(read2))	return ret;	// 1: NULL read, -1: N limit is exceeded

	(this->*pAddReadInfo)(pos, len);		// add additional info to the Read name

	if(_fqFile1)	_fqFile1->AddRead(read1, false, 1),
					_fqFile2->AddRead(read2, true,  2);
	if(_bedFile)	_bedFile->AddRead(pos, false, 1),
					_bedFile->AddRead(pos2, true, 2);
	if(_samFile)	_samFile->AddTwoReads(read1, read2, pos, pos2, len);
	return 0;
}

// Prints output file formats and sequencing mode
//	@signOut: output marker
//	@predicate: 'output' marker
void OutFiles::OutFile::PrintFormat(const char* signOut, const char* predicate) const
{
	if(HasFormat(ofFQ)) {
		cout << signOut << predicate << "sequence: " << _fqFile1->FileName();
		if(Seq::IsPE())		cout << SepCm << _fqFile2->FileName();
		cout << EOL;
	}
	if(HasFormat(ofBED, ofSAM)) {
		cout << signOut << predicate << "alignment: ";
		if(HasFormat(ofBED)) {
			cout << _bedFile->FileName();
			if(HasFormat(ofSAM))	cout << SepCm;
		}
		if(HasFormat(ofSAM))	cout << _samFile->FileName();
		cout << EOL;
	}
	if(HasFormat(ofWIG)) {
		cout << signOut << predicate << "coverage: ";
		_wigFile->PrintNames();
		cout << EOL;
	}
}

/************************ class OutFile: end ************************/

/************************ class OutFiles ************************/

string	OutFiles::MapQual;				// the mapping quality
bool	OutFiles::MultiThread;			// true if program is executed in a  multi thread
bool	OutFiles::RandomReverse = true;	// true if Read should be reversed randomly
OutFiles::oFormat OutFiles::Format;		// output formats

// Prints item title ("reads|fragments") accordingly file formats
void OutFiles::PrintItemTitle()
{
	const char* frags = "fragments";
	if(HasWigOnly())			cout << frags;
	else if(!HasFormat(ofWIG))	cout << FT::ItemTitle(FT::ABED);
	else cout << FT::ItemTitle(FT::ABED) << '/' << frags;
	cout << COLON;
}

// Initializes static members
//	@fFormat: types of output files
//	@mapQual: the mapping quality
//	@wigStrand: true if wigs with different strands should be generated
//	@strandErrProb: the probability of strand error
//	@zipped: true if output files should be zipped
void OutFiles::Init(oFormat fFormat, BYTE mapQual, bool wigStrand, float strandErrProb, bool zipped)
{
	Format = fFormat;
	MapQual = BSTR(mapQual);
	WigOutFiles::IsStrands = !Seq::IsPE() && wigStrand;
	ExtOutFile::Zipped = zipped;
	OutFile::Init(strandErrProb);
	ReadName::Init();
}

// Creates new instance for writing.
//	@fName: common file name without extention
//	@control: if true, then control ('input') is generated
//	@cSizes: chrom sizes, or NULL
//	@cFiles: chrom files
//	@rqPattFName: name of valid file with Read quality pattern, or NULL
//	@cmLine: command line
OutFiles::OutFiles(
	const string& fName, bool control, const ChromSizesExt& cSizes,
	const char* rqPattFName, const string& cmLine)
	: _freq(NULL), _gMode(GM::Test)
{
	_gMode = GM::Test;
	ReadQualPattern qPatt(rqPattFName);

	memset(_oFiles, 0, 2*sizeof(OutFile*));	// in case of an incomplete constructor calls
	_oFiles[0] = new OutFile(fName, cSizes, qPatt, cmLine);
	if(control)	_oFiles[1] = new OutFile(fName+"_input", cSizes, qPatt, cmLine);

	_rqTempl = qPatt.Templ();	// must be after OutFile constructor because it fills the pattern
	if(HasFormat(ofFREQ)) {
		_fFreqName = fName + ".freq"
#ifdef OS_Windows
			+".txt"
#endif
		;
		_freq = new FragFreq();
	}
}

// Clone constructor for multithreading.
//	@file: original instance
OutFiles::OutFiles(const OutFiles& file) : _freq(NULL), _gMode(file._gMode)
{
	_oFiles[0] = new OutFile(*file._oFiles[0]);
	_oFiles[1] = file._oFiles[1] ? new OutFile(*file._oFiles[1]) : NULL;
}
	
OutFiles::~OutFiles()
{
	for(int i=0; i<2; i++)
		if(_oFiles[i])	delete _oFiles[i];
	if(_freq) {
		ofstream fs;
		fs.open(_fFreqName.c_str());
		fs	<< "Fragment's distribution:: lognormal: sigma=" << DistrParams::lnSigma
			<< ", mean=" << DistrParams::lnMean
			<< "; size selection normal: ";
		if(DistrParams::IsSS())
			fs << "sigma=" << DistrParams::ssSigma << ", mean=" << DistrParams::lnMean << EOL;
		else
			fs << "OFF\n";
		_freq->Print(fs);
		fs.close();
		delete _freq;
	}
}

// Starts recording chrom
void OutFiles::BeginWriteChrom(chrid cID)
{
	_oFiles[0]->BeginWriteChrom(cID);
	if(_oFiles[1])	_oFiles[1]->BeginWriteChrom(cID);
}

// Stops recording chrom
void OutFiles::EndWriteChrom(chrid cID)
{
	_oFiles[0]->EndWriteChrom(cID);
	if(_oFiles[1])	_oFiles[1]->EndWriteChrom(cID);
}

// Prints output file formats and sequencing mode
//	@signOut: output marker
void OutFiles::PrintFormat	(const char* signOut) const
{
	const char* output = "Output ";

	_oFiles[0]->PrintFormat(signOut, output);
	if(HasFormat(ofFREQ))
		cout << signOut << output << "frequency distribution: " << _fFreqName << EOL;
	if(_oFiles[1])	cout << signOut << output << "control supplied\n";
}

// Prints Read quality settins
//	@signOut: output marker
void OutFiles::PrintReadQual (const char* signOut) const 
{
	if(HasWigOnly() || Format == ofFREQ)	return;
	cout << signOut << "Read quality" << SepDCl;
	if(HasFormat(ofSAM, ofFQ)) {
		cout << "pattern: ";
		if(_rqTempl)	{
			string templ(_rqTempl, Read::Len);
			cout << '\'' << templ << '\'';
		}
		else	Read::PrintSeqQuality();
		if(HasFormat(ofSAM, ofBED))	cout << SepSCl;
	}
	if(HasFormat(ofSAM, ofBED))		cout << "mapping = " << MapQual;
	cout << EOL;
}

/************************ class OutFiles: end ************************/
