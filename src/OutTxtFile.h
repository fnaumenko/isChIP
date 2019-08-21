/**********************************************************
OutTxtFile.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 03.04.2019
-------------------------
Provides output text files functionality
***********************************************************/
#pragma once
#include <fstream>		// Frag freguency ofstream
#include "Data.h"
//#include <time.h>		// if random number gen is in separate file


//#define RAND_STD			// rand() (Windows) or rand_r(int *seed) (Linux)
// Mersenne Twister by Agner Fog, 2008-11-16 http://www.agner.org/random/
//#define RAND_MT
// Xorshift by George Marsaglia http://en.wikipedia.org/wiki/Xorshift
#define RAND_XORSHIFT

#ifdef RAND_STD
	#define RAND_MAX_	RAND_MAX
#ifdef __unix__
	#define rand()	rand_r(&_seed)
#endif
#else
	#define RAND_MAX_	0xFFFFFFFF	//(65536.*65536.)	// 2^32
#endif

// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
// If this doesn't work then insert compiler-specific definitions here:
#if defined(__GNUC__) || (defined(_MSC_VER) && _MSC_VER >= 1600)
  // Compilers supporting C99 or C++0x have stdint.h defining these integer types
  #include <stdint.h>
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
  // 16 bit systems use long int for 32 bit integer.
  typedef   signed long int int32_t;
  typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
  // Older Microsoft compilers have their own definition
  typedef   signed __int32  int32_t;
  typedef unsigned __int32 uint32_t;
  typedef   signed __int64  int64_t;
  typedef unsigned __int64 uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#else
  // This works with most compilers
  typedef signed int          int32_t;
  typedef unsigned int       uint32_t;
  typedef long long           int64_t;
  typedef unsigned long long uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#endif

  // 'Random' encapsulates random number generator.
class Random
{
#ifdef RAND_MT
// Version of Mersenne Twister: MT11213A or MT19937
#if 0	// Constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else   // Constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
#endif

private:
	static int Seed;

#ifdef RAND_STD
	int _seed;
#elif defined RAND_MT
	int mti;					// Index into mt
	uint32_t y, mt[MERS_N];		// State vector
	
	// Generates 32 random bits
	uint32_t rand();
#elif defined RAND_XORSHIFT
	uint32_t x, y, z, w;

	// Generates 32 random bits
	uint32_t rand();
#endif

protected:
	// Generates random double number in the interval 0 <= x < 1
	inline double DRand() { return (double)rand() / RAND_MAX_; }

public:
	// Sets and returns real seed
	//	@seed: if 0, random seed
	static	int SetSeed(UINT seed);

	Random();

	// Returns random integer within interval [1, max]
	int	Range(int max);
	
	// Returns true with given likelihood
	//	@sample: probability of returning true; from 0.0. to 1.0
	bool Sample(float sample);

	// Returns random true ot false with probability 0.5
	inline bool Boolean() { return rand() & 0x1; }
};

// 'GM" defines generation modes
static struct GM
{
private:
	static const char* title[];		// title: printed member's name
public:
	enum Mode { Test /*test only*/, Control	/*test and control*/ };

	// Gets title: first title letter
	inline static const char* Title(Mode m)	{ return title[m]; }
} gm;

// 'DistrParams' keeps fragment initial lognormal and normal size selection (SS) distributions params
// It is defined here but not in Imitator.h becuase of using in 
static struct DistrParams
{
	static float lnMean;	// mean of initial lognormal distribution
	static float lnSigma;	// sigma of initial lognormal distribution
	static float ssMean;	// mean of size selection normal distribution,
							// or 0 if SS is off, or mean of lognormal distribution by default
	static int	 ssSigma;	// sigma of size selection normal distribution

	// Initializes lognorm fragment and size selection distribution values;
	//	@ln: expectation & standard deviation of frag normal distribution
	//	@ss: expectation & standard deviation of size sel normal distribution
	//	@isSSset: true if size sel is applied
	static void Init(const pairVal& ln, const pairVal& ss, bool isSSset) {
		lnMean = ln.first;
		lnSigma = ln.second;
		if(isSSset) {
			ssMean = ss.first == vUNDEF ? LnMode() : ss.first;
			ssSigma = ss.second;
		}
		else
			ssMean = 0;		// size sel OFF
		//ssMean = ss.first;
		//if(!(ssSigma = ss.second))	ssMean = 0;		// size sel OFF
		//else if(ssMean == vUNDEF)	ssMean = LnMode();
	}

	// Returns mean of lognormal distribtion
	static float LnMean() { return exp(lnMean + lnSigma*lnSigma/2); }

	// Returns mode of lognormal distribtion
	static float LnMode() { return exp(lnMean - lnSigma*lnSigma); }

	static inline bool IsSS() { return bool(ssMean); }

} distrParams;

// 'Seq' implements ChIP sequencing attributes
static struct Seq
{
public:
	// sequencing modes
	enum sMode { SE, PE, Undef };

private:
	static sMode	mode;		// sequencing mode
	static ULLONG	maxFragCnt;	// up limit of saved fragments

public:
	// Initializes data
	//	@smode: sequence mode
	//	@rLimit: maximum possible number of writed Reads
	static inline void Init(int smode, ULONG rLimit) {
		maxFragCnt = rLimit >> (mode = sMode(smode));	// reduce twice in case if PE
	}

	// Returns stated sequencing mode
	static inline sMode Mode()	{ return mode; }

	// Returns true if PE mode is stated
	static inline bool IsPE()	{ return mode; }

	// Returns maximum possible number of recorded frags
	static inline ULLONG FragsLimit() { return maxFragCnt; }

	// Returns maximum possible number of recorded Reads
	static inline float	ReadsLimit() { return float(maxFragCnt << mode); }	// increase twice for PE

	static void Print() {
		cout << "Sequencing: " << (IsPE() ? "paired" : "single") << "-end"
			 << SepSCl << FT::ItemTitle(FT::ABED) << " limit = " << ReadsLimit() << EOL;
	}

} seq;

// 'ReadQualPattern' sets the Read quality pattern according to the external external or in its absence
class ReadQualPattern
{
	char* _rqPatt;	// Read quality pattern buffer, or NULL if no external Read quality
	char* _rqTempl;	// pointer to Read quality within FQ|SAM external template,
					// or NULL if non of such files
public:
	// Creates Read quality pattern buffer and fills it by first valid line from file.
	//	@rqPattFName: name of valid file with a quality line
	ReadQualPattern(const char* rqPattFName);

	~ReadQualPattern() { if(_rqPatt) delete[] _rqPatt, _rqPatt = NULL; }

	// Fills external FQ|SAM template by Read quality pattern and remembers it
	//	@templ: pointer to external FQ|SAM Read quality template
	void Fill(char* templ);

	// Returns Read quality pattern buffer
	inline const char* Pattern() const { return _rqPatt; }

	// Returns pointer to Read quality within FQ|SAM external template
	inline const char* Templ() const { return _rqTempl; }
};

// 'ReadName' represents Read's name in output files.
struct ReadName
{
private:
	typedef void (ReadName::*tpAddRInfo)(chrlen, fraglen);
	static const char NmDelimiter;		// basic name delimiter

	static tpAddRInfo	pAddInfo;	// pointer to the 'Adds info to the name' method
	static ULLONG	rCnt;			// total Read counter
	static BYTE		len;			// Maximum length of Read name

	char*	_name;		// Read's name
	BYTE	_headLen;	// length of the constant head part containing program's title and chrom's title
	BYTE	_chrLen;	// length of the constant head plus current mark
	BYTE	_len;		// total length of Read's name

	inline ULLONG CountIncr() {	return MultiThread ? InterlockedIncrement(&rCnt) : ++rCnt; }
	void AddNumb(chrlen, fraglen);
	void AddPos	(chrlen pos, fraglen);
	void Add2Pos(chrlen pos, fraglen len);

public:
	static bool	MultiThread;		// true if program is executed in a multi thread

	static void Init();

	// Returns maximum length of Read's name
	inline static BYTE MaxLength() { return len; }

	// Creates instance and filles it by constant part of Read name
	ReadName();

	ReadName(const ReadName& rName);

	inline ~ReadName() { delete[] _name; }

	// Gets Read name
	inline char* Name() const { return _name; }

	// Gets length of Read name
	inline BYTE	Length() const { return _len; }

	// Sets current chrom's mark
	void SetChrom(chrid cID);

	// Adds info in Read name
	inline void AddInfo(chrlen pos, fraglen len = 0) { (this->*pAddInfo)(pos, len); }
};

// 'ExtOutFile' is extention of 'TxtOutFile'.
//	Implies additional methods for writing Read [extended] name,
//	and methods for writing data to line buffer back (left from current position)
class ExtOutFile : public TxtOutFile
{
public:
	static const BYTE MateLen;	// The length of Mate suffix
	static bool	Zipped;			// true if filed should be zippped

private:
	static const char* Mate[];	// Array of Mate suffixes

	const ReadName& _rName;	// Read's name

	// Copies block of chars before the current position in the line write buffer.
	//	@src:  pointer to the block of chars
	//	@num: number of chars
	//	adds delimiter before string and decreases current position.
	void LineAddCharsBack(const char* src, size_t num, bool addDelim);
	
protected:
	// Initializes Read name
	//	@headrName: conditionally constant head of Read name
	//	@tailrName: dynamic part of Read's name contained chroms number and Read's number|position,
	//	or NULL for empty name
	//void InitReadName(const ReadName* headrName, const ReadName* tailrName)
	//{ headRName = headrName; tailRName = tailrName; }

	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name with extention
	//	@rName: Read's name
	//	@trName: dynamic part of Read's name contained chroms number and Read's number|position,
	//	or NULL for empty name
	//	@commLine: command line to add as a comment in the first line
	inline ExtOutFile(const string& fName, const ReadName& rName, char delim=TAB)
		: _rName(rName), TxtOutFile(fName, delim) {}

	// Copy constructor for multithreading
	//	@rName: Read's name
	inline ExtOutFile(const ExtOutFile& file, const ReadName& rName) 
		: _rName(rName), TxtOutFile(file) {}

	// Returns a pointer to the line write buffer at current position.
	inline char* LineCurrPosBuf() const	{ return _lineBuff + _lineBuffOffset; }

	// Returns a pointer to the line write buffer at given position.
	//	@shift: given relative pointer position 
	//inline char* LineCurrPosBuf(rowlen shift) const	{ return _lineBuff + shift;	}

	// Moves current position of the line write buffer by shift.
	//inline void SlipOffset(int shift)	{ _buffLineOffset += shift; }

	// Decreases current position of the line write buffer by one.
	inline void LineDecreaseOffset() {	_lineBuffOffset--; }

	// Increases current position of the line write buffer by one.
	//inline void LineIncreaseOffset() { _buffLineOffset++; }

	// Sets byte in the line write buffer to the specified value.
	//	@offset: shift of buffer start position
	//	@ch: byte to be set
	inline void LineSetChar(rowlen offset, char ch)	{ _lineBuff[offset] = ch; }

	// Copies block of chars to the current position of the line write buffer.
	//	@src: pointer to the block of chars
	//	@len: number of chars to be copied
	inline void LineCopyChars(const char* src, size_t len) const {
		memcpy(_lineBuff + _lineBuffOffset, src, len);
	}

	// Adds byte to the current position in the line write buffer, adds delimiter after byte
	//	and increases current position.
	//	@ch: value to be set
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddChar(char ch, bool addDelim = false);

	// Adds byte before the current position of the line write buffer, adds delimiter before byte
	//	and decreases current position.
	//	@chr: value to be set
	//	@addDelim: if true then adds delimiter before and decreases current position
	inline void LineAddCharBack(char ch)	{ _lineBuff[--_lineBuffOffset] = ch; }

	// Copies the string before the current position in the line write buffer,
	//	adds delimiter after string and decreases current position.
	//	@str: string to be copied
	//	@addDelim: if true then adds delimiter after string and decreases current position
	inline void LineAddStrBack(const string& str, bool addDelim = true) {
		LineAddCharsBack(str.c_str(), str.length(), addDelim);
	}

	// Copies the string before current position in the line write buffer
	//	and decreases current position.
	//	@str: C string to be copied
	//	@len: length of string to be copied
	//	@addDelim: if true then adds delimiter after string and decreases current position
	inline void LineAddStrBack(const char* str, int len, bool addDelim = true) {
		LineAddCharsBack(str, len, addDelim);
	}

	// Copies default Read name and pos extention after current position in the line write buffer,
	//	and decreases current position.
	//	@mate: mate number for PE Reads, or 0 for SE Read
	void LineAddReadName(BYTE mate);

	// Copies default Read name and pos extention  before current position in the line write buffer,
	//	and decreases current position.
	//	@addDelim: if true then adds delimiter after Read Name and decreases current position
	void LineAddReadNameBack(bool addDelim = true)
	{ 
		LineAddStrBack(_rName.Name(), _rName.Length(), addDelim); 
	}

	// Adds last part of the line write buffer (from current position to the end)
	//	to the file write buffer.
	inline void LineBackToBuffer()
	{ RecordToIOBuff(_lineBuff + _lineBuffOffset, _lineBuffLen - _lineBuffOffset); }

	// Fills line by Read quality pattern
	//	@pos: start position
	//	@qPatt: Read quality pattern
	inline void LineFillReadPatt(rowlen pos, ReadQualPattern& qPatt)
	{ qPatt.Fill(_lineBuff + pos); }
};

// 'BedROutFile' implements methods for writing BaseItems file
class BedROutFile : public ExtOutFile
{
private:
	rowlen _offset;		// the length of the constant chrom name in the line buffer

public:
	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name without extention
	//	@rName: Read's name
	//	@commLine: command line to add as a comment in the first line
	BedROutFile(const string& fName, const ReadName& rName, const string& commLine);

	// Clone constructor for multithreading
	inline BedROutFile(const BedROutFile& file, const ReadName& rName) : ExtOutFile(file, rName) {}

	// Sets treated chrom's name to line write buffer
	void SetChrom(chrid cID);

	// Adds Read to the line's write buffer.
	//	@pos: valid Read's start position
	//	@reverse: if true then set reverse strand, otherwise set forward
	//	@mate: mate number for PE Reads, or 0 for SE Read
	void AddRead(chrlen pos, bool reverse, BYTE mate=0);
};

// 'FqOutFile' implements methods for writing FQ file
class FqOutFile : public ExtOutFile
// 'public' to allow implicit conversion in ExtOutFile(const FqOutFile&) invoke
{
private:
	static rowlen ReadStartPos;		// constant Read field start position

public:
	// Creates new instance for writing
	//	@fName: file name without extention
	//	@rName: Read's name
	//	@trName: dynamic part of Read's name contained chroms number and Read's number|position,
	//	or NULL for empty name
	//	@qPatt: external Read quality pattern
	FqOutFile(const string& fName, const ReadName& rName, ReadQualPattern& qPatt);

	// Forms Read from fragment and adds it to the file.
	//	@read: valid read
	//	@reverse: if true then complement added read 
	void AddRead(const char* read, bool reverse);

	// Prints output file names separated by comma
	//void PrintNames() const;

};

// 'SamOutFile' implements methods for writing SAM file
class SamOutFile : public ExtOutFile
{
private:
	static rowlen ReadStartPos;		// constant Read field start position
	static string Fields5_6;		// combined value from 5 to 6 field: defined in constructor
	static const string Fields7_9;	// combined value from 7 to 9 field for SE mode: predefined
	static const char RNEXT;		// Ref. name of the mate/next read
	static const char* FLAG[];		// FLAG value for SE: 00000=0 (+), 10000=16 (-)
									// FLAG value for PE: 01100011=99 (+), 10010011=147 (-)
	static const BYTE tagCLlen = 3;	// length of the header line tag 'CL:'

	char	_POS[3 * CHRLEN_CAPAC];	// POS field
	string	_cName;					// current chrom's name

	// Adds Read with prepared 7-9 fields in local buffer to the line's write buffer.
	//	@fields7_9: prepared 7-9 fields (RNEXT,PNEXT,TLEN)
	//	@len7_9: lengths of 7-9 fields string
	//	@read: valid Read
	//	@flag: FLAG field value
	//	@pos: valid start position of mate Read (PE) or Read (SE)
	void AddRead(const char* fields7_9, int len7_9, const char* read, const char* flag, chrlen pos);

	inline int AddFields7_9(chrlen pos, int fLen)
	{ return sprintf(_POS, "%c\t%d\t%d", RNEXT, pos, fLen); }

public:
	// Creates new instance for writing, initializes line write buffer writes header.
	//	@fName: file name without extention
	//	@rName: Read's name
	//	@hrName: conditionally constant part of the Read name
	//	@trName: dynamic part of Read's name contained chroms number and Read's number|position,
	//	or NULL for empty name
	//	@cmLine: command line
	//	@cSizes: chrom sizes
	//	@qPatt: external Read quality pattern
	SamOutFile(const string& fName, const ReadName& rName,
		const string& cmLine, const ChromSizes& cSizes, ReadQualPattern& qPatt);

	// Clone constructor for multithreading
	SamOutFile(const SamOutFile& file, const ReadName& rName)
		: ExtOutFile(file, rName) {}

	// Sets current treated chrom
	inline void SetChrom(chrid cID) { _cName = Chrom::AbbrName(cID); }

	// Adds Read to the line's write buffer.
	//	@read: valid Read
	//	@pos: valid Read's start position
	//	@reverse: if true then set reverse strand, otherwise set forward
	inline void AddRead(const char* read, chrlen pos, bool reverse) {
		AddRead(Fields7_9.c_str(), Fields7_9.length(), read, FLAG[reverse], ++pos);
	}

	// Adds two mate Reads to the line's write buffer.
	//	@read1: valid first mate Read
	//	@read2: valid second mate Read
	//	@pos1: valid first mate Read's start position
	//	@pos2: valid second mate Read's start position
	//	@fLen: fragment's length
	void AddTwoReads(const char* read1, const char* read2, chrlen pos1, chrlen pos2, int fLen);
};

#include <map>
typedef short	covr;	// type fragment coverage

// 'Coverage' represets chrom's fragment coverage data
//	and implements a single method for gradual filling (incrementing) coverage
class Coverage : public map<chrlen,covr>
{
	typedef map<chrlen,covr>::iterator iter;
public:
	typedef map<chrlen,covr>::const_iterator citer;

	// Adds fragment to coverage
	//	@pos: frag's position
	//	@len: frag's length
	void AddFrag(chrlen pos, fraglen len);

#ifdef DEBUG
	void Print() const
	{
		cout << "pos\tval\n";
		for(citer it=begin(); it!=end(); ++it)
			cout << it->first << TAB << it->second << EOL;
	}

	// Prints output in BedGraph format
	void Output() const
	{
		citer it=begin();
		chrlen start = it->first;
		covr val = it->second;

		cout << "start\tend\tval\n";
		for(++it; it!=end(); ++it) {
			if(val)
				cout << start << TAB << it->first << TAB << val << EOL;
			start = it->first;
			val = it->second;
		}
	}
#endif	// _DEBUG
};

// 'Coverages' represets 3 chrom's fragment coìerage data pointers,
//	keeping indepentently in WigOutFiles duplicates
class Coverages
{
private:
	static const BYTE Count = 3;

	Coverage* _covers[Count];

public:
	inline Coverages() { memset(_covers, 0, Count * sizeof(Coverage*));	}
	
	inline void SetCover(Coverage* cover, BYTE i) { _covers[i] = cover; }

	// Adds SE fragment to coverage
	//	@pos: frag's position
	//	@len: frag's length
	//	@reverse: true if read is reversed (neg strand)
	void AddFrag(chrlen pos, fraglen len, bool reverse);

	// Adds PE fragment to coverage
	//	@pos: frag's position
	//	@len: frag's length
	inline void AddFrag(chrlen pos, fraglen len) { _covers[Count-1]->AddFrag(pos, len); }
};

// 'WigOutFiles' implements methods for writing 1 or 3 WIG files
class WigOutFiles
{
	// Coverage representative (factory): 
	//	implements methods for Coverage creating, deleting and icnrement control
	struct CoverProxy
	{
		Coverage*	Cover;	// Coverage; initialized only while given chrom treatment
		bool	Closed;		// true if coverage is closed for increment
		bool	Saved;		// true if coverage is recorded

		inline CoverProxy() : Cover(NULL) { Closed = Saved = false; }

		inline Coverage* CreateCover() { return Cover = new Coverage(); }

		inline void RemoveCover() { if(Cover) delete Cover, Cover = NULL; }

		// Returns true if coverage is closed but unsaved
		inline bool IsUnsaved() { return !Saved && Closed; }
	};

	// 'WigOutFile' implements methods for writing BedGraph file
	class WigOutFile : public Chroms<CoverProxy>, public TxtOutFile
	{
		// Writes coverage for given chrom to file
		void WriteCover(chrid cID, CoverProxy& proxy);

		// Initializes line write buffer, adds command and definition lines
		//	@cl: command line to add as a comment in the first file line
		//	@strand: string denoted strand, or empty string
		void Init(const string& fName, const string& cl, const string& strand="");

	public:
		// Creates new instance for writing and initializes line write buffer.
		//	@fName: file name without extention
		//	@cmLine: command line to add as a comment in the first line
		//	@cSizes: chrom sizes
		WigOutFile(const string& fName, const string& cmLine, const ChromSizesExt& cSizes);

		// Creates new strand-separated instance for writing and initializes line write buffer.
		//	@fName: file name without extention
		//	@cmLine: command line to add as a comment in the first line
		//	@strand: string denoted strand, or empty string
		//	@wig: pattern replacing basic Chrom container
		WigOutFile(const string& fName, const string& cmLine, const string& strand, const WigOutFile& wig);

		~WigOutFile();

		// Creates and returns coverage for given chrom
		inline Coverage* OpenCover(chrid cID) { return At(cID).Data.CreateCover(); }

		// Close coverage for given chrom for increment.
		//	@cID: chrom ID
		//	@i:	index of wig file
		void CloseCover(chrid cID, BYTE i);
	};

	static const BYTE Count = 3;

	WigOutFile*	_files[Count];

public:
	static bool IsStrands;	// true if wigs with different strands should be generated

	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name with extention
	//	@cmLine: command line to add as a comment in the first line
	//	@cSizes: chrom sizes
	WigOutFiles(const string& fName, const string& cmLine, const ChromSizesExt& cSizes);

	~WigOutFiles();

	// Starts accumalating coverage for given chrom
	//	@cID: chrom
	//	@covrs: extern Coverage pointers
	void StartChrom(chrid cID, Coverages& covrs);

	// Stops accumalating coverage for given chrom
	void StopChrom(chrid cID);

	// Prints output file names separated by comma
	//	@signOut: output marker
	void PrintNames() const;
};

// 'Output' wraps test and control output files
class Output
{
public:
	// Output file formats
	enum oFormat {	// values are masks in oFormat variable
		ofFQ	= 0x01,
		ofBED	= 0x02,
		ofSAM	= 0x04,
		ofWIG	= 0x08,
		ofFREQ	= 0x10
	};

	static string	MapQual;	// the mapping quality
	static bool	RandomReverse;	// true if Read should be reversed randomly
private:
	static oFormat Format;		// output formats

	// 'OutFile' wraps output files: FQ, SAM, BED, WIG.
	class OutFile
	{
		typedef int	 (OutFile::*tpAddRead)(const RefSeq&, chrlen, fraglen, Gr::Type);
	
		// pointer to the 'add Read' method
		static tpAddRead pAddRead;
		static float	 StrandErrProb;		// the probability of strand error

		FqOutFile	* _fqFile1;	// mate1 or single output FQ
		FqOutFile	* _fqFile2;	// mate2 output FQ
		BedROutFile	* _bedFile;	// output BED
		SamOutFile	* _samFile;	// output SAM
		WigOutFiles	* _wigFile;	// output WIG
		Coverages	_covers;	// Coverage pointers for WIG; each own for each clone
		ReadName	_rName;		// Read's name; local for clone independence
		Random	_rng;			// random number generator; local for clone independence
		bool	_primer;		// true if file is primer (not clone); only for WigOutFile

		// Adds one SE Read
		//	@seq: cutted reference chromosome
		//	@pos: current fragment's position
		//	@len: length of current fragment
		//	@g: FG or BG; needs for strand error imitation
		//	return:	1: fragment is out of range (end of chrom)
		//			0: Read is added successfully
		//			-1: N limit is exceeded
		int AddReadSE (const RefSeq& seq, chrlen pos, fraglen len, Gr::Type g);
	
		// Adds two PE Reads
		//	@seq: cutted reference chromosome
		//	@pos: current fragment's position
		//	@len: length of current fragment
		//	@g: FG or BG; needs for strand error imitation; not used
		//	return:	1: fragment is out of range (end of chrom)
		//			0: Reads are added successfully
		//			-1: N limit is exceeded
		int AddReadPE (const RefSeq& seq, chrlen pos, fraglen len, Gr::Type g);

		// Empty (trial) method
		inline int NotAddRead (const RefSeq&, chrlen, fraglen, Gr::Type)	{ return 0; }
		//inline int NotAddRead (const RefSeq*, chrlen, fraglen, bool, const Featr*)	{ return 0; }

	public:
		// Initializes static members
		//	@singleThread: true if single thread is set
		//	@sErrProb: the probability of strand error
		inline static void Init(float sErrProb) { StrandErrProb = sErrProb; }

		// Sets sequense mode.
		//	@trial: if true, then set empty mode, otherwise current working mode
		inline static void SetSeqMode(bool trial)	{ 
			pAddRead = trial ? &Output::OutFile::NotAddRead :
				(Seq::IsPE() ? &Output::OutFile::AddReadPE : &Output::OutFile::AddReadSE);
		}

		// Creates and initializes new instance for writing.
		//	@fName: common file name without extention
		//	@cSizes: chrom sizes
		//	@rqPatt: external Read quality pattern
		//	@cmLine: command line
		OutFile(const string& fName, const ChromSizesExt& cSizes,
			ReadQualPattern& rqPatt, const string& cmLine);

		// Clone constructor for multithreading
		//	@oFile: original instance
		//	@threadNumb: number of thread
		OutFile(const OutFile& oFile);
	
		~OutFile();

		// Start recording chrom
		void BeginWriteChrom(chrid cID);

		// Stop recording chrom
		inline void EndWriteChrom(chrid cID) { if(_wigFile)	_wigFile->StopChrom(cID); }

		// Adds read(s) to output file
		//	@seq: cutted reference chromosome
		//	@pos: current fragment's position
		//	@len: current fragment's length
		//	@g: FG or BG; needs for strand error imitation
		//	return:	1: fragment is out of range (end chrom)
		//			0: Read(s) is(are) added, or nothing (trial)
		//			-2: N limit is exceeded; Read(s) is(are) not added
		inline int AddRead (const RefSeq& seq, chrlen pos, fraglen len, Gr::Type g) {
			return (this->*pAddRead)(seq, pos, len, g);
		}
		//inline int AddRead (const RefSeq* seq, chrlen pos, fraglen fLen, bool primer, const Featr* ftr) {
		//	return (this->*callAddRead[_mode])(seq, pos, fLen, primer, ftr);

		// Prints output file formats and sequencing mode
		//	@signOut: output marker
		//	@predicate: 'output' marker
		void PrintFormat(const char* signOut, const char* predicate) const;
	};

	// Returns true if format f is defined
	static inline bool HasFormat(oFormat f) { return Format & f; }
	// Returns true if formats f1 or f2 are defined
	static inline bool HasFormat(oFormat f1, oFormat f2) { return Format & (f1 | f2); }
	// Returns true if only WIG format is defined
	static inline bool HasWigOnly() { return Format == ofWIG; }

	OutFile*	_oFiles[2];	// test [0] and control [1] OutFile objects
	string		_fFreqName;	// file name of freq's distribution
	FragFreq*	_freq;		// frag's frequency
	const char*	_rqTempl;	// pointer to Read quality within FQ|SAM template,
							// or NULL if external Read quality pattern is not set;
							// needed only for PrintReadQual() (printing prog params)
	GM::Mode	_gMode;		// current generation mode,
							// corresponding to the index in the files and counters arrays
public:
	// Returns true if SAM type is assigned.
	static inline bool IsSamSet() { return HasFormat(ofSAM); }

	// Prints item title ("reads|fragments") accordingly file formats
	static void PrintItemTitle();

	// Initializes static members
	//	@fFormat: types of output files
	//	@mapQual: the mapping quality
	//	@wigStrand: true if wigs with different strands should be generated
	//	@strandErrProb: the probability of strand error
	//	@zipped: true if output files should be zipped
	static void Init(int fFormat, BYTE mapQual, bool wigStrand, float strandErrProb, bool zipped);

	// Sets sequense mode.
	//	@trial: if true, then set empty mode, otherwise current working mode
	inline static void SetSeqMode(bool trial)	{ OutFile::SetSeqMode(trial); }

	// Creates new instance for writing.
	//	@fName: common file name without extention
	//	@control: if true, then control ('input') is generated
	//	@cSizes: chrom sizes, or NULL
	//	@cFiles: chrom files
	//	@rqPattFName: name of valid file with Read quality pattern, or NULL
	//	@cmLine: command line
	Output(const string& fName, bool control, const ChromSizesExt& cSizes,
		const char* rqPattFName, const string& cmLine);

	// Clone constructor for multithreading.
	//	@file: original instance
	Output(const Output& file);
	
	~Output();

	// Set generation mode
	//	@testMode: if true, set Test mode, otherwhise Control mode
	inline void SetGMode(GM::Mode gm) { _gMode = gm; }

	// Starts recording chrom
	void BeginWriteChrom(chrid cID);

	// Stops recording chrom
	void EndWriteChrom(chrid cID);

	// Adds read(s) to output file
	//	@seq: cutted reference chromosome
	//	@pos: current fragment's position
	//	@len: length of current fragment
	//	@g: FG or BG; needs for strand error imitation
	//	return:	1: fragment is out of range (end chrom)
	//			0: Read(s) is(are) added, or nothing (trial)
	//			-1: N limit is exceeded; Read(s) is(are) not added
	int AddRead (const RefSeq& seq, chrlen pos, fraglen len, Gr::Type g) {
		if(_freq)	_freq->AddFrag(len);
		return _oFiles[_gMode]->AddRead(seq, pos, len, g);
	}

	// Prints output file formats and sequencing mode
	//	@signOut: output marker
	void PrintFormat	(const char* signOut) const;

	// Prints Read quality settins
	//	@signOut: output marker
	void PrintReadQual	(const char* signOut) const;
};
