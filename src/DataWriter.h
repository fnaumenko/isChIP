/**********************************************************
DataWriter.h
Provides bioinfo writers functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 11/23/2023
***********************************************************/
#pragma once

#include <memory>		// smart ptr
#include <string>
#include "ChromSeq.h"
#include "ChromSizesExt.h"
#include "Distrib.h"
#include "OrderedData.h"
#include "RandomGen.h"

// 'GM" defines generation modes
static struct GM
{
private:
	static const char* title[];		// title: printed member's name
public:
	enum class eMode {
		Test,	// test only
		Control	// test and control
	};

	// Gets title: first title letter
	static const char* Title(eMode m) { return title[int(m)]; }
} gm;


// 'SeqMode' represents sequencing modes
static class SeqMode
{
public:
	// sequencing modes
	enum eMode { SE, PE, Undef };	// using 'enum class' in this case is inconvenient

private:
	static eMode	mode;				// sequencing mode
	static ULLONG	maxFragCnt;			// up limit of saved fragments
	static const char* modeTitles[];	// printed modes's name

public:
	// Initializes data
	//	@smode: sequence mode
	//	@rLimit: maximum possible number of writed Reads
	static void Init(int smode, ULONG rLimit) {
		maxFragCnt = rLimit >> (mode = eMode(smode));	// reduce twice in case if PE
	}

	// Returns stated sequencing mode
	static eMode Mode() { return mode; }

	// Returns true if PE mode is stated
	static bool IsPE() { return mode; }

	// Returns maximum possible number of recorded frags
	static ULLONG FragsLimit() { return maxFragCnt; }

	// Returns maximum possible number of recorded Reads
	static float	ReadsLimit() { return float(maxFragCnt << mode); }	// increase twice for PE

	// Prints sequencing modes
	//	@signOut: output marker
	static void Print(const char* signOut);
} seq;

// 'ReadName' keeps and manages the qualified name of Read
//	'qualified' means included app, chrom, unique number, [position], [mate]
class ReadName
{
private:
	typedef void (ReadName::* tfAddRInfo)(const Region&);

	static tfAddRInfo	fAddInfo;	// pointer to the 'Adds info to the name' method
	static BYTE		len;			// Maximum length of Read name

	char* _name = NULL;		// Read's name
	BYTE _headLen;			// length of the constant head part containing program's title and chrom's title
	BYTE _headChrLen = 0;	// length of the constant head plus current chrom mark
	BYTE _len = 0;			// total length of Read's name
	ULLONG& _rCnt;			// external Read counter

	ULLONG CountIncr() { return InterlockedIncrement(&_rCnt); }

	// Adds to Read name its position as string
	//void AddPos(const string& s);

	// methods called by pointer from AddInfo(...)
	void AddNumb(const Region&);
	void AddPosSE(const Region& frag);
	void AddPosPE(const Region& frag);

public:
	static void Init();

	// Returns maximum length of Read's name
	static BYTE MaxLength() { return len; }

	// Initializes instance by constant part of Read name
	//	@rCnt: external Read counter
	ReadName(ULLONG& rCnt);

	// Copy constructor
	//ReadName(const ReadName& rName);

	~ReadName() { delete[] _name; }

	// Gets qualified Read's name
	char* Name() const { return _name; }

	// Gets length of Read name
	BYTE	Length() const { return _len; }

	// Sets external primary read counter
	void SetReadCounter(ULLONG& cnt) { _rCnt = cnt; }

	// Sets current chrom's mark
	void SetChrom(const string&& cMark);

	// Adds info into Read name
	// Calls AddNumb() or AddPos() or AddPosPE()
	void AddInfo(const Region& frag) { (this->*fAddInfo)(frag); }
};


// 'ReadWriter' is base class for data out files containing reads.
//	Implies additional methods for writing Read extended name,
//	and methods for writing data to line buffer back (left from current position)
class ReadWriter : public TxtWriter
{
	// 'ReadQualPattern' represents Read quality pattern storage and methods
	class ReadQualPattern
	{
		unique_ptr<char> _rqPatt;	// Read quality pattern buffer, or nullptr if no external Read quality
		readlen _defLen = 0;		// length of defined by user pattern

	public:
		// Creates Read quality pattern buffer and fills it by first valid line from file.
		//	@rqPattFName: name of valid file with a quality line
		ReadQualPattern(const char* rqPattFName);

		// Fills external FQ|SAM template by fixed Read quality pattern
		//	@templ: pointer to external FQ|SAM Read quality template
		void Fill(char* templ)	const { memcpy(templ, _rqPatt.get(), Read::FixedLen); }

		// Fills external FQ|SAM template by variable Read quality pattern
		//	@templ: pointer to external FQ|SAM Read quality template
		void Fill(char* templ, readlen rlen) const;

		// Returns Read quality pattern buffer to print
		const void Print() const;
	};

	static const char* Mate[];		// Array of Mate suffixes
	static string sReadConstLen;	// constant string "length=XX"
	static unique_ptr<ReadQualPattern> RqPattern;	// Read quality pattern

	using tfAddReadName = void(ReadWriter::*)(BYTE);
	// 'Add qualified Read name' methods: [0] - empty method, [1] - with mate extention
	static tfAddReadName fAddReadNames[2];

	const ReadName& _rName;			// Read's name

	void AddReadNameEmpty(BYTE) {}

	// Adds Read mate to name after current position in the line write buffer,
	//	and increases current position.
	//	@mate: mate number for PE Reads, or 0 for SE Read
	void LineAddReadMate(BYTE mate) { LineAddChars(Mate[mate - 1], MateLen); }

	// Copies block of chars before the current position in the line write buffer.
	//	@src:  pointer to the block of chars
	//	@len: number of chars
	void LineAddCharsBack(const char* src, size_t len);

protected:
	static const string ReadLenTitle;				// string " length="

	// Creates new instance for writing and initializes line write buffer.
	//	@ftype: type of file
	//	@fName: file name with extention
	//	@rName: Read's name
	//	@trName: dynamic part of Read's name contained chroms number and Read's number|position,
	//	or NULL for empty name
	//	@commLine: command line to add as a comment in the first line
	ReadWriter(FT::eType ftype, const string& fName, const ReadName& rName, char delim = TAB)
		: _rName(rName), TxtWriter(ftype, fName, delim) {}

	// Copy constructor for multithreading
	//	@rName: Read's name
	ReadWriter(const ReadWriter& file)
		: _rName(file._rName), TxtWriter(file) {}

	// Returns a pointer to the line write buffer at current position.
	char* LineCurrPosBuf() const { return _lineBuff + _lineBuffOffset; }

	// Sets byte in the line write buffer to the specified value.
	//	@offset: shift of buffer start position
	//	@ch: byte to be set
	void LineSetChar(reclen offset, char ch) { _lineBuff[offset] = ch; }

	// Adds byte before the current position of the line write buffer and decreases current position.
	//	@chr: value to be set
	void LineAddCharBack(char ch) { _lineBuff[--_lineBuffOffset] = ch; }

	// Copies the string before the current position in the line write buffer,
	//	adds delimiter after string and decreases current position.
	//	@str: string to be copied
	void LineAddStrBack(const string& str) { LineAddCharsBack(str.c_str(), str.length()); }

	// Copies qualified Read name and positions after current position in the line write buffer,
	//	and increases current position.
	//	@addDelim: if true then adds delimiter after string and increases current position
	void LineAddReadName(bool addDelim = true) { LineAddChars(_rName.Name(), _rName.Length(), addDelim); }

	// Copies qualified Read name and positions after current position in the line write buffer,
	//	and increases current position.
	//	@mate: mate number for PE Reads, or 0 for SE Read
	//	Invoked in BedOutFile.
	void LineAddReadName(BYTE mate);

	// Copies qualified Read name started with '@' and read variable legth after current position
	//	in the line write buffer, and increases current position
	//	@len: Read length
	//	Invoked in FqWriter.
	void LineAddReadVarName(readlen len);

	// Copies qualified Read name before current position in the line write buffer,
	//	adds delimiter after Read Nameand decreases current position.
	void LineAddReadNameBack() { LineAddCharsBack(_rName.Name(), _rName.Length()); }

	// Copies qualified Read name started with '@', positions and Read constant length before current position
	//	in the line write buffer, adds delimiter after Read Name and decreases current position.
	//	Invoked in FqWriter.
	void LineAddReadConstNameBack();

	// Adds last part of the line write buffer (from current position to the end) to the file write buffer.
	void LineBackToBuffer() { RecordToIOBuff(_lineBuff + _lineBuffOffset, _lineBuffLen - _lineBuffOffset); }

	// Fills line by Read const quality pattern from the specified position
	void LineFillReadConstPatt(reclen pos) const { RqPattern->Fill(_lineBuff + pos); }

	// Fills line by Read variable quality pattern from the current position and increases current position
	//	@rlen: Read's length
	void LineFillReadVarPatt(readlen rlen);

#ifdef _DEBUG
	const void PrintBuff(const char* pos, size_t len) const {
		for (int i = 0; i < len; i++)	cout << *(pos + i);	cout << LF;
	}
#endif

public:
	static bool	MultiThread;		// true if program is executed in a multi thread
	static const BYTE MateLen;		// The length of Mate suffix

	// Sets Read quality pattern and const Read length
	//	@rqPattFName: valid file name of the Read quality pattern
	static void SetReadQualityPatt(const char* rqPattFName) {
		RqPattern.reset(new ReadQualPattern(rqPattFName));
		if (!DistrParams::IsRVL())	sReadConstLen = ReadLenTitle + to_string(Read::FixedLen);
	}

	static void Init() { ReadName::Init(); }

	// Prints Read quality pattern
	static void PrintReadQualPatt() { RqPattern->Print(); }
};

// 'RBedWriter' implements methods for writing BED alignment file
class RBedWriter : public ReadWriter
{
	reclen _offset = 0;		// the length of the constant chrom name in the line buffer

public:
	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name without extention
	//	@rName: Read's name
	//	@commLine: command line to add as a comment in the first line, or nullptr
	RBedWriter(const string& fName, const ReadName& rName, const string* commLine = nullptr);

	// Clone constructor for multithreading
	RBedWriter(const RBedWriter& file) : ReadWriter(file) {}

	// Sets treated chrom's name to line write buffer
	//void SetChrom(chrid cID);
	void SetChrom(const string& chr);

	// Adds Read to the line's write buffer.
	//	@read: valid Read
	//	@reverse: if true then set reverse strand, otherwise set forward
	//	@mate: mate number for PE Reads, or 0 for SE Read
	void AddRead(const Read& read, bool reverse, BYTE mate = 0);
};

// 'FqWriter' implements methods for writing FQ file
class FqWriter : public ReadWriter
	// 'public' to allow implicit conversion in ReadWriter(const FqWriter&) invoke
{
	static reclen ReadStartPos;		// fixed Read field start position

	typedef void(FqWriter::* fAddRead)(const Read&, bool);
	// Current 'add read' method: with fixed or variable length 
	static fAddRead addRead;

	// Adds Read with fixed length to the line's write buffer.
	//	@read: valid Read
	//	@reverse: if true then add complemented read 
	void AddFLRead(const Read& read, bool reverse);

	// Adds Read with variable length
	void AddVLRead(const Read& read, bool reverse);

public:
	// initializes static members
	static void Init() { addRead = DistrParams::IsRVL() ? &FqWriter::AddVLRead : &FqWriter::AddFLRead; }

	// Creates new instance for writing
	//	@fName: file name without extention
	//	@rName: Read's name
	FqWriter(const string& fName, const ReadName& rName);

	// Clone constructor for multithreading
	FqWriter(const FqWriter& file) : ReadWriter(file) {}


	// Forms Read from fragment and adds it to the file.
	//	@read: valid Read
	//	@reverse: if true then add complemented read 
	void AddRead(const Read& read, bool reverse) { (this->*addRead)(read, reverse); }
};

// 'SamWriter' implements methods for writing SAM file
class SamWriter : public ReadWriter
{
	static reclen ReadStartPos;			// fixed Read field start position
	static string Fld_5_6;				// combined value from 5 to 6 field: initialised in constructor
	static string FLAG[];				// FLAG value for SE: 00000=0 (+), 10000=16 (-)
										// FLAG value for PE: 01100011=99 (+), 10010011=147 (-)
	const string Fld_7_9 = "*\t0\t0";	// combined value from 7 to 9 field for SE mode: predefined
	const BYTE tagCLlen = 3;			// length of the header line tag 'CL:'
	const char CIGAR_M = 'M';			// CIGAR marker

	typedef void(SamWriter::* tfAddRead)(const Read&, const string&, const string&);
	// Current 'add read' method: with fixed or variable length 
	static tfAddRead fAddRead;

	string	_cName;						// current chrom's name

	// Adds Read with fixed length to the line's write buffer.
	//	@read: valid Read
	//	@fld_7_9: prepared 7-9 fields (RNEXT,PNEXT,TLEN)
	//	@flag: FLAG field value to the line's write buffer.
	void AddFLRead(const Read& read, const string& fld_7_9, const string& flag);

	// Adds Read with variable length
	void AddVLRead(const Read& read, const string& fld_7_9, const string& flag);

public:
	// initializes static members
	static void Init() {
		fAddRead = DistrParams::IsRVL() ? &SamWriter::AddVLRead : &SamWriter::AddFLRead;
	}

	// Creates new instance for writing, initializes line write buffer writes header.
	//	@fName: file name without extention
	//	@rName: Read's name
	//	@cSizes: chrom sizes
	SamWriter(const string& fName, const ReadName& rName, const ChromSizes& cSizes);

	// Clone constructor for multithreading
	SamWriter(const SamWriter& file) : ReadWriter(file) {}

	// Sets current treated chrom
	void SetChrom(const string& chr) { _cName = chr; }

	// Adds Read to the line's write buffer.
	//	@read: valid Read
	//	@pos: Read's start position
	//	@len: Read's length
	//	@reverse: if true then add complemented read
	void AddRead(const Read& read, bool reverse) { (this->*fAddRead)(read, Fld_7_9, FLAG[reverse]); }

	// Adds two mate Reads to the line's write buffer.
	//	@read1: valid first mate Read
	//	@read2: valid second mate Read
	//	@fLen: fragment's length
	void AddTwoReads(const Read& read1, const Read& read2, int fLen);
};

// 'DataWriter' wraps writers, including test and control
class DataWriter
{
public:
	// DataWriter file formats
	enum class eFormat {	// values are masks in eFormat variable
		UNDEF = 0,
		FG = 1,
		BED = FG << 1,
		SAM = FG << 2,
		BGR = FG << 3,		// bedGraph
		FDENS = FG << 4,
		RDENS = FG << 5,
		FDIST = FG << 6,
		RDIST = FG << 7,
	};

private:
	static const string* commLine;	// command line to add at the first (commented) line of file
									// Is initizlies only in DataWriter constructor while the string is on the stack
	static int	Format;					// output formats as int
	static bool	inclReadName;			// true if Read name is included into output data
	static bool isStrand;				// true if coverage should be saved by strands as well
	static const char* entityTitles[];	// entity titles for printing
	static const BYTE ND = 2;			// number of distribution files

	// 'BioWriters' wraps writers for FQ, SAM, BED, BG, WIG
	class BioWriters
	{
		using OrderedCover = OrderedData<AccumCover, BedGrWriter>;
		using OrderedFreq = OrderedData<Freq, VarWigWriter>;

		typedef int	 (BioWriters::* tfAddRead)(const Region&, readlen, bool);
		static tfAddRead fAddRead;
		static float	 StrandErrProb;	// the probability of strand error

		mutable ULLONG	_rCnt = 0;		// total Read counter; managed by _rName, not used by clones
		const ChromSeq* _seq{};			// current reference sequence
		ReadName	 _rName{ _rCnt };	// Read's name; local for clone independence by setting different chroms
		unique_ptr<FqWriter>		_fqFile1{};		// FQ mate1 or single output
		unique_ptr<FqWriter>		_fqFile2{};		// FQ mate2 output 
		unique_ptr<RBedWriter>		_bedFile{};		// BED output
		unique_ptr<SamWriter>		_samFile{};		// SAM output
		unique_ptr<OrderedCover>	_bgFiles{};		// bedGraph output
		unique_ptr<OrderedFreq>		_fragWgFile{};	// frag density wig output
		unique_ptr<OrderedFreq>		_readWgFile{};	// read density wig output

		// Adds one SE Read
		//	@frag: added fragment
		//	@rLen: Read's length
		//	@reverse: if true then add complemented read
		//	return:	1: fragment is out of range (end of chrom)
		//			0: Read is added successfully
		//			-1: N limit is exceeded
		int AddReadSE(const Region& frag, readlen rlen, bool reverse);

		// Adds two PE Reads
		//	@frag: added fragment
		//	@rLen: Read's length
		//	@reverse: not used
		//	return:	1: fragment is out of range (end of chrom)
		//			0: Reads are added successfully
		//			-1: N limit is exceeded
		int AddReadPE(const Region& frag, readlen rlen, bool reverse);

		// Empty (trial) method
		//int AddReadEmpty (DataWriter*, const Region&, /*Gr::eType,*/ bool)	{ return 0; }
		int AddReadEmpty(const Region&, readlen, bool) { return 0; }

	public:
		// Initializes static members
		//	@singleThread: true if single thread is set
		//	@sErrProb: the probability of strand error
		static void Init(float sErrProb) { StrandErrProb = sErrProb; }

		// Sets sequense mode.
		//	@trial: if true, then set empty mode, otherwise current working mode
		static void SetSeqMode(bool trial) {
			fAddRead = trial ? &DataWriter::BioWriters::AddReadEmpty :
				(SeqMode::IsPE() ? &DataWriter::BioWriters::AddReadPE : &DataWriter::BioWriters::AddReadSE);
		}

		// Creates and initializes new instance for writing.
		//	@fName: common file name without extention
		//	@cSizes: chrom sizes
		BioWriters(const string& fName, const ChromSizesExt& cSizes);

		// Clone constructor for multithreading
		//	@param[in] primer: primer instance
		BioWriters(const BioWriters& primer);

		// Start recording chrom
		void BeginWriteChrom(const ChromSeq& seq);

		// Stop recording chrom
		void EndWriteChrom() const;

		// Adds read(s) to output file
		//	@frag: added fragment
		//	@rLen: Read's length
		///	@g: FG or BG; needs for strand error imitation
		//	@reverse: if true then add complemented read
		//	return:	1: fragment is out of range (end chrom)
		//			0: Read(s) is(are) added, or nothing (trial)
		//			-2: N limit is exceeded; Read(s) is(are) not added
		int AddRead(const Region& frag, readlen rlen, /*Gr::eType g,*/ bool reverse) {
			return (this->*fAddRead)(frag, rlen, reverse);
		}

		// Prints output file formats and sequencing mode
		//	@signOut: output marker
		//	@predicate: 'output' marker
		void PrintFormat(const char* signOut, const char* predicate) const;
	};

	// 'DistrWriters' manages two distribution: fragments and reads
	class DistrWriters
	{
		static const string fExt[];				// distrib file extentios
		static const char* entityAdjust[];		// entity title adjustment for printing

		const string _fName = strEmpty;			// common part of frag's/read's distribution file name
		Distrib* _dist[ND]{ nullptr,nullptr };	// distributions: fragments (0) and reads (1)
		function<void(fraglen flen)> _fAddFrag = [](fraglen) {};	// 'add fragment to distribution' function
		function<void(readlen rlen)> _fAddRead = [](readlen) {};	// 'add read to distribution' function

		// Returns file name for fragments (ind=0) and reads (ind=1) distributions
		const string FileName(BYTE ind) const {
			return _fName + fExt[ind]
#ifdef OS_Windows
				+ ".txt"
#endif
				;
		}

	public:
		DistrWriters(const string& fName, bool isFragDist, bool isReadDist);

		// Writes distributions to files and delete them
		~DistrWriters();

		// Adds frag/read length to statistics
		void AddFrag(fraglen flen, readlen rlen);

		// Prints output file formats
		//	@signOut: output marker
		//	@predicate: output common title
		void PrintFormat(const char* signOut, const char* predicate) const;
	};

	// Returns true if format f is defined
	static bool HasFormat(eFormat f) { return Format & int(f); }

	// Returns true if formats f1 or f2 or f3 are defined
	static bool HasFormat(eFormat f1, eFormat f2, eFormat f3 = eFormat::UNDEF)
	{
		return Format & (int(f1) | int(f2) | int(f3));
	}

	// Returns true if Read name is included into output data
	static bool InclReadName() { return inclReadName; }

	// Returns true if contiguous format shifted by i from the base f is defined
	static bool HasContigFormat(eFormat f, BYTE i) { return HasFormat(eFormat(BYTE(f) * 1 << i)); }

	// Returns true if both formats is defined
	static bool HasBothFormats(eFormat f1, eFormat f2) { return OnesCount(Format & (int(f1) | int(f2))) == 2; }

	unique_ptr<BioWriters> _writers[2];	// test [0] and control [1] writers
	shared_ptr<DistrWriters> _dists;	// frag's & read's distribution; common for duplicates
	Random	_rng;					// random generator; needed for Read variable length generation
	BYTE	_gMode;					// current generation mode as int,
									// corresponding to the index to call test/control files

	// Prints item title or count
	template <typename T>
	static void PrintItemsSummary(T t1, T t2) {
		bool isReadFormat = InclReadName() || HasFormat(eFormat::RDENS, eFormat::RDIST);
		if (isReadFormat)
			cout << t1;
		if (HasFormat(eFormat::BGR, eFormat::FDENS, eFormat::FDIST)) {
			if (isReadFormat)	cout << '/';
			cout << t2;
		}
	}

public:
	static string	MapQual;	// the mapping quality

	static const string* CommLine() { return commLine; }

	// Initializes static members
	//	@fFormat: types of output files
	//	@cmLine: command line to add as a comment in the first file line
	//	@mapQual: the mapping quality
	//	@bgStrand: true if bedGraphs with different strands should be generated
	//	@strandErrProb: the probability of strand error
	//	@zipped: true if output files should be zipped
	static void Init(int fFormat, BYTE mapQual, bool bgStrand, float strandErrProb, bool zipped);

	// Sets sequense mode.
	//	@trial: if true, then set empty mode, otherwise current working mode
	static void SetSeqMode(bool trial) { BioWriters::SetSeqMode(trial); }

	// Sets Read quality pattern by valid file name.
	static void SetReadQualPatt(const char* rqPattFName) { ReadWriter::SetReadQualityPatt(rqPattFName); }

	// Prints item title ("reads/fragments") according to output formats
	static void PrintItemTitle();

	// Prints item count (reads/fragments) according to output formats
	//	@fCnt: number of fragments
	static void PrintItemCount(ULLONG fCnt);

	// Creates new instance for writing.
	//	@fName: common file name without extention
	//	@control: if true, then control ('input') is generated
	//	@cmLine: command line to add as a comment in the first file line
	//	@cSizes: chrom sizes, or NULL
	DataWriter(const string& fName, bool control, const string& cmLine, const ChromSizesExt& cSizes);

	// Clone constructor for multithreading.
	//	@file: original instance
	DataWriter(const DataWriter& file);

	// Set generation mode
   //	@testMode: if true, set Test mode, otherwhise Control mode
	void SetGMode(GM::eMode gm) { _gMode = BYTE(gm); }

	// Starts recording chrom
	void BeginWriteChrom(const ChromSeq& seq);

	// Stops recording chrom
	void EndWriteChrom();

	// Adds read(s) to output file
	//	@pos: current fragment's position
	//	@flen: length of current fragment
	///	@g: FG or BG; needs for strand error imitation
	//	@reverse: if true then add complemented read
	//	return:	1: fragment is out of range (end chrom)
	//			0: Read(s) is(are) added, or nothing (trial)
	//			-1: N limit is exceeded; Read(s) is(are) not added
	int AddRead(chrlen pos, fraglen flen, /*Gr::eType g,*/ bool reverse);

	// Prints output file formats and sequencing mode
	//	@signOut: output marker
	void PrintFormat(const char* signOut) const;

	// Prints Read quality settins
	//	@signOut: output marker
	void PrintReadQual(const char* signOut) const;
};
