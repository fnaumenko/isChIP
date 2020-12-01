/**********************************************************
OutTxtFile.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 26.11.2020
-------------------------
Provides output text files functionality
***********************************************************/
#pragma once
#include <fstream>		// Frag freguency ofstream
#include <memory>		// smart ptr
#include "Data.h"
#include <map>	// Density
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
	inline static const char* Title(eMode m)	{ return title[int(m)]; }
} gm;

// 'Seq' represents seq mode and 
static class Seq
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
	static inline void Init(int smode, ULONG rLimit) {
		maxFragCnt = rLimit >> (mode = eMode(smode));	// reduce twice in case if PE
	}

	// Returns stated sequencing mode
	static inline eMode Mode()	{ return mode; }

	// Returns true if PE mode is stated
	static inline bool IsPE()	{ return mode; }

	// Returns maximum possible number of recorded frags
	static inline ULLONG FragsLimit() { return maxFragCnt; }

	// Returns maximum possible number of recorded Reads
	static inline float	ReadsLimit() { return float(maxFragCnt << mode); }	// increase twice for PE

	// Prints sequencing modes
	//	@signOut: output marker
	static void Print(const char* signOut) {
		cout << signOut << "Sequencing: " << modeTitles[IsPE()] << "-end"
			 << SepSCl << FT::ItemTitle(FT::eType::ABED) << " limit = " << ReadsLimit() << LF;
	}

} seq;

// 'ReadName' keeps and manages the qualified name of Read
//	'qualified' means the name included app name, chrom, unique number, [position], [mate]
class ReadName
{
private:
	typedef void (ReadName::*tpAddRInfo)(const Region&);

	static tpAddRInfo	pAddInfo;	// pointer to the 'Adds info to the name' method
	static ULLONG	rCnt;			// total Read counter
	static BYTE		len;			// Maximum length of Read name

	char*	_name;			// Read's name
	BYTE	_headLen;		// length of the constant head part containing program's title and chrom's title
	BYTE	_headChrLen;	// length of the constant head plus current chrom mark
	BYTE	_len;			// total length of Read's name

	inline ULLONG CountIncr() { return InterlockedIncrement(&rCnt); }

	// Adds to Read name its number followed be delimitor
	//	@shift: starting position in buffer
	void AddNumb(BYTE shift);
#ifdef OS_Windows
	// Adds to Read name its position
	//	@pos: added position
	void AddPos(chrlen pos) { _len = _headChrLen + PrintNumbToBuff(_name + _headChrLen, pos); }
#endif
	// methods called by pointer from AddInfo(...)
	inline void AddNumb(const Region&) { AddNumb(_headChrLen); }
	void AddPosSE(const Region& frag);
	void AddPosPE(const Region& frag);

public:
	static bool	MultiThread;		// true if program is executed in a multi thread

	static void Init();

	// Returns maximum length of Read's name
	inline static BYTE MaxLength() { return len; }

	// Initializes instance by constant part of Read name
	ReadName();

	ReadName(const ReadName& rName);

	inline ~ReadName() { delete[] _name; }

	// Gets qualified Read's name
	inline char* Name() const { return _name; }

	// Gets length of Read name
	inline BYTE	Length() const { return _len; }

	// Sets current chrom's mark
	void SetChrom(chrid cID);

	// Adds info in Read name
	// Calls AddNumb() or AddPos() or AddPosPE()
	inline void AddInfo(const Region& frag) { (this->*pAddInfo)(frag); }
};

// 'DataOutFile' is base class for data out files. In fact it's extention of 'TxtOutFile'.
//	Implies additional methods for writing Read [extended] name,
//	and methods for writing data to line buffer back (left from current position)
class DataOutFile : public TxtOutFile
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
		inline void Fill(char* templ)	const { memcpy(templ, _rqPatt.get(), Read::FixedLen); }

		// Fills external FQ|SAM template by variable Read quality pattern
		//	@templ: pointer to external FQ|SAM Read quality template
		void Fill(char* templ, readlen rlen) const;

		// Returns Read quality pattern buffer to print
		const void Print() const;
	};

public:
	static const BYTE MateLen;		// The length of Mate suffix
	static const string* CommLine;	// command line to add at the first (commented) line of file
									// Is initizlies only during Output constructor execution
									// while the string is on the stack
private:
	static const char* Mate[];		// Array of Mate suffixes
	static string sReadConstLen;					// constant string "length=XX"
	static unique_ptr<ReadQualPattern> RqPattern;	// Read quality pattern

	const ReadName& _rName;	// Read's name

	// Copies block of chars before the current position in the line write buffer.
	//	@src:  pointer to the block of chars
	//	@len: number of chars
	//	@addDelim: if true then adds delimiter before string and decreases current position.
	void LineAddCharsBack(const char* src, size_t len, bool addDelim);

protected:
	static const string ReadLenTitle;				// string " length="

	// Initializes Read name
	//	@headrName: conditionally constant head of Read name
	//	@tailrName: dynamic part of Read's name contained chroms number and Read's number|position,
	//	or NULL for empty name
	//void InitReadName(const ReadName* headrName, const ReadName* tailrName)
	//{ headRName = headrName; tailRName = tailrName; }

	// Creates new instance for writing and initializes line write buffer.
	//	@ftype: type of file
	//	@fName: file name with extention
	//	@rName: Read's name
	//	@trName: dynamic part of Read's name contained chroms number and Read's number|position,
	//	or NULL for empty name
	//	@commLine: command line to add as a comment in the first line
	inline DataOutFile(FT::eType ftype, const string& fName, const ReadName& rName, char delim=TAB)
		: _rName(rName), TxtOutFile(ftype, fName, delim) {}

	// Copy constructor for multithreading
	//	@rName: Read's name
	inline DataOutFile(const DataOutFile& file, const ReadName& rName) 
		: _rName(rName), TxtOutFile(file) {}

	// Returns a pointer to the line write buffer at current position.
	inline char* LineCurrPosBuf() const	{ return _lineBuff + _lineBuffOffset; }

	// Returns a pointer to the line write buffer at given position.
	//	@shift: given relative pointer position 
	//inline char* LineCurrPosBuf(rowlen shift) const	{ return _lineBuff + shift;	}

	// Moves current position of the line write buffer by shift.
	//inline void SlipOffset(int shift)	{ _buffLineOffset += shift; }

	// Decreases current position of the line write buffer by one.
	inline void LineDecreaseOffset() { _lineBuffOffset--; }

	// Increases current position of the line write buffer by one.
	//inline void LineIncreaseOffset() { _buffLineOffset++; }

	// Sets byte in the line write buffer to the specified value.
	//	@offset: shift of buffer start position
	//	@ch: byte to be set
	inline void LineSetChar(rowlen offset, char ch)	{ _lineBuff[offset] = ch; }

	// Adds byte to the current position in the line write buffer, adds delimiter after byte
	//	and increases current position.
	//	@ch: value to be set
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddChar(char ch, bool addDelim = false);

	// Adds byte before the current position of the line write buffer, adds delimiter before byte
	//	and decreases current position.
	//	@chr: value to be set
	//	@addDelim: if true then adds delimiter before and decreases current position
	//inline void LineAddCharBack(char ch)	{ _lineBuff[--_lineBuffOffset] = ch; }

	// Copies the string before the current position in the line write buffer,
	//	adds delimiter after string and decreases current position.
	//	@str: string to be copied
	//	@addDelim: if true then adds delimiter after string and decreases current position
	inline void LineAddStrBack(const string& str, bool addDelim = true)
	{ LineAddCharsBack(str.c_str(), str.length(), addDelim);	}

	// Copies qualified Read name and positions after current position in the line write buffer,
	//	and increases current position.
	//	@addDelim: if true then adds delimiter after string and increases current position
	inline void LineAddReadName(bool addDelim = true)
	{ LineAddStr(_rName.Name(), _rName.Length(), addDelim); }

	// Copies qualified Read name and positions after current position in the line write buffer,
	//	and increases current position.
	//	@mate: mate number for PE Reads, or 0 for SE Read
	//	Invoked in BedOutFile.
	void LineAddReadName(BYTE mate);

	// Copies qualified Read name started with '@' and read variable legth after current position
	//	in the line write buffer, and increases current position
	//	@sRLen: Read length as a string
	//	Invoked in FqOutFile.
	void LineAddReadVarName(const string& sRLen);

	// Copies qualified Read name before current position in the line write buffer,
	//	adds delimiter after Read Nameand decreases current position.
	inline void LineAddReadNameBack()
	{ LineAddCharsBack(_rName.Name(), _rName.Length(), true);	}

	// Copies qualified Read name started with '@', positions and Read constant length before current position
	//	in the line write buffer, adds delimiter after Read Name and decreases current position.
	//	Invoked in FqOutFile.
	inline void LineAddReadConstNameBack();
	
	// Adds last part of the line write buffer (from current position to the end)
	//	to the file write buffer.
	inline void LineBackToBuffer()
	{ RecordToIOBuff(_lineBuff + _lineBuffOffset, _lineBuffLen - _lineBuffOffset); }

	// Fills line by Read const quality pattern from the specified position
	inline void LineFillReadConstPatt(rowlen pos) const { RqPattern->Fill(_lineBuff + pos);	}

	// Fills line by Read variable quality pattern from the current position and increases current position
	//	@rlen: Read's length
	void LineFillReadVarPatt(readlen rlen);

#ifdef _DEBUG
	const void PrintBuff(const char* pos, size_t len) const { 
		for (int i = 0; i < len; i++)	cout << *(pos + i);	cout << LF;
	}
#endif

public:
	// Sets Read quality pattern and const Read length
	//	@rqPattFName: valid file name of the Read quality pattern
	static void Init(const char* rqPattFName) {
		RqPattern.reset(new ReadQualPattern(rqPattFName));
		if (!DistrParams::IsRVL())	sReadConstLen = ReadLenTitle + to_string(Read::FixedLen);
	}

	// Prints Read quality pattern
	inline static void PrintReadQualPatt() { RqPattern->Print(); }
};

// 'BedROutFile' implements methods for writing BaseItems file
class BedROutFile : public DataOutFile
{
private:
	rowlen _offset;		// the length of the constant chrom name in the line buffer

public:
	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name without extention
	//	@rName: Read's name
	BedROutFile(const string& fName, const ReadName& rName);

	// Clone constructor for multithreading
	inline BedROutFile(const BedROutFile& file, const ReadName& rName) : DataOutFile(file, rName) {}

	// Sets treated chrom's name to line write buffer
	void SetChrom(chrid cID);

	// Adds Read to the line's write buffer.
	//	@read: valid Read
	//	@reverse: if true then set reverse strand, otherwise set forward
	//	@mate: mate number for PE Reads, or 0 for SE Read
	void AddRead(const Read& read, bool reverse, BYTE mate=0);
};

// 'FqOutFile' implements methods for writing FQ file
class FqOutFile : public DataOutFile
// 'public' to allow implicit conversion in DataOutFile(const FqOutFile&) invoke
{
private:
	static rowlen ReadStartPos;		// fixed Read field start position

public:
	// Creates new instance for writing
	//	@fName: file name without extention
	//	@rName: Read's name
	FqOutFile(const string& fName, const ReadName& rName);

	// Forms Read from fragment and adds it to the file.
	//	@read: valid Read
	//	@reverse: if true then add complemented read 
	void AddRead(const Read& read, bool reverse);
};

// 'SamOutFile' implements methods for writing SAM file
class SamOutFile : public DataOutFile
{
private:
	static rowlen ReadStartPos;			// fixed Read field start position
	static string Fld_5_6;				// combined value from 5 to 6 field: initialised in constructor
	static string FLAG[];				// FLAG value for SE: 00000=0 (+), 10000=16 (-)
										// FLAG value for PE: 01100011=99 (+), 10010011=147 (-)
	const string Fld_7_9 = "*\t0\t0";	// combined value from 7 to 9 field for SE mode: predefined
	const BYTE tagCLlen = 3;			// length of the header line tag 'CL:'
	const char CIGAR_M = 'M';			// CIGAR marker
	
	string	_cName;						// current chrom's name

	// Adds Read with prepared 7-9 fields in local buffer to the line's write buffer.
	//	@read: valid Read
	//	@fld_7_9: prepared 7-9 fields (RNEXT,PNEXT,TLEN)
	//	@flag: FLAG field value
	void AddRead(const Read& read, const string& fld_7_9, const string& flag);

	// Returns string containing initialized PE feilds 7-9
	string GetPeFld_7_9(chrlen pos, int fLen)
	{ return "=\t" + to_string(pos) + TAB + to_string(fLen); }

public:
	// Creates new instance for writing, initializes line write buffer writes header.
	//	@fName: file name without extention
	//	@rName: Read's name
	//	@cSizes: chrom sizes
	SamOutFile(const string& fName, const ReadName& rName, const ChromSizes& cSizes);

	// Clone constructor for multithreading
	inline SamOutFile(const SamOutFile& file, const ReadName& rName) : DataOutFile(file, rName) {}

	// Sets current treated chrom
	inline void SetChrom(chrid cID) { _cName = Chrom::AbbrName(cID); }

	// Adds Read to the line's write buffer.
	//	@read: valid Read
	//	@pos: Read's start position
	//	@len: Read's length
	//	@reverse: if true then add complemented read
	inline void AddRead(const Read& read, bool reverse) { AddRead(read, Fld_7_9, FLAG[reverse]); }

	// Adds two mate Reads to the line's write buffer.
	//	@read1: valid first mate Read
	//	@read2: valid second mate Read
	//	@fLen: fragment's length
	void AddTwoReads(const Read& read1, const Read& read2, int fLen);
};

typedef UINT coval;		// coverage value

// 'Density' represents chrom's fragment coverage data
//	and implements a single method for gradual filling (incrementing) coverage
class Density : public map<chrlen, coval>
{
	typedef map<chrlen, coval>::iterator iter;

	// Adds item start position to statistics
	inline void AddPos(chrlen pos) { (*this)[pos /*+ pos%2*/]++; }

public:
	typedef map<chrlen, coval>::const_iterator citer;

	// Default constructor
	inline Density() {}

	// Adds Read to accumulate the density
	//	@tag: added Read
	//	@reverse: if true then add complemented read
	void AddRead(const Read& tag, bool reverse) { AddPos(reverse ? tag.End() : tag.Start()); }

	// Adds fragment to accumulate the density
	//	@frag: added frag
	inline void AddRead(const Region& frag) { AddPos(frag.Centre()); }

	// Adds fragment to accumulate the coverage
	void AddFragCov(const Region& frag);

#ifdef _DEBUG
	void Print() const
	{
		cout << "pos\tval\n";
		for(citer it=begin(); it!=end(); ++it)
			cout << it->first << TAB << it->second << LF;
	}

	// Prints output in BedGraph format
	void Output() const
	{
		citer it=begin();
		chrlen start = it->first;
		coval val = it->second;

		cout << "start\tend\tval\n";
		for(++it; it!=end(); ++it) {
			if(val)
				cout << start << TAB << it->first << TAB << val << LF;
			start = it->first;
			val = it->second;
		}
	}
#endif	// _DEBUG
};

// 'Coverages' represets 3 chrom's fragment coverage data pointers,
//	keeping indepentently in BedGrOutFiles duplicates
class Coverages
{
	static const BYTE Count = 3;

	Density* _covers[Count];	// binding coverages: [0] - total, [1] - pos strend, [2] - neg strand

public:
	inline Coverages() { memset(_covers, 0, Count * sizeof(Density*));	}
	
	// Bind a coverage instance by index
	//	@cover: binding coverage
	//	@i: coverage index: [0] - total, [1] - pos strend, [2] - neg strand
	inline void SetCover(Density* cover, BYTE i) { _covers[i] = cover; }

	// Adds SE fragment to coverage
	//	@frag: added fragment
	//	@reverse: true if read is reversed (neg strand)
	void AddFrag(const Region& frag, bool reverse);

	// Adds PE fragment to coverage
	//	@frag: added fragment
	inline void AddFrag(const Region& frag) { _covers[Count-1]->AddFragCov(frag); }
};

// Density representative (factory): 
//	implements methods for Density creating, deleting and icnrement control
struct DensFactory
{
	Density* Dens;	// Density; initialized only while given chrom treatment
	bool	Closed;		// true if coverage is closed for increment
	bool	Saved;		// true if coverage is recorded

	inline DensFactory() : Dens(NULL) { Closed = Saved = false; }

	inline Density* CreateDens() { return Dens = new Density(); }

	inline void RemoveDens() { if (Dens) delete Dens, Dens = NULL; }

	// Returns true if coverage is closed but unsaved
	inline bool IsUnsaved() { return !Saved && Closed; }
};

// 'WigOutFile' is a base class for BedGrOutFile (BedGraph format) and DensOutFile (variableStep format)
class WigOutFile : public Chroms<DensFactory>, public TxtOutFile
{
	static const string WigFormats[];

	Mutex::eType _mType;

	// Fill IO buffer by chrom's data
	virtual void PutChromData(chrid cID, DensFactory& df) {}

	// Writes data for given chrom to file
	void WriteChromData(chrid cID, DensFactory& df);

	//virtual void AssignData(WigOutFile& file) {}

	// Initializes line write buffer, adds command and definition lines
	//	@fName: file name without extention
	//	@declDescr: brief file description in declaration line
	//	@strand: string denoted strand, or empty string
	void Init(FT::eType ftype, const string& fName, const string& declDescr, const string& strand = strEmpty);

public:
	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name without extention
	//	@descr: brief file description in declaration line
	//	@mtype: mutex locker type
	//	@cSizes: chrom sizes
	WigOutFile(FT::eType ftype, const string& fName, const string& descr, Mutex::eType mtype, const ChromSizesExt& cSizes);

	// Creates new strand-separated BedGraph instance for writing and initializes line write buffer.
	//	@strandInd: strand index: 0 - positive, 1 - negative
	//	@fName: file name without extention
	//	@strand: string denoted strand, or empty string
	WigOutFile(int strandInd, const string& fName, const string& strand) :
		_mType(Mutex::eType(int(Mutex::eType::WR_BG) + strandInd)),
		TxtOutFile(FT::eType::BGRAPH, fName, TAB)
	{ Init(FT::eType::BGRAPH, fName, strand); }

	~WigOutFile();

	// Creates and returns data container for given chrom
	inline Density* OpenChromData(chrid cID) { return At(cID).Data.CreateDens(); }

	// Close data container for given chrom for increment.
	//	@cID: chrom ID
	void CloseChromData(chrid cID);
};

// 'BedGrOutFiles' implements methods for writing 1 or 3 WIG files
class BedGrOutFiles
{
	// 'BedGrOutFile' implements methods for writing BedGraph file
	class BedGrOutFile : public WigOutFile
	{
		static const string DeclDescr;		//	brief file description in declaration line
		static const string StrandTitles[];

		// Fill IO buffer by chrom data
		void PutChromData(chrid cID, DensFactory& df);

	public:
		// Creates new instance for writing and initializes line write buffer.
		//	@fName: file name without extention
		//	@cSizes: chrom sizes
		BedGrOutFile(const string& fName, const ChromSizesExt& cSizes)
			: WigOutFile(FT::eType::BGRAPH, fName, DeclDescr, Mutex::eType::WR_BG, cSizes) {}

		// Creates new strand-separated BedGraph instance for writing and initializes line write buffer
		//	@strandInd: strand index: 0 - positive, 1 - negative
		//	@fName: file name without extention
		//	@strand: string denoted strand, or empty string
		//	@file: united strands file
		BedGrOutFile(int strandInd, const string& fName, const BedGrOutFile& file)
			: WigOutFile(strandInd, fName, StrandTitles[strandInd])
		{ Assign(file); }
	};

	static const BYTE Count = 3;

	BedGrOutFile *_files[Count];

public:
	static bool IsStrands;			// true if wigs with different strands should be generated

	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name with extention
	//	@cSizes: chrom sizes
	BedGrOutFiles(const string& fName, const ChromSizesExt& cSizes);

	~BedGrOutFiles();

	// Starts accumalating coverage for given chrom
	//	@cID: chrom
	//	@covrs: extern Density pointers
	void StartChrom(chrid cID, Coverages& covrs);

	// Stops accumalating coverage for given chrom
	void StopChrom(chrid cID);

	// Prints output file names separated by comma
	//	@signOut: output marker
	void PrintNames() const;
};

// 'DensOutFile' implements methods for writing density in WIG (wiggle_0) file
class DensOutFile : public WigOutFile
{
	// Fill IO buffer by chrom data
	void PutChromData(chrid cID, DensFactory& df);

public:
	// Creates new instance for writing and initializes line write buffer.
	//	@fName: file name without extention
	//	@cSizes: chrom sizes
	DensOutFile(const string& fName, bool rdens, const ChromSizesExt& cSizes)
	: WigOutFile(
		FT::eType::WIG_VAR,
		fName + (rdens ? ".rdens" : ".fdens"),
		rdens ? "read density" : "frag density",
		Mutex::eType::WR_RDENS,
		cSizes
	) {}
};

// 'Output' wraps test and control output files
class Output
{
public:
	// Output file formats
	enum class eFormat {	// values are masks in eFormat variable
		UNDEF	= 0,
		FG		= 0x01,
		BED		= 0x02,
		SAM		= 0x04,
		BGR		= 0x08,		// bedGraph
		FDENS	= 0x10,
		RDENS	= 0x20,
		FDIST	= 0x40,
		RDIST	= 0x80,
	};

	static string	MapQual;	// the mapping quality
	//static bool	RandomReverse;	// true if Read should be reversed randomly
private:
	static eFormat Format;				// output formats
	static const char* entityTitles[];	// entity titles for printing
	static const BYTE ND = 2;			// number of distribution/density files

	// 'OutFile' wraps output files: FQ, SAM, BED, WIG.
	class OutFile
	{
		typedef int	 (OutFile::* tpAddRead)(const Region&, readlen, bool);

		// pointer to the 'add Read' method
		static tpAddRead pAddRead;
		static float	 StrandErrProb;		// the probability of strand error

		const RefSeq* _seq;
		FqOutFile	* _fqFile1;		// mate1 or single output FQ
		FqOutFile	* _fqFile2;		// mate2 output FQ
		BedROutFile	* _bedFile;		// output BED
		SamOutFile	* _samFile;		// output SAM
		BedGrOutFiles* _bgFile;		// output BedGraph
		DensOutFile	* _densFile[ND];// output density WIG
		Coverages	_covers;		// coverage data pointers for BedGraph; each own for each clone
		Density*	_dens[ND];		// output density
		ReadName	_rName;			// Read's name; local for clone independence
		bool		_primer;		// true if file is primer (not clone); only for BedGrOutFile

		// Adds one SE Read
		//	@frag: added fragment
		//	@rLen: Read's length
		///	@g: FG or BG; needs for strand error imitation
		//	@reverse: if true then add complemented read
		//	return:	1: fragment is out of range (end of chrom)
		//			0: Read is added successfully
		//			-1: N limit is exceeded
		int AddReadSE (const Region& frag, readlen rlen, /*Gr::eType g,*/ bool reverse);

		// Adds two PE Reads
		//	@frag: added fragment
		//	@rLen: Read's length
		///	@g: FG or BG; needs for strand error imitation; not used
		//	@reverse: not used
		//	return:	1: fragment is out of range (end of chrom)
		//			0: Reads are added successfully
		//			-1: N limit is exceeded
		int AddReadPE (const Region& frag, readlen rlen, /*Gr::eType g,*/ bool reverse);

		// Empty (trial) method
		//inline int NotAddRead (Output*, const Region&, /*Gr::eType,*/ bool)	{ return 0; }
		inline int NotAddRead(const Region&, readlen, /*Gr::eType,*/ bool) { return 0; }

	public:
		// Initializes static members
		//	@singleThread: true if single thread is set
		//	@sErrProb: the probability of strand error
		inline static void Init(float sErrProb) { StrandErrProb = sErrProb; }

		// Sets sequense mode.
		//	@trial: if true, then set empty mode, otherwise current working mode
		static void SetSeqMode(bool trial)	{ 
			pAddRead = trial ? &Output::OutFile::NotAddRead :
				(Seq::IsPE() ? &Output::OutFile::AddReadPE : &Output::OutFile::AddReadSE);
		}

		// Creates and initializes new instance for writing.
		//	@fName: common file name without extention
		//	@cSizes: chrom sizes
		OutFile(const string& fName, const ChromSizesExt& cSizes);

		// Clone constructor for multithreading
		//	@oFile: original instance
		//	@threadNumb: number of thread
		OutFile(const OutFile& oFile);
	
		~OutFile();

		// Start recording chrom
		void BeginWriteChrom(const RefSeq& seq);

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
		inline int AddRead(const Region& frag, readlen rlen, /*Gr::eType g,*/ bool reverse) {
			return (this->*pAddRead)(frag, rlen, reverse);
		}

		// Prints output file formats and sequencing mode
		//	@signOut: output marker
		//	@predicate: 'output' marker
		void PrintFormat(const char* signOut, const char* predicate) const;
	};

	// 'DistrFiles' manages two distribution: fragments and reads
	class DistrFiles
	{
		static const string fExt[];			// distrib file extentios
		static const char* entityAdjust[];	// entity titló adjustment for printing

		const string _fName = strEmpty;		// common part of frag's/read's distribution file name
		LenFreq* _dist[ND] = {NULL, NULL};	// distributions: fragments (0) and reads (1)

		// Returns file name for fragments (ind=0) and reads (ind=1) distributions
		const string FileName(BYTE ind) const {
			return _fName + fExt[ind]
#ifdef OS_Windows
				+ ".txt"
#endif
				;
		}

	public:
		DistrFiles(const string& fName, bool isFragDist, bool isReadDist);

		// Writes distributions to files and delete them
		~DistrFiles();

		// Adds frag/read length to statistics
		void AddFrag(fraglen flen, readlen rlen);

		// Prints output file formats
		//	@signOut: output marker
		//	@predicate: output common title
		void PrintFormat(const char* signOut, const char* predicate) const;
	};

	// Returns true if formats f1 or f2 or f3 are defined
	static bool HasFormat(eFormat f1, eFormat f2 = eFormat::UNDEF, eFormat f3 = eFormat::UNDEF)
	{ return int(Format) & (int(f1) | int(f2) | int(f3)); }

	// Returns true if at least one of detailed data formats (BG, BED, SAM) is defined
	static bool HasDataFormat()
	{ return int(Format) & (int(eFormat::FG) | int(eFormat::BED) | int(eFormat::SAM)); }

	// Returns true if contiguous format shifted by i from the base f is defined
	static bool HasContigFormat(eFormat f, BYTE i)	{ return HasFormat(eFormat(BYTE(f) * 1 << i)); }
	
	// Returns true if both formats is defined
	static bool HasBothFormats(eFormat f1, eFormat f2) { return OnesCount(BYTE(Format) & (BYTE(f1) | BYTE(f2))) == 2;	}
	
	// Returns true if only WIG format is defined
	//static inline bool HasWigOnly() { return Format == eFormat::BGR; }

	unique_ptr<OutFile> _oFiles[2];	// test [0] and control [1] OutFile objects
	shared_ptr<DistrFiles> _dists;	// frag's & read's distribution
	Random	_rng;					// random generator; needed for Read variable length generation
	GM::eMode	_gMode;				// current generation mode,
									// corresponding to the index in the files and counters arrays

	// Return number of bits in an integer
	static int OnesCount(int n);

	//static void PrintPairFormat(eFormat f1, eFormat f2, 
	//	const char* signOut, const char* predicate, const char* type, const string titles[], 
	//	const string const Name(BYTE i) )
	//{
	//	if (HasFormat(f1, f2)) {
	//		std::cout << signOut << predicate << type << SepDCl;
	//		for (BYTE i = 0; i < ND; i++) {
	//			if (HasContFormat(FDENS, i))
	//				std::cout << signOut << predicate[i] << SepCl << Name(i);
	//			if (i == 0 && HasBothFormats(FDENS, RDENS))
	//				std::cout << SepSCl;
	//		}
	//		std::cout << LF;
	//	}
	//}

	// Prints item title or count
	template <typename T>
	static void PrintItemsSummary(T t1, T t2) {
		bool isReadFormat = HasDataFormat() || HasFormat(eFormat::RDENS, eFormat::RDIST);
		if (isReadFormat)
			cout << t1;
		if (HasFormat(eFormat::BGR, eFormat::FDENS, eFormat::FDIST)) {
			if (isReadFormat)	cout << '/';
			cout << t2;
		}
	}

public:
	// Returns true if SAM type is assigned.
	//static inline bool IsSamSet() { return HasFormat(eFormat::SAM); }

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
	inline static void SetSeqMode(bool trial)	{ OutFile::SetSeqMode(trial); }

	// Sets Read quality pattern by valid file name.
	inline static void SetReadQualPatt(const char* rqPattFName) { 
		DataOutFile::Init(rqPattFName);
	}

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
	Output(const string& fName, bool control, const string& cmLine, const ChromSizesExt& cSizes);

	// Clone constructor for multithreading.
	//	@file: original instance
	Output(const Output& file);
	
	 // Set generation mode
	//	@testMode: if true, set Test mode, otherwhise Control mode
	inline void SetGMode(GM::eMode gm) { _gMode = gm; }

	// Starts recording chrom
	void BeginWriteChrom(const RefSeq& seq);

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
	void PrintFormat	(const char* signOut) const;

	// Prints Read quality settins
	//	@signOut: output marker
	void PrintReadQual	(const char* signOut) const;
};
