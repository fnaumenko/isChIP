#pragma once
#include "Data.h"

typedef short fraglen;	// type length of fragment

class FqFile : public TxtFile
/*
 * 'FqFile' implements writing files in FQ format.
 */
{
#ifdef _FQSTATN

private:
	enum eLineLenIndex { HEADER1, READ, HEADER2, QUAL };

public:	
	// Creates new instance by file name and type : general constructor
	inline FqFile(const string & fileName, eType type)
		: TxtFile(fileName, type == NEW ? TxtFile::WRITE : TxtFile::READ, 4)
	{ Init(); }

	// Returns checked length of current readed Read.
	readlen ReadLength() const;
	
	// Gets checked Read from readed Sequence.
	char* GetCurrRead() const;

	// Returns checked Sequence.
	char*	GetSequence	();
	
#endif	// _FQSTATN

#ifdef _FILE_WRITE

private:
	static string Ext;				// extention of output file (.fq or .fq.gz)
	static rowlen ReadStartPos;		// constant Read field start position

	// Presets line write buffer
	void InitBuffer();

public:
	// Creates new instance for writing
	//	@fName: file name without extention
	//	@zipExt: zip extention if file should be zipped; otherwise empty string
	inline FqFile(const string& fName, const string& zipExt)
		: TxtFile(fName + Ext + zipExt, TxtFile::WRITE, 4) {}

	// Creates new instance - mate - with predefined file's name and extention
	//	@fName: file name without extention
	//	@zipExt: zip extention if file should be zipped; otherwise empty string
	//	@mateNumb: mate number: 1 or 2
	inline FqFile(const string& fName, const string& zipExt, BYTE mateNumb)
		: TxtFile(fName + USCORE + BSTR(mateNumb) + Ext + zipExt, TxtFile::WRITE, 4) {}

#ifdef _MULTITHREAD
	// Creates a clone of existed @file for writing only.
	//	@threadNumb: used for forming buffer with differ sizes for anisochronous writing
	inline FqFile(const FqFile& file, threadnumb threadNumb) : TxtFile(file, threadNumb)
	{ InitBuffer(); }
#endif

	// Initializes line write buffer
	void InitToWrite ();

	// Forms Read from fragment and adds it to the file.
	//	@rName: Read's name
	//	@read: valid read
	//	@reverse: if true then complement added read 
	//	@mate: SINGLE for one-side sequensing, or MATE_FIRST | MATE_SECOND for paired-end
	void AddRead	(const string& rName, const char* read, bool reverse, eMate mate);

#endif	// _FILE_WRITE
};

class BedRFile : public TxtFile
/*
 * 'BedRFile' implements Bed file for writing.
 */
{
private:
	static string Ext;			// Extention of output BED file .bed or .bed.gz)
	static const string Score;	// Read's max score

	rowlen _offset;				// current writing position in line write buffer

public:
	// Creates new instance for writing.
	//	@fName: file name without extention
	//	@zipExt: zip extention if file should be zipped; otherwise empty string
	inline BedRFile(const string& fName, const string& zipExt)
		: TxtFile(fName + Ext + zipExt, WRITE, 1) {}
	 
#ifdef _MULTITHREAD
	// Creates a clone of existed instance for writing.
	//	@file: original instance
	//	@threadNumb: number of thread
	inline BedRFile(const BedRFile& file, threadnumb threadNumb)
		: TxtFile(file, threadNumb) {}
#endif

	// Initializes line write buffer; only for master, clones are initialized by master
	//	@commandLine: command line to add as comment on the first line
	void InitToWrite (const string& commandLine);

	// Sets chrom's name to line write buffer
	void BeginWriteChrom(chrid cID);

	// Adds Read to the line's write buffer.
	//	@rName: Read's name
	//	@pos: valid Read's start position
	//	@reverse: if true then set reverse strand, otherwise set forward
	void AddRead(const string& rName, chrlen pos, bool reverse);

	// Adds two mate Reads to the line's write buffer.
	//	@rName: Read's name
	//	@pos1: valid first mate Read's start position
	//	@pos2: valid second mate Read's start position
	void AddTwoReads(const string& rName, chrlen pos1, chrlen pos2);
};

class SamFile : public TxtFile
/*
 * 'SamFile' implements SAM file for writing.
 */
{
private:
	static string Ext;				// Extention of output SAM file (.sam or .sam.gz)
	static rowlen ReadStartPos;		// constant Read field start position
	static string Comb5_6;			// combined value from 5 to 6 field: defined in constructor
	const static string Comb7_9;	// combined value from 7 to 9 field for SE mode: predefined
	static rowlen Offset5_9;		// length of predefined values from 5 to 8; 0 for PE mode
	static string Flag[];			// FLAG value for SE: 01100011->99 (+), 10010011->147 (-)

	string	_cName;			// current chrom's name
	BYTE	_headLineCnt;	// number of written lines in header

	// Writes header line to line write buffer.
	//	@tag0: line tag
	//	@tag1: first subtag
	//	@val1: first value
	//	@tag2: second subtag
	//	@val2: second value
	//	@closeLine: true if line is ended
	void SetHeaderLine(
		const char* tag0,
		const char* tag1, const string& val1,
		const char* tag2=NULL, const string& val2=StrEmpty,
		bool closeLine=true
		);

	// Generates and writes SAM header
	//	@cSizes: chrom sizes
	//	@commandLine: command line
	void CreateHeader(const ChromSizes& cSizes, const string& commandLine);

	// Adds full-defined Read to the line's write buffer.
	//	@rName: Read's name
	//	@read: valid Read
	//	@flag: FLAG field value
	//	@pos1: valid start position of mate Read (PE) or Read (SE)
	//	@pos2: valid start position of mate Read (PE) or -1 (SE)
	//	@fLen: +|- fragment's length (PE) or 0 (SE)
	void AddStrongRead(const string& rName, const char* read, const string& flag,
		chrlen rPos1, chrlen rPos2 = 0, int fLen = 0);

	// Creates and initialize line write buffer.
	void InitBuffer();

public:
	// Creates new instance for writing.
	//	@fName: file name without extention
	//	@zipExt: zip extention if file should be zipped; otherwise empty string
	SamFile(const string& fName, const string& zipExt)
		: _headLineCnt(0), TxtFile(fName + Ext + zipExt, WRITE, 1) {}

#ifdef _MULTITHREAD
	// Creates a clone of existed instance for writing.
	//	@file: original instance
	//	@threadNumb: number of thread
	SamFile(const SamFile& file, threadnumb threadNumb)	: TxtFile(file, threadNumb)
	{ InitBuffer();	}
#endif

	// Initializes line write buffer and header and makes ready for writing
	//	@cSizes: chrom sizes
	//	@commandLine: command line
	void InitToWrite(const ChromSizes& cSizes, const string& commandLine);

	// Sets current chrom
	void BeginWriteChrom(chrid cID) { _cName = Chrom::AbbrName(cID); }

	// Adds Read to the line's write buffer.
	//	@rName: Read's name
	//	@read: valid Read
	//	@pos: valid Read's start position
	//	@reverse: if true then set reverse strand, otherwise set forward
	void AddRead(const string& rName, const char* read, chrlen pos, bool reverse);

	// Adds two mate Reads to the line's write buffer.
	//	@rName: name of Read
	//	@read1: valid first mate Read
	//	@read2: valid second mate Read
	//	@pos1: valid first mate Read's start position
	//	@pos2: valid second mate Read's start position
	//	@fLen: fragment's length
	void AddTwoReads(const string& rName, const char* read1, const char* read2,
		chrlen pos1, chrlen pos2, int fLen);

	// Gets count of data lines.
	inline ULONG Count() const { return RecordCount() - _headLineCnt; }
};


class OutFile
/*
 * 'OutFile' wraps output files: FQ|SAM|BED.
 */
{
public:
	// Output file formats
	enum eFormat {
		ofFQ	= 0x1,
		ofBED	= 0x2,
		ofSAM	= 0x4
	};

	// Output file formats
	enum eMode {
		mSE	= 0,
		mPE	= 1,
		mEmpty = 2
	};

private:
	static eMode	Mode;	// 0 in case of one-side sequencing, 1 in case of paired-end
	
	typedef int	(OutFile::*AddReads)(string&, const Nts&, ULONG, chrlen, fraglen, bool);
	static AddReads callAddRead[];	// 0: 'add SE Read' method,
									// 1: 'add PE Read' method,
									// 3: empty method (for trial cutting)
	eMode		_mode;
	FqFile	*	_fqFile1;	// mate1 or single FQ
	FqFile	*	_fqFile2;	// mate2 FQ
	BedRFile*	_bedFile;
	SamFile	*	_samFile;

	// Adds one SE Read
	int AddReadSE (string& rName, const Nts& nts,
		ULONG rNumb, chrlen pos, fraglen fragLen, bool reverse);
	
	// Adds two PE Reads
	int AddReadPE (string& rName, const Nts& nts,
		ULONG rNumb, chrlen pos, fraglen fragLen, bool reverse);

	// Empty (trial) method
	inline int NoAddRead (string& rName, const Nts& nts,
		ULONG rNumb, chrlen pos, fraglen fragLen, bool reverse) { return 2; }

public:
	// Returns 0 in case of one-side sequencing, 1 in case of paired-end
	static inline BYTE	PairedEnd()	{ return Mode; }

	// Creates new instance for writing.
	//	@fName: common file name without extention
	//	@outType: types of output files
	//	@mode: SE or PE
	//	@isZipped: true if output files should be zipped
	//	@commandLine: command line
	OutFile(const string& fName, eFormat outType, eMode mode, bool isZipped);

#ifdef _MULTITHREAD
	// Creates a clone of existed instance for writing.
	//	@oFile: original instance
	//	@threadNumb: number of thread
	OutFile(const OutFile& oFile, threadnumb threadNumb);
#endif
	
	~OutFile();

	// Initializes buffers and makes ready for writing
	//	@cSizes: chrom sizes or NULL
	//	@commandLine: command line
	void Init(const ChromSizes* cSizes, const string& commandLine);

	// Sets/clears empty mode.
	// In empty mode no output is produced.
	//	@val: if true, than set empty mode, otherwise working mode
	inline void SetEmptyMode(bool val)	{ _mode = val ? mEmpty : Mode; }

	// Returns count of writed Reads.
	ULONG Count() const;

	// Returns true if SAM type is assigned.
	inline bool IsSamSet()	{ return (bool)_samFile; }

	// Adds read(s) to output file
	//	@cName: chrom's name
	//	@nts: cutted chromosome
	//	@rNumb: current number of writed Read|pair of Reads, or 0 if Read named by position
	//	@pos: current fragment's position
	//	@fragLen: length of current fragment
	//	@reverse: true if read is reversed (neg strand), otherwise read is forward (pos strand)
	//	return:	-1 if fragment is out of range,
	//			0 if limitN is exceeded,
	//			1 if Read(s) is(are) added,
	//			2 if not produce the output file
	int AddRead (const string& cName, const Nts& nts,
		ULONG rNumb, chrlen pos, fraglen fragLen, bool reverse)
	{
		string rName = Read::Name() + COLON + cName;
		return (this->*callAddRead[int(_mode)])(	// AddReadSE, AddReadPE or NoAddRead
			rName,
			//Read::Name() + COLON + cName,			// doesn't allowed by g++
			nts, rNumb, pos, fragLen, reverse);
	}

	// Sets chrom's name for writing.
	void BeginWriteChrom(chrid cID) const {
		if(_bedFile)	_bedFile->BeginWriteChrom(cID);
		if(_samFile)	_samFile->BeginWriteChrom(cID);
	}

	// Finishes writing to file and close it.
	void Write() const;

	// Prints output files as head info
	//	@signOut: output marker
	void Print(const char* signOut) const;
};

class SamFragDistr
{
public:
	SamFragDistr(const string& fname);
};
