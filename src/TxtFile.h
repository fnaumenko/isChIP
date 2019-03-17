#pragma once
#include "common.h"

#define INT_CAPACITY	10		// maximal number of digits in integer

// Number of basics file's reading|writing buffer blocks.
// Should be less than 2047 because of ULONG type of block size variable.
// Otherwise the behaviour is unpredictable.

typedef short rowlen;	// type: length of row in TxtFile

// 'TabFilePar' keeps basic parameters for TabFile
struct TabFilePar {
	const BYTE	MinFieldCnt;	// minimum number of feilds in file line; these fields are checked during initialization
	const BYTE	MaxFieldCnt;	// maximum possible number of feilds in data line; these fields are checked by a call
	const char Comment;			// char indicates that line is comment
	const char* LineSpec;		// substring on which each data line is beginning

	inline TabFilePar(BYTE	minTabCnt, BYTE maxTabCnt, char comm, const char* lSpec) :
		MinFieldCnt(minTabCnt),
		MaxFieldCnt(maxTabCnt<minTabCnt ? minTabCnt : maxTabCnt),
		Comment(comm),
		LineSpec(lSpec)
		{}
};

// 'File Type' implements bioinformatics file type routines 
static class FT
{
private:
	struct fType {
		const char* Extens;			// file extension
		const string Item;			// item title
		const string ItemPl;		// item title in plural
		TabFilePar FileParam;		// TabFile parameters, defined feilds
	};
	static const fType	Types[];
	static const BYTE	Count;

	// Validates file format
	//	@fName: file name (with case insensitive extension and [.gz])
	//	@t: file type
	//	return: true if file extension correspondes to file type
	//static bool CheckType(const char* fName, eTypes t) { 
	//	return GetType(fName) == (t == ABED ? BED : t); 
	//}

public:
	enum eTypes {
		UNDEF,	// undefined type
		ABED,	// alignment bed
		BED,	// ordinary bed
		WIG,	// wiggle
		SAM,	// sam
		FQ,		// fastQ
		FA		// fasta
	};

	// Gets file format
	//	@fName: file name (with case insensitive extension)
	static eTypes GetType(const char* fName);

	// Validates file extension
	//	@fName: file name (with case insensitive extension and [.gz])
	//	@t: file type
	//	@printfName: true if file name should be ptinted
	//	@throwExc: true if throw exception, otherwise throw warning
	//	return: true if file extension correspondes to file type
	static bool CheckType(const char* fName, eTypes t, bool printfName=true, bool throwExc=true);

	// Gets file extension
	//	@t: file type
	//	@isZip: true if add ".gz"
	static const string Ext(eTypes t, bool isZip = false)	{ return Types[t].Extens; }
#ifdef _ISCHIP
	// Gets file extension, beginning at DOT and adding .gz if needed
	//	@t: file type
	//	@isZip: true if add ".gz"
	static const string RealExt(eTypes t, bool isZip = false);
#endif
	// Gets an item's title
	//	@t: file type
	//	@pl: true if plural form
	static const string& ItemTitle(eTypes t, bool pl = false)
	{ return pl ? Types[t].ItemPl : Types[t].Item; }

	// Gets TabFile params
	//	@t: file type
	static inline TabFilePar FileParams(eTypes t) { return Types[t].FileParam; }

} fformat;

class TxtFile
/*
 * Basic class 'TxtFile' implements a fast buffered serial (stream) reading/writing text files 
 * (in Windows and Linux standart).
 * Supports reading/writing zipped (gz) files.
 * Optimised for huge files.
 * Restriction: the default size of buffer, setting as NUMB_BLK * BasicBlockSize, 
 * should be bigger than the longest line in file. Otherwise file become invalid,
 * and less than UINT.
 * If size of reading files is less than default buffer's size, 
 * the buffer's size sets exactly to be sufficient to read a whole file.
 * The reading/writing unit is a 'record', which in common case is a predifined set of lines.
 * For FQ files 'record' is a set of 4 lines.
 * For common text files 'record' is identical to 'line'.
 * Empty lines are skipping by reading and therefore will not be writing by 'cloning' a file.
 */
{
public:
	enum eAction { 
		READ,		// reads only existing file
		WRITE,		// creates file if it not exist and writes to it; file is cleared before
		READ_ANY	// creates file if it not exist and reads it
	};

//private:
//protected:
public:
	enum eFlag {			// signs of file
		// The first two right bits are reserved for storing the length of the EOL marker: 1 or 2 
		// The first bit is set to 1 in the constructor.
		// If CR symbol is found at the end of line,
		//	the second bit is raised to 1, the first turn down to 0
		EOLSZ		= 0x003,	// mask for the 2 first bits
		EOLCHECKED	= 0x004,	// the presence of a symbol CR is checked; for Reading mode
		ZIPPED		= 0x008,	// file is zipped
		ABORTING	= 0x010,	// invalid file should be completed by throwing exception; for Reading mode
		ENDREAD		= 0x020,	// last call of GetRecord() has returned NULL; for Reading mode
		PRNAME		= 0x040,	// print file name in exception's message; for Reading mode
		MTHREAD		= 0x080,	// file in multithread mode: needs to be locked while writing
		CLONE		= 0x100,	// file is a clone
		//CONSTIT		= 0x200	// file is a constituent of an aggregate file
	};

private:
	const static int BasicBlockSize = 1024 * 1024;	// basic block 1 Mb
	const static int NUMB_BLK = 4;	// size of basic block in Mb
	
	LLONG	_fSize;			// the length of uncompressed file; for zipped file more than
							// 4'294'967'295 its unzipped length is unpredictable
	string	_fName;			// file's name
	mutable short _flag;	// bitwise storage for signs included in eFlag

protected:
	void *	_stream;			// FILE* (for unzipped file) or gzFile (for zipped file)
	char *	_buff;				// basic I/O (read/write) buffer
	UINT	_buffLen;			// the length of the basic I/O buffer
	mutable UINT _currRecPos;	// start position of the last readed/writed record in the block
	mutable ULONG _recCnt;		// local counter of readed/writed records
	mutable Err::eCode	_errCode;
	//Stopwatch	_stopwatch;
	
	inline void RaiseFlag(eFlag f) const	{ _flag |= f; }
	inline void SetFlag	(eFlag f, bool val)	{ val ? _flag |= f : _flag &= ~f; }
	inline bool IsFlag	(eFlag f)	const	{ return (_flag & f) != 0; }	// != 0 to avoid warning C4800
	inline bool IsZipped()			const	{ return IsFlag(ZIPPED); }

	// *** 3 methods used by TxtInFile only

	// Gets the number of characters corresponded to EOL
	// In Windows EOL matches '\r\n', in Linux EOL matches '\n',
	//	return: in Windows always 1, in Linux: 2 for file created in Windows, 1 for file created in Linux
	inline UINT EOLSize() const	{ return _flag & EOLSZ; }

	// Returns true if EOL size is not defined
	inline bool IsEOLnotDef() const { return (_flag & EOLCHECKED) == 0; }

	// Establishes the presence of CR symbol at the end of line.
	//	@c: if c is CR then the second bit is raised to 1, the first turn down to 0,
	//	so the value return by EOLSZ mask is 2, otherwise remains in state 1
	inline void SetEOL(char c)	{ if(c==CR)	_flag ^= EOLSZ;	RaiseFlag(EOLCHECKED); }

private:
	// Initializes instance variables, opens a file, sets a proper error code.
	//	@fName: valid full name of file
	//	@mode: opening mode
	//	@flStream: clonable file stream or NULL
	//	return: true is success, otherwise false.
	bool SetBasic(const string& fName, eAction mode, void* flStream);
	
	// Allocates memory for the I/O buffer with checking.
	//	return: true if successful
	bool CreateIOBuff();

	// Returns file name if name is printable, otherwise empty string.
	inline const string& FileNameToExcept() const { return IsFlag(PRNAME) ? _fName : strEmpty; }

protected:
	// Constructs an TxtFile instance: allocates I/O buffer, opens an assigned file.
	//	@fName: valid full name of assigned file
	//	@mode: opening mode
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	//	@rintName: true if file name should be printed in the exception's message
	TxtFile(const string& fName, eAction mode, bool abortInvalid=true, bool printName=true);

#ifdef _MULTITHREAD
	// Constructs a clone of an existing instance.
	// Clone is a copy of opened file with its own separate I/O buffer.
	// Used for multithreading file recording
	//	@file: opened file which is cloned
	TxtFile	(const TxtFile& file);
#endif

	// Destructs the instance: close assigned file, releases I/O buffer.
	~TxtFile();

	// Sets error code and throws exception if it is allowed.
	void SetError(Err::eCode errCode) const;

	// Returns true if instance is a clone.
	inline bool IsClone() const { return IsFlag(CLONE); }

	// Returns true if instance is valid.
	inline bool IsGood() const	{ return _errCode == Err::NONE; }

	// Returns true if instance is invalid.
	inline bool IsBad()	const	{ return _errCode != Err::NONE; }

	// Returns current error code.
	//inline Err::eCode ErrCode() const	{ return _errCode; }

	// Throws exception
	//	@msg: exception message
	//	@isExc: if true then throw exception, otherwise warning
	inline void ThrowExcept(const string& msg, bool isExc = true) const {
		Err(msg, FileNameToExcept()).Throw(isExc);
	}

public:
	// Gets file name.
	inline const string& FileName() const { return _fName; }

	// Gets size of uncompressed file,
	// or size of compressed file if its uncompressed length is more than UINT_MAX.
	// So it can be used for estimatation only.
	inline LLONG Length() const	{ return _fSize; }

	// Gets number of readed/writed records.
	ULONG RecordCount() const { return _recCnt; }
};

class TxtInFile : public TxtFile
{
private:
	UINT	_recLen;		// the length of record with EOL marker
	UINT	_readedLen;		// number of actually readed chars in block
	UINT *	_linesLen;		// array of lengths of lines in a record
	BYTE	_linesInRecCnt;	// number of lines in a record

	// Raises ENDREAD sign an return NULL
	inline char* ReadingEnded()		  { RaiseFlag(ENDREAD); return NULL; }

	// Reads next block.
	//	@offset: shift of start reading position
	//	return: 1 if file is not finished; 0 if it is finished; -1 if unsuccess reading
	int ReadBlock(const UINT offset);

	// Returns file name or empty string depends on if name is printing
	const string& FileNameToExcept() const { return IsFlag(PRNAME) ? FileName() : strEmpty; }

	// Gets string containing file name and current record number.
	//inline const string RecordNumbToStr() const
	//{ return LineNumbToStr(0); }

	// Gets number of line.
	//	@lineInd: index of line in a record
	inline ULONG LineNumber	(BYTE lineInd) const { 
		return (_recCnt-1)*_linesInRecCnt + lineInd + 1;
	}

protected:
	// Constructs an TxtInFile instance: allocates buffers, opens an assigned file.
	//	@fName: valid full name of assigned file
	//	@mode: opening mode: READ or READ_ANY
	//	@cntRecLines: number of lines in a record
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	//	@rintName: true if file name should be printed in an exception's message
	TxtInFile(const string& fName, eAction mode, BYTE cntRecLines,
		bool abortInvalid=true, bool printName=true);

	inline ~TxtInFile() { if(_linesLen)	delete [] _linesLen; }

	
	// Gets length of current reading record including EOL marker.
	inline UINT	RecordLength() const	{ return _recLen; }

	// Gets current record.
	inline const char* Record() const { 
		return IsFlag(ENDREAD) ? NULL : _buff + _currRecPos - _recLen;
	}

	// Gets next record.
	//	return: point to the next record in a buffer.
	inline const char* NextRecord() const { return _buff + _currRecPos; }

	// Sets _currLinePos to the beginning of next non-empty line inread/write buffer.
	//	@counterN: if not NULL, adds to counterN the number of 'N' in a record. Used in Fa() only.
	//	@posTab: if not NULL, sets TABs positions in line to this array
	//	@cntTabs: if 'posTab' is not NULL, the length of 'posTab' array
	//	return: point to line or NULL if no more lines
	char*	GetRecord(chrlen* const counterN=NULL, short* const posTab=NULL, BYTE cntTabs=0);
	
	// Gets string containing file name and current line number.
	//	@lineInd: index of line in a record; if 0, then first line
	const string LineNumbToStr(BYTE lineInd = 0) const;

#ifdef _FQSTATN
	// Gets length of line.
	//	@lineInd: index of line in a record
	//	@withEOL: if true then with EOL marker length, ptherwise without it
	inline UINT	LineLengthByInd	(BYTE lineInd, bool withEOL = false) const { 
		return _linesLen[lineInd] - (withEOL ? 0 : EOLSize());
	}

	// Throw exception if no record is readed.
	void CheckGettingRecord() const { 
		if( IsEOLnotDef() )
			Err("attempt to get record's info without reading record", FileName().c_str()).Throw();
	}
#endif

public:
	// Throws exception with message as file name followed by exception message
	//	@msg: exception message
	inline void ThrowExcept(const string& msg) const {
		Err(msg, FileNameToExcept()).Throw(true);
	}

	// Prints file name followed by message 
	//	@msg: message
	//	@eol: if true then carriage should be return
	inline void PrintMsg(const string& msg, bool eol) const {
		Err(msg, FileNameToExcept()).Throw(false, eol);
	}

	// Throws exception occurred in the current reading line 
	//	@code: exception code
	inline void ThrowLineExcept(Err::eCode code) const {
		Err(_errCode=code, LineNumbToStr().c_str()).Throw();
	}

	// Throws exception occurred in the current reading line 
	//	@msg: exception message
	inline void ThrowLineExcept(const string& msg) const {
		Err(msg, LineNumbToStr().c_str()).Throw();
	}

	// Throws warning occurred in the current reading line 
	//	@msg: exception message
	//	@warnMsg: warning message
	inline void ThrowLineWarning(const string& msg, const string& warnMsg) const {
		Err(msg, LineNumbToStr().c_str()).Warning(warnMsg);
	}
	
	// Gets length of current line without EOL marker: only for single-line record!
	inline chrlen LineLength()	const { return RecordLength() - EOLSize(); }
};

#ifdef _FILE_WRITE

// 'TxtOutFile' represents TxtFile for writing
class TxtOutFile : public TxtFile
{
#ifdef _ISCHIP
	friend class ExtOutFile;
#endif

private:
	char	_delim;			// field delimiter in line
	// === line write buffer
	char*	_lineBuff;		// line write buffer; for writing mode only
	rowlen	_lineBuffLen;	// length of line write buffer in writing mode, otherwise 0
	rowlen	_lineBuffOffset;// current shift from the _buffLine; replaced by #define!!!
#ifdef _MULTITHREAD
	// === total counter of writed records
	ULONG*	_totalRecCnt;	// pointer to total counter of writed records; for clone only
#endif

	// Allocates memory for write line buffer with checking.
	//	return: true if successful
	bool CreateLineBuff(rowlen len);

	//	Adds delimiter on the given shift and increases current position.
	inline void LineAddDelim() { _lineBuff[_lineBuffOffset++] = _delim; }

	// Copies block of chars to the current position in the line write buffer,
	//	adds delimiter after string	and increases current position.
	//	@src: pointer to the block of chars
	//	@num: number of chars
	//	@addDelim: if true then adds delimiter and increases current position
	//	return: new current position
	rowlen LineAddChars(const char* src, size_t num, bool addDelim);

public:
	// Constructs an TxtOutFile instance: allocates I/O buffer, opens an assigned file.
	//	@fName: valid full name of assigned file
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	//	@printName: true if file name should be printed in the exception's message
	inline TxtOutFile(const string& fName, char delim = TAB,
		bool abortInvalid = true, bool printName = true) :
		_lineBuff(NULL), _lineBuffLen(0), _lineBuffOffset(0), _delim(delim),
		TxtFile(fName, WRITE, abortInvalid, printName)
		{
#ifdef _MULTITHREAD
			_totalRecCnt = &_recCnt;	// for atomic increment
#endif
		}		// line buffer will be created in SetLineBuff()

	// Writes nonempty buffer, deletes line buffer and closes file
	~TxtOutFile();

protected:
#ifdef _MULTITHREAD
	// Constructs a clone of an existing instance.
	// Clone is a copy of opened file with its own separate basic I/O and write line buffers.
	// Used for multithreading file recording
	//	@file: opened file which is cloned
	TxtOutFile	(const TxtOutFile& file);
#endif

protected:
	// Sets current position of the line write buffer.
	inline void LineSetOffset(rowlen offset)	{ _lineBuffOffset = offset; }

	// Copies the C string to the current position in the line write buffer,
	//	adds delimiter after string and increases current position.
	//	@str: C string to be copied
	//	@len: length of string to be copied
	//	@addDelim: if true then adds delimiter and increases current position
	//	return: new current position
	inline rowlen LineAddStr(const char* str, int len, bool addDelim=true)	{ 
		return LineAddChars(str, len, addDelim);
	}

	// Copies the string to the current position in the line write buffer,
	//	adds delimiter after string and increases current position.
	//	@str: string to be copied
	//	@addDelim: if true then adds delimiter and increases current position
	//	return: new current position
	inline rowlen LineAddStr(const string& str, bool addDelim=true)	{ 
		return LineAddChars(str.c_str(), str.length(), addDelim);
	}

	// Adds integral value to the current position of the line write buffer,
	//	adds delimiter after value and increases current position.
	//	@val: value to be set
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddInt(LLONG val, bool addDelim=true);

	// Adds two integral values separated by default delimiter to the current position 
	// of the line write buffer, adds delimiter after value and increases current position.
	//	@val1: first value to be set
	//	@val2: second value to be set
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddInts(ULONG val1, ULONG val2, bool addDelim=true);

	// Adds three integral values separated by default delimiter to the current position 
	// of the line write buffer, adds delimiter after value and increases current position.
	//	@val1: first value to be set
	//	@val2: second value to be set
	//	@val3: third value to be set
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddInts(ULONG val1, ULONG val2, ULONG val3, bool addDelim=true);

	// Adds floating point value to the current position of the line write buffer,
	//	adds delimiter after value and increases current position.
	//	@val: value to be set
	//	@ndigit: number of digits to generate
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddFloat(float val, BYTE ndigit, bool addDelim=true);

	// Adds line to the IO buffer (from 0 to the current position).
	//	@num: number of bytes from the beginning of line buffer to be added to file buffer
	//	@offset: start position for the next writing cycle
	//	@closeLine: if true then close line by EOL
	void LineToIOBuff(rowlen offset=0, bool closeLine=true);

	// Adds record to the IO buffer.
	// Generates exception if writing is fall.
	//	@src: record
	//	@len: length of the record
	//	@closeLine: if true then close line by EOL
	void RecordToIOBuff	(const char *src, UINT len, bool closeLine=true);

	// Adds string to IO buffer.
	inline void StrToIOBuff	(const string& str) { RecordToIOBuff(str.c_str(), str.length()); }

	// Adds commented string to the IO buffer.
	//	@str: recorded string
	inline void CommLineToIOBuff(const string& str) { StrToIOBuff("# " + str);	}

public:
	// Allocates line write buffer.
	//	@len: length of buffer
	//	@delim: delimiter between line fields
	void SetLineBuff(rowlen len);

	// Adds to line int and float values and adds line to the IO buff.
	//	@ndigit: number if float digits to write
	//void WriteLine(int val1, float val2, BYTE ndigit);

	// Adds to line int and float values and adds line to the IO buff.
	//	@ndigit: number if float digits to write
	void WriteLine(int val1, float val2, float val3, BYTE ndigit, const char*str);

	// Adds to line 2 int values and adds line to the IO buff.
	void WriteLine(ULONG val1, ULONG val2);

	// Adds to line C string and adds line to the IO buff.
	void WriteLine(const char* str);

	// Adds to line string and int value and adds line to the IO buff.
	void WriteLine(const string& str, int val);

	// Writes thread-safely current block to file.
	void Write	() const;
	
	// Adds content of another file (concatenates)
	//	return: true if successful
	//bool	AddFile(const string fName);
};
#endif	// _FILE_WRITE

#ifndef _FQSTATN
#ifndef _WIGREG

// 'Region' represents a simple region within nucleotides array (chromosome).
struct Region
{
	chrlen Start;	// start position of the region in standard chromosomal coordinates
	chrlen End;		// end position of the region in standard chromosomal coordinates

	inline Region(chrlen start=0, chrlen end=0) : Start(start), End(end) {}

	// Gets length of region.
	// The End is not included in the bases https://genome.ucsc.edu/FAQ/FAQformat.html#format1
	inline chrlen Length() const {	return End - Start; }

	inline chrlen Centre() const {	return Start + (Length()>>1); }

	// Initializes instance
	inline void Init(chrlen start, chrlen end) { Start = start; End = end; }

	inline bool operator==(const Region& r) const { return End == r.End && Start == r.Start; }

	// Returns true if this instance is invalid
	//inline bool Invalid() const { return Start >= End; }

	// Returns true if Region r is covered by this instance.
	//	@r: tested Region; should be sorted by start position
	inline bool Cover(const Region& r) const { return r.End <= End && r.Start >= Start; }

	// Returns true if Region r is adjoined with this instance.
	//	@r: tested Region; should be sorted by start position
	inline bool Adjoin(const Region& r) const { return r.Start == End; }

	// Returns true if Region r is crossed with this instance.
	//	@r: tested Region; should be sorted by start position
	inline bool Cross(const Region& r) const { return r.Start < End && r.End > Start; }

	// Compares two Regions by start position. For sorting a container.
	static inline bool CompareByStartPos(const Region& r1, const Region& r2) {
		return r1.Start < r2.Start;
	}

	// Extends Region with chrom length control.
	// If extended Region starts from negative, or ends after chrom length, it is fitted.
	//	@extLen: extension length in both directions
	//	@cLen: chrom length; if 0 then no check
	void Extend(chrlen extLen, chrlen cLen);

#ifdef DEBUG
	inline void Print() const { cout << Start << TAB << End << EOL; }
#endif
};

class Regions
/*
 * 'Regions' represents a container of defined regions within chromosome.
 * Defined in TxtFile.h since it's used in Fa class.
 */
{
private:
	vector<Region> _regions;

public:
	typedef vector<Region>::const_iterator Iter;

	// Iterator to region Begin position
	inline const Iter Begin()	const { return _regions.begin(); }

	// Iterator to region End position
	inline const Iter End()		const { return _regions.end(); }

	// Default (empty) constructor needed for adding to Chroms collection and to state common chroms
	inline Regions() {}
	
	// Copying constructor
	//inline Regions(const Regions& rgns) { _regions = rgns._regions;	}
	
	// Single region constructor
	inline Regions(chrlen start, chrlen end) { _regions.push_back(Region(start, end)); }
	
	// Gets count of regions.
	inline chrlen Count()	const { return chrlen(_regions.size()); }
	
	// Geta total length of regions.
	chrlen Length() const;
	
	// Gets first start position.
	inline chrlen FirstStart()	const { return _regions.front().Start; }
	
	// Gets last end position.
	inline chrlen LastEnd()		const { return _regions.back().End; }

	inline Regions& operator=(const Regions& rgn)
	{
		_regions = rgn._regions;
		return *this;
	}

	// Reserves container's capacity.
	//	@count: reserved number of regions. The real number may be differ.
	inline void Reserve(chrlen count) { _regions.reserve(count); }
	
	// Clears container.
	inline void Clear()	{ _regions.clear(); }

	// Adds gap beginning of current gap position.
	//	@gapStart: gap's start position
	//	@currGapStart: current gap start position holder
	//	@minGapLen: minimal length which defines gap as a real gap
	void AddGap(chrlen gapStart, chrlen currGapStart, chrlen minGapLen);
	
	inline const Region & operator[](chrlen ind) const { return _regions[ind]; }

	// Copies external Regions to this instance
	inline void Copy(const Regions &regns) { _regions = regns._regions;	}
	
	// Copies subReagions from external Regions.
	inline void Copy(const vector<Region>& regns, chrlen start, chrlen stop) {
		_regions = vector<Region>(regns.begin() + start, regns.begin() + stop + 1);
	}
#ifdef DEBUG
	void Print() const;
#endif
	friend class ShellRegions;

#if defined _DENPRO || defined _BIOCC

	// Returns an iterator referring to the past-the-end element, where end is external
	//	@curr_it: region's const iterator, from which the search is started
	//	@end: external pre-defined end coordinate
	Iter ExtEnd(Iter curr_it, chrlen end) const;

	// Initializes this instance by intersection of two Regions.
	void FillOverlap(const Regions &regn1, const Regions &regn2);

	// Initializes this instance by inverted external Regions.
	//	@regn: external Regions
	//	@masEnd: the maximum possible end-coordinate of region:
	//	the chromosome length in case of nucleotides sequance.
	void FillInvert(const Regions &regn, chrlen maxEnd);

	//// Initializes this instance by external Regions.
	//inline void Copy(const vector<Region>& regns) {	_regions = regns; }
	
	// Adds Region.
	void inline AddRegion(chrlen start, chrlen end)	{ _regions.push_back(Region(start, end)); }

protected:
	// Reads data from file @fname
	//	@fname: name of file
	//	return: written minimal gap length
	short Read(const string & fName);
	
	// Saves data to file
	//	@fname: name of file
	//	@minGapLen: minimal length which defines gap as a real gap
	void Write(const string & fName, short minGapLen) const;

#endif	// _DENPRO, _BIOCC
};

#endif	// _WIGREG

class TabFile : public TxtInFile
/*
 * Class 'TabFile' is a common class for reading/writing text files, 
 * each line of which is composed of number of fields splitted by TAB markers.
 * Only defined in constructor number of fields is processed.
 * Number of fields can be equal 1, in that case ordinary plain text can be processed.
 */
{
private:
	TabFilePar _params;
	char	*_currLine;		// current readed line; for reading only
	short	*_fieldPos;		// array of start positions in _currLine for each field
	USHORT	_lineSpecLen;	// length of line specifier; for reading only

	// Checks if filed valid and throws exception if not.
	//	@fInd: field index
	//	return: true if field is valid
	bool IsFieldValid	(BYTE fInd) const;

	// Initializes new instance.
	//	@mode: action mode (read, write, all)
	void Init(eAction mode);

public:
	// Creates new instance.
	//	@fName: name of file
	//	@minCntFields: obligatory number of fields separated by TAB; these fields are checked during initialization 
	//	@maxCntFields: maximum number of fields separated by TAB; these fields are checked by a call
	//	@mode: action mode (read, write, all)
	//	@lineSpec: specific substring on which each data line is beginning;
	//	other lines are interpreting as comments and would be skipped; for reading only
	//	@comment: char indicates that line is comment; for reading only
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	//	@rintName: true if file name should be printed in exception's message
	TabFile(
		const string& fName,
		eAction mode=READ,
		BYTE minCntFields=1,
		BYTE maxCntFields=1,
		char comment=HASH,
		const char* lineSpec=NULL,
		bool abortInvalid=true,
		bool printName=true
	) : _params(minCntFields, maxCntFields, comment, lineSpec),
		TxtInFile(fName, mode, 1, abortInvalid, printName)
	{	Init(mode); }

	// Creates new instance for reading
	//	@fName: name of file
	//	@params: TabFile params
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	//	@rintName: true if file name should be printed in exception's message
	TabFile(
		const string& fName,
		const TabFilePar& params,
		bool abortInvalid=true,
		bool printName=true
	) : _params(params),
		TxtInFile(fName, TxtFile::READ, 1, abortInvalid, printName)
	{	Init(TxtFile::READ); }


#if defined _ISCHIP && defined  _MULTITHREAD
	// Creates a clone of TabFile class.
	// Clone is a copy of opened file with its own buffer for writing only. Used for multithreading file recording.
	//  @file: opened file which is copied
	//  @threadNumb: number of thread
	inline TabFile(const TabFile& file) :
		_params(file._params),
		_fieldPos(file._fieldPos), 
		_lineSpecLen(file._lineSpecLen),
		TxtInFile(file/*, threadNumb*/) {}
#endif

	inline ~TabFile()	{ if( _fieldPos ) delete [] _fieldPos; }

	// Gets count of lines.
	inline ULONG Count() const { return RecordCount(); }

	// Returns true if invalid instance should completed by throwing exception.
	//	Used in WigMap (bioCC) only
	inline bool IsAborting() const	{ return IsFlag(ABORTING); }

	// Reads first line and set it as current.
	// Throws exception if invalid and aborting file is set
	//	@cntLines: returned estimated count of lines.
	//	It works properly only if lines are sorted by ascending, f.e. in sorted bed-files.
	//	return: current line
	const char*	GetFirstLine(ULONG *cntLines);

	// Reads next line and set it as current.
	//	return: current line.
	const char* GetLine();

	// Reads string by field's index from current line without check up.
	inline const char* StrField	(BYTE fInd)	const {	return _currLine + _fieldPos[fInd]; }

	// Reads string by field's index from current line with check up.
	const char* StrFieldValid	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? StrField(fInd) : NULL;
	}
	// Gets pointer to the short chrom name in current line without check up; for BED-file only
	const char* ChromName()	const {	return StrField(0) + strlen(Chrom::Abbr); }

	// Reads checked integer by field's index from current line without check up.
	inline int IntField	(BYTE fInd)	const {	return atoi(StrField(fInd)); }

	// Reads checked integer by field's index from current line with check up.
	int	 IntFieldValid	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? IntField(fInd) : vUNDEF;
	}
	// Reads float by field's index from current line without check up.
	inline float FloatField	(BYTE fInd)	const {	return float(atof(StrField(fInd))); }

	// Reads float by field's index from current line with check up.
	float FloatFieldValid	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? FloatField(fInd) : vUNDEF;
	}
	// Reads long by field's index from current line without check up.
	inline long LongField	(BYTE fInd)	const {	return long(atol(StrField(fInd))); }

	// Reads long by field's index from current line with check up.
	long LongFieldValid	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? LongField(fInd) : vUNDEF;
	}
};

#if !defined _WIGREG

class FaFile : public TxtInFile
/*
 * 'FaFile' supports reading/writing chromosomes in FA format.
 * Reading/writing based on TxtFile methods.
 * FaFile just adds some additional info support.
 */
{
public:
	// Fa file axtention
	static const string Ext;

	class Pocket	// keeps temporary variables needed for reading constructor only
	{
		chrlen	_countN,		// count of 'N' in chromosome
				_currGapStart,	// current gap start position
				_minGapLen,		// minimal length which defines gap as a real gap
				_cLen,			// length of chromosome
				_currPos;		// current readed nt position
		Regions& _defRgns;		// external chromosome's defined regions

		// Sets the initial state.
		void Clear();

		// Adds some amount of 'N' as a gap to assign defined regions
		//	@gapStart: gap's start position.
		//	@gapLen: gap's length; if 0, complete adding gaps
		void AddN(chrlen shiftGapStart, chrlen gapLen);

	public:
		inline Pocket(Regions& rgns, chrlen minGapLen)
			: _countN(0), _currGapStart(0), _minGapLen(minGapLen), _defRgns(rgns) {}

		// Gets total amount of 'N'
		inline chrlen CountN()		const { return _countN; }

		// Gets chromosome's length
		inline chrlen ChromLength() const { return _cLen; }

		// Closes adding gaps
		inline void CloseAddN() { _defRgns.AddGap(_cLen, _currGapStart, _minGapLen); }

		friend class FaFile;
	};

	// Opens an existing .fa fil and reads first line.
	//	@fName: full .fa file name
	//	@pocket: external temporary variables
	FaFile	(const string & fName, Pocket& pocket);

//#if defined _FILE_WRITE && defined DEBUG 
//	// Creates a new FaFile and writes header line.
//	FaFile	(const string & fName, const char *chrName);
//#endif

	// Gets current readed line.
	inline const char* Line()	const { return Record(); }
	
	// Reads line and set it as current with filling regions
	// Since method scans the line, it takes no overhead expenses to do this.
	//	@pocket: external temporary variables
	const char* GetLine(Pocket& pocket);

	// Reads line and set it as current without filling regions
	inline const char* GetLine() { return GetRecord(); }

//#ifdef _FILE_WRITE
	// Adds line to file
	//inline void AddLine(const char *src, BYTE len)	{ RecordToIOBuff(src, len); }
//#endif	// _FILE_WRITE
};

#endif	// _WIGREG
#else	// defined _FQSTATN

// 'FqFile' implements reading/writing file in FQ format.
class FqFile : public TxtInFile
{
	enum eLineLenIndex { HEADER1, READ, HEADER2, QUAL };

public:	
	// Creates new instance for reading by file name
	inline FqFile(const string& fileName)
		: TxtInFile(fileName, TxtFile::READ, 4, true, false) {}

	// Returns checked length of current readed Read.
	readlen ReadLength() const;
	
	// Gets checked Read from readed Sequence.
	const char* GetCurrRead() const;

	// Returns checked Sequence.
	const char*	GetSequence	();
	
	// Returns count of sequences.
	inline ULONG Count() const { return RecordCount(); }
};

#endif	// _FQSTATN
