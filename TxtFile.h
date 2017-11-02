#pragma once
#include "common.h"

#define INT_CAPACITY	10		// maximal number of digits in integer
// Number of basics file's reading|writing buffer blocks.
// Should be less than 2047 because of ULONG type of block size variable.
// Otherwise the behaviour is unpredictable.
#define NUMB_BLK 32
#define BASE_BLK_SIZE (2 * 1024 * 1024)	// basic block 2 Mb

typedef short rowlen;	// type: length of row in TxtFile

#define _buffLineOffset _readingLen

class TxtFile
/*
 * Basic class 'TxtFile' implements a fast buffered serial (stream) reading/writing text files 
 * (in Windows and Linux standart).
 * Supports reading/writing zipped (gz) files.
 * Optimised for huge files.
 * Restriction: the default size of buffer, setting as NUMB_BLK * BASE_BLK_SIZE, 
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
		READ,	// reads only existing file
		WRITE,	// creates file if it not exist and writes to it; file is cleared before
		ALL		// creates file if it not exist and reads it
	};

private:
	enum eFlag {			// signs of file
		CLONE		= 0x001,	// file is a clone
		CONSTIT		= 0x002,	// file is a constituent of an aggregate file
		ZIPPED		= 0x004,	// file is zipped
		ABORTING	= 0x008,	// invalid file should be completed by throwing exception; for Reading mode
		CRSET		= 0x010,	// CR symbol is present in the file; for Reading mode
		CRCHECKED	= 0x020,	// the presence of a symbol CR is checked; for Reading mode
		ENDREAD		= 0x040,	// last call of GetRecord() has returned NULL; for Reading mode
		PRNAME		= 0x080,	// print file name in exception's message; for Reading mode
		MTHREAD		= 0x100		// file in multithread mode: needs to be locked while writing
	};
	enum eBuff {		// signs of buffer; used in CreateBuffer() only
		BUFF_BASIC,		// basic (block) read|write buffer
		BUFF_LINE		// line write buffer
	};

	string	_fName;			// file's name
	LLONG	_fSize;			// the length of uncompressed file; for zipped file more than
							// 4'294'967'295 its unzipped length is unpredictable
	void *	_stream;		// FILE* (for unzipped file) or gzFile (for zipped file)
	mutable short _flag;	// bitwise storage for signs included in eFlag
	// === basic read|write buffer
	char *	_buff;			// basic accumulative read/write buffer
	UINT *	_linesLen;		// for Reading mode only: array of lengths of lines in a record
	UINT	_buffLen;		// the length of basic buffer
	mutable UINT _currRecPos;// start position of the last readed/writed record in current block
	ULONG	_cntRecords;	// counter of readed/writed records
	BYTE	_cntRecLines;	// number of lines in a record
	UINT	_recLen;		// for Reading mode only: the length of record with EOL marker
	// === basic & line buffer common use
	UINT	_readingLen;	// while file reading: number of actually readed chars in block
							// while line writing: current shift from the _buffLine
	// === line write buffer
	char*	_buffLine;		// line write buffer; for writing mode only
	rowlen	_buffLineLen;	// length of line write buffer in writing mode, otherwise 0
	//rowlen	_buffLineOffset;// current shift from the _buffLine; replacement by #define!!!
protected:
	char	_delim;

private:
	inline void RaiseFlag	(eFlag flag) const	{ _flag |= flag; }
	inline void SetFlag	(eFlag flag, bool val)	{ val ? _flag |= flag : _flag &= ~flag; }
	inline bool IsFlag(eFlag flag)	const	{ return (_flag & flag) != 0; }	// != 0 to avoid warning C4800
	inline bool IsZipped()				const	{ return IsFlag(ZIPPED); }

	// Raises ENDREAD sign an return NULL
	inline char* ReadingEnded()		  { RaiseFlag(ENDREAD); return NULL; }

	// Reads next block.
	//	@offset: shift of start reading position
	//	return: 1 if file is not finished; 0 if it is finished; -1 if unsuccess reading
	int ReadBlock(const UINT offset);

	// Initializes instance variables and opens a file with setting a proper error code.
	//	@fName: valid full name of file
	//	@mode: opening mode
	//	@file: clonable file or NULL
	//	return: true is success, otherwise false.
	bool SetBasic(const string& fName, eAction mode, void* file);
	
	// Allocates memory for file read/write buffer or write line buffer with checking.
	//	return: true if successful
	bool CreateBuffer(eBuff buffType);

	// Returns file name or empty string depends on if name is printing
	const string& FileNameToExcept() const { return IsFlag(PRNAME) ? _fName : strEmpty; }

	// Gets string containing file name and current record number.
	//inline const string RecordNumbToStr() const
	//{ return LineNumbToStr(0); }

	// Gets number of line.
	//	@lineInd: index of line in a record
	inline ULONG LineNumber	(BYTE lineInd) const
	{ return (_cntRecords-1)*_cntRecLines + lineInd + 1; }

protected:
	mutable Err::eCode	_errCode;
	
	// General constructor, always is called first.
	//	@fName: valid full name of file
	//	@mode: opening mode
	//	@cntRecLines: number of lines in a record
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	//	@rintName: true if file name should be printed in the exception's message
	TxtFile(const string& fName, eAction mode, BYTE cntRecLines, bool abortInvalid=true, bool printName=true);

#ifdef _MULTITHREAD
	// Creates new instance with read buffer belonges to aggregated file: constructor for concatenating.
	//	For reading only.
	//	@fName: valid full name of file
	//  aggrFile: file that is aggregated
	//TxtFile	(const string& fName, const TxtFile& aggrFile);
	
	// Creates a clone of existing instance.
	// Clone is a copy of opened file with its own separate write basic and write line buffers.
	// Used for multithreading file recording
	//	@file: opened file which is cloned
	//	@threadNumb: number of thread
	//	Exception: file_error
	TxtFile	(const TxtFile& file, threadnumb threadNumb);
#endif
	
	~TxtFile();

	// Gets the amount of characters corresponded to EOL
	// In file created in Window EOL matches '\r\n', in file created in Linux EOL matches '\n',
	//	return: in Windows always 1, in Linux: 2 for file created in Windows, 1 for file created in Linux
	inline BYTE EOLSize() const	{ return (BYTE)IsFlag(CRSET) + 1; }

	// Returns true if instance is a clone.
	inline bool IsClone() const { return IsFlag(CLONE); }	// return _isClone;

	// Sets error code and throws exception if it is allowed.
	void SetError(Err::eCode errCode) const;


	// Gets length of current reading record including EOL marker.
	inline UINT	RecordLength() const	{ return _recLen; }
	
	// Gets length of line without EOL marker.
	//	@lineInd: index of line in a record
	//inline UINT	LineLengthByInd	(BYTE lineInd) const { return _linesLen[lineInd] - EOLSize(); }
	
	// Gets string containing file name and current line number.
	//	@lineInd: index of line in a record; if 0, then first line
	const string LineNumbToStr(BYTE lineInd = 0) const {
		return (IsFlag(PRNAME) ? (_fName + SepSCl) : SepSCl)
			+ "line " + NSTR(LineNumber(lineInd));
	}

	// Gets current record.
	inline const char* Record() const
	{ return IsFlag(ENDREAD) ? NULL : _buff + _currRecPos - _recLen; }

	// Gets next record.
	//	return: point to the next record in a buffer.
	inline const char* NextRecord() const	{ return _buff + _currRecPos; }

	// Sets _currLinePos to the beginning of next non-empty line inread/write buffer.
	//	@counterN: if not NULL, adds to counterN the number of 'N' in a record. Used in Fa() only.
	//	@posTab: if not NULL, sets TABs positions in line to this array
	//	@cntTabs: if 'posTab' is not NULL, the length of 'posTab' array
	//	return: point to line or NULL if no more lines
	char*	GetRecord(chrlen* const counterN=NULL, short* const posTab=NULL, BYTE cntTabs=0);
	
#ifdef _FQSTATN
	// Throw exception if no record is readed.
	void CheckGettingRecord() const { 
		if( !IsFlag(CRCHECKED) )
			Err("attempt to get record's info without reading record", _fName.c_str()).Throw();
	}
#endif

#ifdef _FILE_WRITE

private:
	//	Adds delimiter on the given shift and increases current position.
	//	return: new current position
	rowlen LineAddDelim();

	// Copies block of chars to the current position in the line write buffer,
	//	adds delimiter after string	and increases current position.
	//	@src: pointer to the block of chars
	//	@num: number of chars
	//	return: new current position
	rowlen LineAddChars(const char* src, size_t num);

	// Copies block of chars before the current position in the line write buffer,
	//	adds delimiter before string and decreases current position.
	//	@src:  pointer to the block of chars
	//	@num: number of chars
	void LineAddCharsBack(const char* src, size_t num);

protected:
	// Returns a pointer to the line write buffer at current position.
	inline char* LineCurrPosBuf() const	{ return _buffLine + _buffLineOffset;	}

	// Returns a pointer to the line write buffer at given position.
	//	@shift: given relative pointer position 
	inline char* LineCurrPosBuf(rowlen shift) const	{ return _buffLine + shift;	}

	// Creates line write buffer.
	//	@len: length of buffer
	//	@delim: delimiter between line fields
	void SetWriteBuffer(rowlen len, char delim);

	// Creates line write buffer and copies its content from mate file
	//	@file: file from which the buffer sould be copied
	//	@createLineBuff: true if line write buffer should be created (for mates), otherwise false (for clones)
	void CopyWriteBuffer(const TxtFile& file, bool createLineBuff);

	// Sets current position of the line write buffer.
	inline void LineSetOffset(rowlen offset)	{ _buffLineOffset = offset; }

	// Sets the same bytes in the line write buffer to the specified value.
	//	@offset: shift of buffer start position of filling bytes
	//	@val: value to be set
	//	@len: number of bytes to be set to the value, or the rest of buffer by default
	inline void LineFill(rowlen offset, char val, rowlen len=0)
	{ memset(_buffLine + offset, val, len ? len : (_buffLineLen-offset)); }

	// Sets some bytes in the line write buffer to the specified value.
	//	@offset: shift of buffer start position of filling bytes
	//	@src: bytes to be set
	//	@len: number of bytes to be set to the value, or the rest of buffer by default
	inline void LineFill(rowlen offset, const char* src, rowlen len)
	{ memcpy(_buffLine + offset, src, len); }


	// Moves current position of the line write buffer by shift.
	//inline void SlipOffset(int shift)	{ _buffLineOffset += shift; }

	// Decreases current position of the line write buffer by one.
	inline void LineDecreaseOffset() {	_buffLineOffset--; }

	// Increases current position of the line write buffer by one.
	//inline void LineIncreaseOffset() { _buffLineOffset++; }

	// Adds byte to the current position in the line write buffer, adds delimiter after byte
	//	and increases current position.
	//	@chr: value to be set
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddChar(char chr, bool addDelim=false);
	
	// Adds byte before the current position of the line write buffer, adds delimiter before byte
	//	and decreases current position.
	//	@chr: value to be set
	//	@addDelim: if true then adds delimiter before and decreases current position
	void LineAddCharBack(char chr, bool addDelim=false);

	// Adds floating point value to the current position of the line write buffer,
	//	adds delimiter after value and increases current position.
	//	@val: value to be set
	//	@ndigit: number of digits to generate
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddFloat(double val, BYTE ndigit, bool addDelim=true);

	// Adds integral value to the current position of the line write buffer,
	//	adds delimiter after value and increases current position.
	//	@val: value to be set
	//	@ndigit: number of digits to generate
	inline void LineAddInt(LLONG val, bool addDelim=true) {
		LineAddFloat(double(val), DigitsCount(val), addDelim);
	}

	// Copies block of chars to the current position of the line write buffer.
	//	@src: pointer to the block of chars
	//	@num: number of chars to be copied
	inline void LineCopyChars(const char* src, size_t num) const {
		memcpy(_buffLine + _buffLineOffset, src, num);
	}

	// Copies the C string to the current position in the line write buffer,
	//	adds delimiter after string and increases current position.
	//	@str: C string to be copied
	//	return: new current position
	//inline rowlen LineAddStr(const char* str)	{ return LineAddStr(str, strlen(str)); }

	// Copies the string to the current position in the line write buffer,
	//	adds delimiter after string and increases current position.
	//	@str: string to be copied
	//	return: new current position
	inline rowlen LineAddStr(const string& str)	{ return LineAddChars(str.c_str(), str.length()); }

	// Copies the string before the current position in the line write buffer,
	//	adds delimiter before string and decreases current position.
	//	@str: string to be copied
	inline void LineAddStrBack(const string& str) {	LineAddCharsBack(str.c_str(), str.length()); }

	// Adds first part of the line write buffer (from 0 to the current position)
	//	to the file write buffer.
	//	@num: number of bytes from the beginning of line buffer to be added to file buffer
	//	@offset: start shift for the next writing cycle
	//	@closeLine: if true then close line by EOL
	void LineToBuffer(rowlen offset=0, bool closeLine=true);

	// Adds last part of the line write buffer (from current position to the end)
	//	to the file write buffer.
	inline void LineBackToBuffer()
	{ AddRecord(_buffLine + _buffLineOffset, _buffLineLen - _buffLineOffset); }

	// Adds record to the file write buffer.
	// Generates exception if writing is fall.
	//	@src: record
	//	@len: length of record
	//	@closeLine: if true then close line by EOL
	void AddRecord	(const char *src, UINT len, bool closeLine=true);

#endif	//  _FILE_WRITE

public:
	// Returns true if instance is valid.
	inline bool IsGood() const	{ return _errCode == Err::NONE; }

	// Returns true if instance is invalid.
	inline bool IsBad()	const	{ return _errCode != Err::NONE; }

	// Returns current error code.
	inline Err::eCode ErrCode() const	{ return _errCode; }

	// Throws exception
	//	@msg: exception message
	//	@isExc: if true then throw exception, otherwise warning
	inline void ThrowExcept(const string& msg, bool isExc = true) const {
		Err(msg, FileNameToExcept()).Throw(isExc);
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

	// Throws exception occurred in the current reading line 
	//	@msg: exception message
	inline void ThrowLineWarning(const string& msg, const string& warnMsg) const {
		Err(msg, LineNumbToStr().c_str()).Warning(warnMsg);
	}

	
	// Gets size of uncompressed file,
	// or size of compressed file if its uncompressed length is more than UINT_MAX.
	// So it can be used for estimatation only.
	inline LLONG Length() const	{ return _fSize; }

	// Gets file name.
	inline const string& FileName() const	{ return _fName; }

	// Gets number of readed/writed records.
	inline ULONG RecordCount() const	{ return _cntRecords; }

	// Gets length of current line without EOL marker: only for single record!
	inline chrlen LineLength()	const { return RecordLength() - EOLSize(); }
	//inline UINT LineLength() const { return LineLengthByInd(0); }

	// Returns true if invalid instance should completed by throwing exception.
	inline bool IsAborting() const	{ return IsFlag(ABORTING); }

#ifdef _FILE_WRITE
	
	// Writes thread-safely current block to file.
	void Write() const;
	
	// Adds content of another file (concatenates)
	//	return: true if successful
	//bool	AddFile(const string fName);

#endif	
};

#ifdef _FILE_WRITE

// 'LineFile' provides methods for writing text files by lines.
class LineFile : public TxtFile
{
public:
	// Creates an instance by file name
	//	@fname: file name
	//	@delim: line fields delimiter
	LineFile(const string& fname, char delim) : TxtFile(fname, TxtFile::WRITE, 1)
	{ _delim = delim; }

	// Initializes line write buffer.
	inline void BeginWrite(rowlen len) { SetWriteBuffer(len, _delim); }

	// Adds to line int and float values and writes line.
	void WriteLine(int val1, float val2, BYTE ndigit);

	// Adds to line 2 int values and writes line.
	void WriteLine(int val1, int val2);

	// Adds to line string and int value and writes line.
	void WriteLine(const string& str, int val);
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

	// Gets length of region
	inline chrlen Length() const {	return End - Start + 1; }

	// Initializes instance
	void Init(chrlen start, chrlen end) { Start = start; End = end; }

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
	//	@cLen: chrom length
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

	inline const Iter Begin()	const { return _regions.begin(); }
	inline const Iter End()		const { return _regions.end(); }

	// Default (empty) constructor needed for adding to Chroms collection and to state common chroms
	inline Regions() {}
	
	// Copying constructor
	inline Regions(const Regions& rgns) { _regions = rgns._regions;	}
	
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

#ifdef DEBUG
	void Print() const;
#endif
	friend class ShellRegions;

#if defined _DENPRO || defined _BIOCC

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

	// Copies external Regions to this instance
	void inline Copy(const Regions &regions) { _regions = regions._regions; }
	
protected:
	// Reads data from file @fname
	//	@fname: name of file
	//	return: written minimal gap length
	short Read(const string & fName);
	
	// Saves data to file
	//	@fname: name of file
	//	@minGapLen: minimal length which defines gap as a real gap
	void Write(const string & fName, short minGapLen) const;

public:
	// Copies subReagions from external Regions.
	inline void Copy(const vector<Region>& regns, chrlen start, chrlen stop) {
		_regions = vector<Region>(regns.begin() + start, regns.begin() + stop + 1);
	}
#endif	// _DENPRO, _BIOCC
};

#endif	// _WIGREG

class TabFile : public TxtFile
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
	bool	_checkFieldCnt;	// true if fields count should be checked; for reading only

	// Checks if filed valid and throws exception if not.
	//	@fInd: field index
	//	return: true if field is valid
	bool IsFieldValid	(BYTE fInd) const;

	// Reads string by field's index from current line without control
	//	@fInd: field index
	inline const char* SField(BYTE fInd) const { return _currLine + _fieldPos[fInd]; }

	// Initializes new instance.
	//	@mode: action mode (read, write, all)
	void Init(eAction mode);

public:
	// Creates new instance.
	//	@fName: name of file
	//	@minCntFields: obligatory number of fields separated by TAB
	//	@maxCntFields: maximum number of fields separated by TAB
	//	@mode: action mode (read, write, all)
	//	@lineSpec: specific substring on which each data line is beginning;
	//	other lines are interpreting as comments and would be skipped; for reading only
	//	@comment: char indicates that line is comment; for reading only
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	//	@rintName: true if file name should be printed in exception's message
	//	@checkFieldCnt: true if fields count should be checked; for reading only
	TabFile(
		const string& fName,
		eAction mode=READ,
		BYTE minCntFields=1,
		BYTE maxCntFields=1,
		char comment=HASH,
		const char* lineSpec=NULL,
		bool abortInvalid=true,
		bool printName=true,
		bool checkFieldCnt=true
	) : _params(minCntFields, (maxCntFields==1 ? minCntFields : maxCntFields) + 1, comment, lineSpec),
		_checkFieldCnt(checkFieldCnt),
		TxtFile(fName, mode, 1, abortInvalid, printName)
	{	Init(mode); }

	// Creates new instance for reading
	//	@fName: name of file
	//	@minCntFields: obligatory number of fields separated by TAB
	//	@maxCntFields: maximum number of fields separated by TAB
	//	@mode: action mode (read, write, all)
	//	@lineSpec: specific substring on which each data line is beginning;
	//	other lines are interpreting as comments and would be skipped; for reading only
	//	@comment: char indicates that line is comment; for reading only
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	//	@rintName: true if file name should be printed in exception's message
	//	@checkFieldCnt: true if fields count should be checked; for reading only
	TabFile(
		const string& fName,
		const TabFilePar& params,
		bool abortInvalid=true,
		bool printName=true,
		bool checkFieldCnt=true
	) : _params(params),
		_checkFieldCnt(checkFieldCnt),
		TxtFile(fName, TxtFile::READ, 1, abortInvalid, printName)
	{	Init(TxtFile::READ); }


#if defined _ISCHIP && defined  _MULTITHREAD
	// Creates a clone of TabFile class.
	// Clone is a copy of opened file with its own buffer for writing only. Used for multithreading file recording.
	//  @file: opened file which is copied
	//  @threadNumb: number of thread
	inline TabFile(const TabFile& file, threadnumb threadNumb) :
		_params(file._params),
		_fieldPos(file._fieldPos), 
		_lineSpecLen(file._lineSpecLen),
		_checkFieldCnt(file._checkFieldCnt),
		TxtFile(file, threadNumb) {}
#endif

	inline ~TabFile()	{ if( _fieldPos ) delete [] _fieldPos; }

	// Gets count of lines.
	inline ULONG Count() const { return RecordCount(); }

	// Reads first line and set it as current.
	// Throws exception if invalid and aborting file is set
	//	@cntLines: returned estimated count of lines.
	//	It works properly only if lines are sorted by ascending, f.e. in sorted bed-files.
	//	return: current line
	const char*	GetFirstLine(ULONG *cntLines);

	// Reads next line and set it as current.
	//	return: current line.
	const char* GetLine();

	// Reads string by field's index from current line.
	const char* StrField	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? SField(fInd) : NULL;
	}
	// Gets pointer to the short chrom name in current line; for BED-file only
	const char* ChromName()	const {
		return IsFieldValid(0) ? (SField(0) + strlen(Chrom::Abbr)) : NULL;
	}
	// Reads integer by field's index from current line.
	int	 IntField	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? atoi(SField(fInd)) : vUNDEF;
	}
	// Reads float by field's index from current line.
	float FloatField	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? float(atof(SField(fInd))) : vUNDEF;
	}
	// Reads long by field's index from current line.
	long LongField	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? long(atol(SField(fInd))) : vUNDEF;
	}
};

#if !defined _WIGREG

class FaFile : public TxtFile
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

#if defined _FILE_WRITE && defined DEBUG 
	// Creates a new FaFile and writes header line.
	FaFile	(const string & fName, const char *chrName);
#endif

	// Gets current readed line.
	inline const char* Line()	const { return Record(); }
	
	// Reads line and set it as current with filling regions
	// Since method scans the line, it takes no overhead expenses to do this.
	//	@pocket: external temporary variables
	const char* GetLine(Pocket& pocket);

	// Reads line and set it as current without filling regions
	inline const char* GetLine() { return GetRecord(); }

#ifdef _FILE_WRITE
	// Adds line to file
	inline void AddLine(const char *src, BYTE len)	{ AddRecord(src, len); }
#endif	// _FILE_WRITE
};

#endif	// _WIGREG
#endif	// _FQSTATN

#if defined _ISCHIP || defined _FQSTATN
// 'FqFile' implements reading/writing file in FQ format.
class FqFile : public TxtFile
{
#ifdef _FQSTATN

private:
	enum eLineLenIndex { HEADER1, READ, HEADER2, QUAL };

public:	
	// Creates new instance for reading by file name
	inline FqFile(const string& fileName)
		: TxtFile(fileName, TxtFile::READ, 4, true, false) {}

	// Returns checked length of current readed Read.
	readlen ReadLength() const;
	
	// Gets checked Read from readed Sequence.
	const char* GetCurrRead() const;

	// Returns checked Sequence.
	char*	GetSequence	();
	
	// Returns count of sequences.
	inline ULONG Count() const { return RecordCount(); }	// << PairEnd; }	// multiply by 2 in case pf
#endif	// _FQSTATN

#ifdef _FILE_WRITE

private:
	static rowlen ReadStartPos;		// constant Read field start position
	//static const char*	QualPattern;	// Read quality pattern

	// Presets line write buffer
	//	@rQualPatt: Read quality pattern, or NULL
	void InitBuffer(const char* rQualPatt);

public:
	// Creates new instance for writing
	//	@fName: file name without extention
	//	@zipExt: zip extention if file should be zipped; otherwise empty string
	//	@mateNumb: mate number: 1 or 2, nothing by default
	inline FqFile(const string& fName, bool isZip, BYTE mateNumb = 0)
		: TxtFile(fName + (mateNumb ? (USCORE + BSTR(mateNumb)) : strEmpty) + FT::RealExt(FT::FQ, isZip),
			TxtFile::WRITE, 4) {}

#ifdef _MULTITHREAD
	// Creates a clone of existed file for writing only.
	//	@threadNumb: used for forming write buffer with differ sizes for anisochronous writing
	inline FqFile(const FqFile& file, threadnumb threadNumb) : TxtFile(file, threadNumb)
	{ CopyWriteBuffer(file, false); }
#endif

	// Creates and initializes line write buffer.
	//	@rQualPatt: Read quality pattern, or NULL
	void InitToWrite (const char* rQualPatt);

	// Initializes line write buffer from mate file
	void inline InitToWrite (const FqFile& file) { CopyWriteBuffer(file, true); }

	// Forms Read from fragment and adds it to the file.
	//	@rName: Read's name
	//	@read: valid read
	//	@reverse: if true then complement added read 
	void AddRead	(const string& rName, const char* read, bool reverse);

	// Gets Read qualuty pattern
	const inline string QualPatt() const { return string(LineCurrPosBuf(ReadStartPos), Read::Len); }

#endif	// _FILE_WRITE
};
#endif	// _ISCHIP || _FQSTATN


