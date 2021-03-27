/**********************************************************
TxtFile.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 23.03.2021
-------------------------
Provides read|write text file functionality
***********************************************************/

#pragma once
#include "common.h"

// Number of basics file's reading|writing buffer blocks.
// Should be less than 2047 because of ULONG type of block size variable.
// Otherwise the behaviour is unpredictable.

typedef short rowlen;	// type: length of row in TxtFile

// 'TabFilePar' keeps basic parameters for TabFile
struct TabFilePar 
{
	static const BYTE BGLnLen;	// predefined BedGraph line length in purpose of reserving container size
	static const BYTE WvsLnLen;	// predefined WIG variable step line length in purpose of reserving container size
	static const BYTE WfsLnLen;	// predefined WIG fixed step line length in purpose of reserving container size

	const BYTE	MinFieldCnt;	// minimum number of feilds in file line; checked during initialization
	const BYTE	MaxFieldCnt;	// maximum possible number of feilds in data line; checked by a call
	const BYTE	LineLen;		// average line length; to reserve container size
	const char	Comment;		// char indicates that line is comment
	const char* LineSpec;		// substring on which each data line is beginning

	inline TabFilePar() : MinFieldCnt(0), MaxFieldCnt(0), LineLen(0), Comment(cNULL), LineSpec(NULL) {}

	TabFilePar(BYTE minTabCnt, BYTE maxTabCnt, BYTE lineLen = 0, char comm = HASH, const char* lSpec = NULL) :
		MinFieldCnt(minTabCnt),
		MaxFieldCnt(maxTabCnt<minTabCnt ? minTabCnt : maxTabCnt),
		LineLen(lineLen),
		Comment(comm),
		LineSpec(lSpec)
		{}

	inline USHORT LineSpecLen() const { return USHORT(LineSpec ? strlen(LineSpec) : 0); }
};

// 'File Type' implements bioinformatics file type routines 
static class FT
{
	struct fTypeAttr {
		const char*	Extens;			// file extension
		const string Item;			// item title
		const string ItemPl;		// item title in plural
		Mutex::eType MtxType;		// type of Mutex lock
		TabFilePar FileParam;		// TabFile parameters, defined feilds
	};
	static const char* bedExt;
	static const char* wigExt;
	static const string Interval;
	static const string Intervals;
	static const string Read;
	static const string Reads;
	static const fTypeAttr	TypeAttrs[];
	static const BYTE	Count;

public:
	// bioinfo file types
	enum class eType {
		UNDEF,	// undefined type
		BED,	// ordinary bed
		ABED,	// alignment bed
		SAM,	// sam
		BAM,	// bam
		BGRAPH,	// bedgraph
		WIG_VAR,// wiggle_0 variable step
		WIG_FIX,// wiggle_0 fixed step
		FQ,		// fastQ
		FA,		// fasta
		CSIZE,	// chrom sizes
		RGN,	// my own region type
		DIST,	// my own distribution type
#ifdef _ISCHIP
		INI		// isChIP ini file type
#endif
	};

	// Gets file type
	//	@fName: file name (with case insensitive extension)
	//	@isABED: if true then returns ABED in case of .bed extention
	static eType GetType(const char* fName, bool isABED = false);

	// Gets file extension, beginning at DOT and adding .gz if needed
	//	@t: file type
	//	@isZip: true if add ".gz"
	static const string Ext(eType t, bool isZip = false);

	// Gets an item's title
	//	@t: file type
	//	@pl: true if plural form
	static const string& ItemTitle(eType t, bool pl = true)
	{ return pl ? TypeAttrs[int(t)].ItemPl : TypeAttrs[int(t)].Item; }

	// Gets TabFile params
	//	@t: file type
	static const TabFilePar& FileParams(eType t) { return TypeAttrs[int(t)].FileParam; }

	// Get Mutex type by file type
	inline static Mutex::eType MutexType(eType t) { return TypeAttrs[int(t)].MtxType; }

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
	enum class eAction { 
		READ,		// reads only existing file
		WRITE,		// creates file if it not exist and writes to it; file is cleared before
		READ_ANY	// creates file if it not exist and reads it
	};

protected:
	enum eFlag {			// signs of file
		// The first two right bits are reserved for storing the length of the LF marker: 1 or 2 
		// The first bit is set to 1 in the constructor.
		// If CR symbol is found at the end of line,
		//	the second bit is raised to 1, the first turn down to 0
		F_CR		= 0x01,	// sign of Curent Return ('\r')
		//EOLSZ		= 0x03,	// mask for the 2 first bits
		EOLCHECKED	= 0x02,	// the presence of a symbol CR is checked; for Reading mode
		ZIPPED		= 0x04,	// file is zipped
		ABORTING	= 0x08,	// invalid file should be completed by throwing exception; for Reading mode
		ENDREAD		= 0x10,	// last call of GetNextRecord() has returned NULL; for Reading mode
		PRNAME		= 0x20,	// print file name in exception's message; for Reading mode
		MTHREAD		= 0x40,	// file in multithread mode: needs to be locked while writing
		CLONE		= 0x80,	// file is a clone
	};

private:
	const static int BlockSize = 4 * 1024 * 1024;	// 4 Mb

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
	inline bool IsFlag	(eFlag f)	const	{ return _flag & f; }
	inline bool IsZipped()			const	{ return IsFlag(ZIPPED); }

	// *** 3 methods used by TxtInFile only

	// Gets the number of characters corresponded to LF
	// In Windows LF matches '\r\n', in Linux LF matches '\n',
	//	return: in Windows always 1, in Linux: 2 for file created in Windows, 1 for file created in Linux
	//inline UINT EOLSize() const	{ return _flag & EOLSZ; }
	inline UINT EOLSize() const	{ 	return 2 - (_flag & F_CR); }

	// Returns true if LF size is not defined
	inline bool IsEOLundef() const { return !bool(_flag & EOLCHECKED); }

	// Establishes the presence of CR symbol at the end of line.
	//	@c: if c is CR then the second bit is raised to 1, the first turn down to 0,
	//	so the value return by EOLSZ mask is 2, otherwise remains in state 1
	inline void SetEOL(char c)	{ if(c==CR)	 RaiseFlag(F_CR); RaiseFlag(EOLCHECKED); }

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

	void operator=(const TxtFile&); //assignment prohibition

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
	//	@errCode: error code
	//	@senderSpec: sender specificator (string added to file name)
	void SetError(Err::eCode errCode, const string& senderSpec = strEmpty, const string& spec = strEmpty) const;

	// Returns true if instance is a clone.
	inline bool IsClone() const { return IsFlag(CLONE); }

	// Returns true if instance is valid.
	inline bool IsGood() const	{ return _errCode == Err::NONE; }

	// Returns true if instance is invalid.
	inline bool IsBad()	const	{ return _errCode != Err::NONE; }

	// Returns current error code.
	//inline Err::eCode ErrCode() const	{ return _errCode; }

public:
	// Gets conditional file name: name if it's printable, otherwise empty string.
	inline const string& CondFileName() const { return IsFlag(PRNAME) ? _fName : strEmpty; }

public:
	// Gets file name.
	inline const string& FileName() const { return _fName; }

	// Gets size of uncompressed file,
	// or size of compressed file if its uncompressed length is more than UINT_MAX.
	// So it can be used for estimatation only.
	inline LLONG Length() const	{ return _fSize; }

	// Gets number of readed/writed records.
	ULONG RecordCount() const { return _recCnt; }

	// Throws exception
	//	@msg: exception message
	//	@isExc: if true then throw exception, otherwise warning
	inline void ThrowExcept(const string& msg, bool isExc = true) const {
		Err(msg, CondFileName()).Throw(isExc);
	}
};

// 'TxtInFile' represents TxtFile for reading.
//	Ñorrectly reads records even if the latter does not end with LF (CR,LF) character (s)
class TxtInFile : public TxtFile
{
private:
	UINT	_recLen;		// the length of record with LF marker
	UINT	_readedLen;		// number of actually readed chars in block
	UINT *	_linesLen;		// array of lengths of lines in a record
	BYTE	_recLineCnt;	// number of lines in a record

	// Raises ENDREAD sign an return NULL
	//inline char* ReadingEnded()		  { RaiseFlag(ENDREAD); return NULL; }

	// Returns record without control
	inline char* RealRecord() const { return _buff + _currRecPos - _recLen; }

	// Reads next block.
	//	@offset: shift of start reading position
	//	return: 0 if file is finished; -1 if unsuccess reading; otherwhise number of readed chars
	int ReadBlock(const UINT offset);

	// Returns true if block is complete, move remainder to the top
	bool CompleteBlock(UINT currLinePos, UINT blankLineCnt);

	// Fills I/O buffer with 0, beginning from @offset position
	//inline void ClearBuff(UINT offset = 0) { memset(_buff + offset, 0, _buffLen - offset); }

	// Gets string containing file name and current record number.
	//inline const string RecordNumbToStr() const { return LineNumbToStr(0); }

	// Gets number of line.
	//	@lineInd: index of line in a record
	inline ULONG LineNumber	(BYTE lineInd) const { return (_recCnt-1)*_recLineCnt + lineInd + 1; }

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
public:
	// Gets length of current reading record including LF marker.
	//	Returns 0 after RollBackLastRecord() invoke.
	inline UINT	RecordLength() const	{ return _recLen; }
protected:
	// Gets current record.
	inline const char* Record() const	{ return IsFlag(ENDREAD) ? NULL : RealRecord(); }

	// Reads record
	//	return: pointer to line or NULL if no more lines
	const char*	GetNextRecord();

	// Reads N-controlled record
	//	@counterN: counter of 'N'
	//	return: pointer to line or NULL if no more lines
	const char*	GetNextRecord(chrlen& counterN);

	// Reads tab-controlled record
	//	@tabPos: TAB's positions array that should be filled
	//	@cntTabs: maximum number of TABS in TAB's positions array
	//	return: pointer to line or NULL if no more lines
	char*	GetNextRecord(short* const tabPos, const BYTE tabCnt);

#if defined _CALLDIST || defined _FQSTATN
	// Gets next record
	//	return: point to the next record in a buffer.
	inline const char* NextRecord() const { return _buff + _currRecPos; }

	// Gets length of line.
	//	@lineInd: index of line in a record
	//	@withoutEOL: if true then without LF marker length, otherwise without it
	inline UINT	LineLengthByInd	(BYTE lineInd, bool withoutEOL = true) const { 
		return _linesLen[lineInd] - withoutEOL * EOLSize();
	}

	// Throw exception if no record is readed.
	void CheckGettingRecord() const { 
		if( IsEOLundef() )
			Err("attempt to get record's info without reading record", FileName().c_str()).Throw();
	}
#endif

	// Returns the read pointer to the beginning of the last read line and decreases line counter.
	//	Zeroes length of current reading record!
	//	@sep: field separator character
	void RollBackLastRecord(char sep);

public:
	// Gets current readed line.
	inline const char* Line()	const { return Record(); }

	// Gets string containing file name and current line number.
	//	@code: code of error occurs
	//	@lineInd: index of line in a record; if 0, then first line
	const string LineNumbToStr(Err::eCode code = Err::EMPTY, BYTE lineInd = 0) const;

	// Throws exception occurred in the current reading line 
	//	@code: exception code
	//inline void ThrowLineExcept(Err::eCode code) const {
	//	Err(_errCode=code, LineNumbToStr().c_str()).Throw();
	//}

	// Throws exception occurred in the current reading line 
	//	@msg: exception message
	//inline void ThrowLineExcept(const string& msg) const {
	//	Err(msg, LineNumbToStr().c_str()).Throw();
	//}

	// Throws warning occurred in the current reading line 
	//	@msg: exception message
	//	@warnMsg: warning message
	//inline void ThrowLineWarning(const string& msg, const string& warnMsg) const {
	//	Err(msg, LineNumbToStr().c_str()).Warning(warnMsg);
	//}
	
	// Gets length of current line without LF marker: only for single-line record!
	inline chrlen LineLength()	const { return RecordLength() - EOLSize(); }
};

#ifdef _FILE_WRITE

// 'TxtOutFile' represents TxtFile for writing
class TxtOutFile : public TxtFile
{
#ifdef _ISCHIP
	friend class DataOutFile;
#endif
public:
	static bool	Zipped;				// true if filed should be zippped

private:
	char	_delim;			// field delimiter in line
	// === line write buffer
	char*	_lineBuff;		// line write buffer; for writing mode only
	rowlen	_lineBuffLen;	// length of line write buffer in writing mode, otherwise 0
	rowlen	_lineBuffOffset;// current shift from the _buffLine; replaced by #define!!!
#ifdef _MULTITHREAD
	// === total counter of writed records
	ULONG*	_totalRecCnt;	// pointer to total counter of writed records; for clone only
	Mutex::eType _mtype;
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
	rowlen LineAddChars(const char* src, rowlen num, bool addDelim);

protected:
	// Returns current buffer write position
	inline rowlen CurrBuffPos() const { return _lineBuffOffset; }

public:
	// Constructs an TxtOutFile instance: allocates I/O buffer, opens an assigned file.
	//	@ftype: type of file
	//	@fName: valid file name without extention
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	//	@printName: true if file name should be printed in the exception's message
	inline TxtOutFile(FT::eType ftype, const string& fName, char delim = TAB,
		bool abortInvalid = true, bool printName = true) :
		_lineBuff(NULL), _lineBuffLen(0), _lineBuffOffset(0), _delim(delim), _mtype(FT::MutexType(ftype)),
		TxtFile(fName + FT::Ext(ftype, Zipped), eAction::WRITE, abortInvalid, printName)
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

	// Sets current position of the line write buffer.
	inline void LineSetOffset(rowlen offset)	{ _lineBuffOffset = offset; }

	// Increases current position by the specified len
	inline void LineIncrOffset(rowlen len) { _lineBuffOffset += len; }

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
		return LineAddChars(str.c_str(), rowlen(str.length()), addDelim);
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
	void LineAddInts(ULONG val1, ULONG val2, ULONG val3, bool addDelim= true);

	// Adds floating point value to the current position of the line write buffer,
	//	adds delimiter after value and increases current position.
	//	@val: value to be set
	//	@ndigit: number of digits to generate
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddFloat(float val, BYTE ndigit, bool addDelim=true);

	// Adds line to the IO buffer (from 0 to the current position).
	//	@num: number of bytes from the beginning of line buffer to be added to file buffer
	//	@offset: start position for the next writing cycle
	//	@closeLine: if true then close line by LF
	void LineToIOBuff(rowlen offset=0, bool closeLine=true);

	// Adds record to the IO buffer.
	// Generates exception if writing is fall.
	//	@src: record
	//	@len: length of the record
	//	@closeLine: if true then close line by LF
	void RecordToIOBuff	(const char *src, UINT len, bool closeLine=true);

	// Adds string to IO buffer.
	inline void StrToIOBuff	(const string& str) { RecordToIOBuff(str.c_str(), str.length()); }

	// Adds commented string to the IO buffer.
	//	@str: recorded string
	inline void CommLineToIOBuff(const string& str) { StrToIOBuff("# " + str);	}

public:
	// Allocates line write buffer.
	//	@len: length of buffer
	void SetLineBuff(rowlen len);

	// Adds to line int and float values and adds line to the IO buff.
	//	@ndigit: number if float digits to write
	//void WriteLine(int val1, float val2, BYTE ndigit);

	// Adds to line int and float values and adds line to the IO buff.
	//	@ndigit: number if float digits to write
	//void WriteLine(int val1, float val2, float val3, BYTE ndigit, const char*str);

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

// 'TabFile' represents TxtInFile consisting of lines with tab-separated filds
class TabFile : public TxtInFile
/*
 * Only number of fields defined in constructor is processed.
 * Number of fields can be equal 1, in that case ordinary plain text can be processed.
 */
{
	FT::eType	_fType;
	USHORT	_lineSpecLen;	// length of line specifier; for reading only
	short* _fieldPos;		// array of start positions in _currLine for each field
	char*	_currLine;		// current readed line; for reading only
	ULONG _estLineCnt = vUNDEF;	// estimated number of items

	// Checks if filed valid and throws exception if not.
	//	@fInd: field index
	//	return: true if field is valid
	bool IsFieldValid	(BYTE fInd) const;

	// Initializes new instance.
	//	@type: file bioinfo type
	void Init(FT::eType type);

	// Frees allocated memory
	void Release() { if (_fieldPos) delete[] _fieldPos; _fieldPos = NULL; }

protected:
	// Initializes instance by a new type, correct estimated number of lines if it's predefined
	void ResetType(FT::eType type);

public:
	// Creates new instance for reading
	//	@fName: name of file
	//	@type: file bioinfo type
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	//	@rintName: true if file name should be printed in exception's message
	TabFile(
		const string& fName,
		FT::eType type = FT::eType::UNDEF,
		eAction	mode = eAction::READ,
		bool abortInvalid = true,
		bool printName = true
	) : _fType(type), _lineSpecLen(0), _fieldPos(NULL), _currLine(NULL),
		TxtInFile(fName, mode, 1, abortInvalid, printName)
	{	if (mode != eAction::WRITE && IsGood())	Init(_fType);	}

#if defined _ISCHIP && defined  _MULTITHREAD
	// Creates a clone of TabFile class.
	// Clone is a copy of opened file with its own buffer for writing only. Used for multithreading file recording.
	//  @file: opened file which is copied
	//  @threadNumb: number of thread
	inline TabFile(const TabFile& file) :
		_fType(file._fType),
		_lineSpecLen(file._lineSpecLen),
		_fieldPos(file._fieldPos),
		_currLine(NULL),
		TxtInFile(file/*, threadNumb*/) {}
#endif

	inline ~TabFile()	{ Release(); }

	// Gets file bioinfo type
	inline FT::eType Type() const { return _fType; }

	// Gets count of readed lines
	inline ULONG Count() const { return RecordCount(); }

	// Returns estimated number of items
	inline ULONG EstCount() const { return _estLineCnt; }

	// Skip commented lines and returns estimated number of uncommented lines
	//ULONG GetUncommLineCount();

	// Reads first line and set it as current.
	// Throws exception if invalid
	//	@cntLines: returned estimated count of lines.
	//	It works properly only if lines are sorted by ascending, f.e. in sorted bed-files.
	//	return: current line
	//const char*	GetFirstLine(ULONG& cntLines);

	// Reads next line and set it as current.
	//	@checkTabs: it true then check the number of incoming fields (tabs)
	//	return: current line.
	const char* GetNextLine(bool checkTabs = true);

	// Reads string by field's index from current line without check up.
	//inline char* StrField(BYTE fInd) { return _currLine + _fieldPos[fInd]; }

	// Reads string by field's index from current line without check up.
	inline const char* StrField	(BYTE fInd)	const {	return _currLine + _fieldPos[fInd]; }

	// Reads string by field's index from current line with check up.
	const char* StrFieldValid	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? StrField(fInd) : NULL;
	}
	
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
	inline long LongField	(BYTE fInd)	const {	return atol(StrField(fInd)); }

	// Reads long by field's index from current line with check up.
	long LongFieldValid	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? LongField(fInd) : vUNDEF;
	}

	// Returns the read pointer to the beginning of the last read line and decreases line counter
	inline void RollBackLine() { RollBackLastRecord(TAB); }
};

// 'DataInFile' is a common program interface (PI) of bed/bam input files
class DataInFile
{
protected:
	int	_cID = vUNDEF;			// current readed chrom ID; int for BAM PI compatibility

public:
	// Returns estimated number of items
	virtual ULONG EstItemCount() const = 0;
	
	// Sets the next chromosome as the current one if they are different
	//	@shift: shift constant chrom mark position to the right
	//	@return: true, if new chromosome is set as current one
	virtual bool GetNextChrom(BYTE shift = 0) = 0;

	// Returns current chrom
	inline chrid GetChrom() const { return _cID; }

	// Retrieves next item's record
	virtual bool GetNextItem() = 0;
	
	// Returns current item start position
	virtual chrlen ItemStart() const = 0;
	
	// Returns current item end position
	virtual chrlen ItemEnd() const = 0;
	
	// Returns current item end position
	inline readlen ItemLength()	const { return readlen(ItemEnd() - ItemStart()); }

	// Returns true if alignment part of paired-end read
	virtual bool ItemIsPaired() const  = 0;
	
	// Returns current item value
	virtual float ItemValue() const  = 0;

	// Returns current item name
	virtual const char* ItemName() const = 0;
	
	// Gets string containing file name and current line number.
	//	@code: code of error occurs
	virtual const string LineNumbToStr(Err::eCode code = Err::EMPTY) const = 0;

	// Gets conditional file name: name if it's printable, otherwise empty string.
	virtual const string CondFileName() const = 0;

	// Returns true if item contains the strand sign
	//	Is invoked in the Feature constructor only.
	virtual bool IsItemHoldStrand() const = 0;

	// Returns current item strand: true - positive, false - negative
	virtual bool ItemStrand() const = 0;
};

// 'BedInFile' represents unified PI for reading bed file
class BedInFile : public DataInFile, public TabFile
{
	const BYTE StrandFieldInd = 5;	// index of strand field

	//char _cMark[Chrom::MaxMarkLength + 1];	// chrom mark buffer
	BYTE _scoreInd;		// index of 'score' filed (used for FBED and all WIGs)
	BYTE _chrMarkPos;	// chrom's mark position in line (BED, BedGraph) or definition line (wiggle_0)
	 
public:
	// Creates new instance for reading and open file
	//	@fName: name of file
	//	@type: file type; not used
	//	@scoreNumb: number of 'score' filed (for FBED and BedGraph)
	//	@abortInval: true if invalid instance should be completed by throwing exception
	//	@prName: true if file name should be printed in exception's message
	BedInFile(const string& fName, FT::eType type, BYTE scoreNumb, bool abortInval, bool prName) :
		_scoreInd(scoreNumb ? scoreNumb-1 : 4),
		_chrMarkPos(BYTE(strlen(Chrom::Abbr))),
		TabFile(fName, type, eAction::READ, abortInval, prName)
	{}

	// Gets pointer to the chrom mark in current line without check up
	//	for WIG only
	inline const char* ChromMark() const { return StrField(0) + _chrMarkPos; }

	// Returns estimated number of items
	inline ULONG EstItemCount() const { return  EstCount(); }

	// Sets the next chromosome as the current one if they are different
	//	@shift: shift constant chrom mark position to the right
	//	@return: true, if new chromosome is set as current one
	bool GetNextChrom(BYTE shift = 0);

	// Retrieves next item's record
	inline bool GetNextItem()	{ return GetNextLine(); }

	// Returns current item start position
	inline chrlen ItemStart()	const { return LongField(1); }

	// Returns current item end position
	inline chrlen ItemEnd()		const { return LongField(2); }

	// Returns true if alignment part of paired-end read
	inline bool ItemIsPaired()	const { return strrchr(ItemName(), '/'); }

	// Returns current item value
	inline float ItemValue()	const { return FloatFieldValid(_scoreInd); }

	// Returns current item name
	inline const char* ItemName() const { return StrFieldValid(3); }

	// Gets string containing file name and current line number.
	//	@code: code of error occurs
	inline const string LineNumbToStr(Err::eCode code = Err::EMPTY) const {
		return TxtInFile::LineNumbToStr(code);
	}

	// Gets conditional file name: name if it's printable, otherwise empty string.
	inline const string CondFileName() const { return TxtFile::CondFileName(); }

	// Returns true if item contains the strand sign
	//	Is invoked in the Feature constructor only.
	bool IsItemHoldStrand() const;

	// Returns current item strand: true - positive, false - negative
	inline bool ItemStrand() const	{ return *StrFieldValid(StrandFieldInd) == PLUS; }

#ifdef _WIG
	// Reset WIG type, score index, chrom mark position offset and estimated number of lines
	void ResetWigType(FT::eType type, BYTE scoreInd, BYTE cMarkPosOffset);

	// Inserts '0' after chrom in current line and returns point to the next decl parameter if exists
	//const char* SplitLineOnChrom();
#endif //_WIG
};

#ifdef _BAM
#ifdef _BIOSTAT		// defined in bioStat makefile
#include "../bam/BamReader.h"	// path in bioStat package
#else
#include "bam/BamReader.h"
#endif	// _BIOSTAT
#include <algorithm>
using namespace BamTools;

// 'BamInFile' represents unified PI for reading bam file
class BamInFile : public DataInFile
{
	// BamTools: http://pezmaster31.github.io/bamtools/struct_bam_tools_1_1_bam_alignment.html

	bool	_prFName;
	BamReader		_reader;
	BamAlignment	_read;
	mutable string	_rName;
	ULONG _estItemCnt = vUNDEF;	// estimated number of items

public:
	// Creates new instance for reading and open file
	//	@fName: name of file
	//	@prName: true if file name should be printed in exception's message
	BamInFile(const string& fName, bool prName);

	// Returns estimated number of items
	inline ULONG EstItemCount() const { return _estItemCnt; }

	// returns chroms count
	inline chrid ChromCount() const { return _reader.GetReferenceCount(); }

	// Sets the next chromosome as the current one if they are different
	//	@shift: shift constant chrom mark position to the right
	//	@return: true, if new chromosome is set as current one
	bool GetNextChrom(BYTE shift = 0);

	// Retrieves next item's record
	bool GetNextItem()	{ 
		return _reader.GetNextAlignmentCore(_read)
			&& _read.Position >= 0;	// additional check because of bag: GetNextAlignment() doesn't
									// return false after last read while reading the whole genome
	}

	// Returns current item start position
	inline chrlen ItemStart()	const { return _read.Position; }

	// Returns current item end position
	inline chrlen ItemEnd()		const { return _read.Position + _read.Length; }

	// Returns true if alignment part of paired-end read
	inline bool ItemIsPaired()	const { return _read.IsPaired(); }

	// Returns current item value
	inline float ItemValue()	const { return _read.MapQuality; }

	// Returns current item name
	inline const char* ItemName() const { return (_rName = _read.Name).c_str(); }

	// Gets string containing file name and current line number.
	//	@code: code of error occurs
	inline const string LineNumbToStr(Err::eCode code = Err::EMPTY) const {	return strEmpty; }

	// Gets conditional file name: name if it's printable, otherwise empty string.
	inline const string CondFileName() const { return _prFName ? _reader.GetFilename() : strEmpty; }

	// DataInFile method empty implementation.
	inline bool IsItemHoldStrand() const { return true; }

	// Returns current item strand: true - positive, false - negative
	inline bool  ItemStrand() const	{ return !_read.IsReverseStrand(); }

	// Returns SAM header data
	inline const string GetHeaderText() const { return _reader.GetHeaderText(); }
};
#endif	// _BAM

#ifndef _WIGREG

// 'Region' represents a simple region within nucleotides array (chromosome).
struct Region
{
	chrlen Start;	// start position of the region in standard chromosomal coordinates
	chrlen End;		// end position of the region in standard chromosomal coordinates

	inline Region(chrlen start=0, chrlen end=0) : Start(start), End(end) {}

	// Gets length of region.
	// The End is not included in the bases https://genome.ucsc.edu/FAQ/FAQformat.html#format1
	inline chrlen Length()	const {	return End - Start; }

	inline bool Empty()		const { return !End; }

	inline chrlen Centre()	const {	return Start + (Length()>>1); }

	// Initializes instance
	inline void Set(chrlen start, chrlen end) { Start = start; End = end; }

	inline bool operator==(const Region& r) const { return End == r.End && Start == r.Start; }

	// Returns true if this instance is invalid
	inline bool Invalid() const { return Start >= End; }

	// Returns true if Region r is covered by this instance.
	//	@r: tested Region; should be sorted by start position
	inline bool BaseCover(const Region& r) const { return r.End <= End && r.Start >= Start; }

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
	inline void Print() const { cout << Start << TAB << End << LF; }
#endif
};

// 'Regions' represents a container of defined regions within chromosome
// Defined in TxtFile.h since it's used in Fa class
class Regions
{
protected:
	vector<Region> _regions;

public:
	typedef vector<Region>::const_iterator Iter;

	// Iterator to region Begin position
	inline const Iter Begin()	const { return _regions.begin(); }

	// Iterator to region End position
	inline const Iter End()		const { return _regions.end(); }

	// Default (empty) constructor to form Chroms collection
	inline Regions() {}
	
	// Single region constructor
	inline Regions(chrlen start, chrlen end) { _regions.emplace_back(start, end); }

	// Copying constructor
	//inline Regions(const Regions& rgns) { _regions = rgns._regions;	}
	
	// Gets total length of regions.
	//chrlen Length() const;
	
	// Gets count of regions.
	inline chrlen Count()		const { return chrlen(_regions.size()); }
	
	// Gets first start position.
	inline chrlen FirstStart()	const { return _regions.front().Start; }
	
	// Gets last end position.
	inline chrlen LastEnd()		const { return _regions.back().End; }

	// Gets conditionally defined length: distance between first start and last end
	inline chrlen DefLength()	const { return LastEnd() - FirstStart(); }

	//Regions& operator=(const Regions& rgn);

	inline const Region& operator[](chrlen ind) const { return _regions[ind]; }

	// Reserves container's capacity.
	//	@count: reserved number of regions. The real number may be differ.
	inline void Reserve(chrlen count) { _regions.reserve(count); }
	
	// Clears container.
	inline void Clear()	{ _regions.clear(); }

	inline void Add	(const Region& rgn)		{ _regions.push_back(rgn); }

	// Copies subregions
	inline void Copy(const vector<Region>& source, chrlen start, chrlen stop) {
		_regions = vector<Region>(source.begin() + start, source.begin() + stop + 1);
	}

#if defined _READDENS || defined _BIOCC

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

	// Initializes this instance by external Regions.
	inline void Copy(const Regions& rgns)	{ _regions = rgns._regions; }
	
	inline void Add(chrlen start, chrlen end){ Add(Region(start, end)); }

#endif	// _READDENS, _BIOCC
#ifdef DEBUG
	void Print() const;
#endif
};

// 'ChromDefRegions' represents a container of defined regions within chromosome
// Defined in TxtFile.h since it's used in Fa class
class ChromDefRegions : public Regions
{
	// 'Combiner' consistently combines regions separated by a space less than the minimum allowed
	class Combiner
	{
		chrlen	_gapLen;	// minimum allowable gap between regions
		Region	_rgn;		// current region

	public:
		// Constructor accepting the minimum allowable gap
		inline Combiner(chrlen minGapLen) : _gapLen(minGapLen) {}

		// Returns false if the next region is separated from the current one
		// by a value less than the established minimum gap.
		// Otherwise returns true and combined region as parameter.
		bool ExceptRegion(Region& rgn);

		// Returns last excepted region
		inline const Region& LastRegion() const { return _rgn; }
	};

	static const int DefCapacuty = 12;	// default container capacity

	string	_fName;
	chrlen	_gapLen;	// total length of gaps
	mutable bool _new;

public:
	static const string Ext;	// regions file extension

	// Creates new instance and initializes it from file if one exist
	//	@cfName: chrom file name without extension
	//	@minGapLen: length, gaps less than which are ignored when reading; if 0 then read all regions 
	ChromDefRegions(const string& cfName, chrlen minGapLen=2);
	
	// Returns true if instance is empty
	inline bool Empty() const { return _new; }

	// Returns total length of gaps
	inline chrlen GapLen() const { return _gapLen; }

	// Increments total length of gaps
	inline void IncrGapLen(chrlen gapLen) { _gapLen += gapLen; }

	// Adds def region beginning of current gap position.
	//	@rgn: def region
	//	@minGapLen: minimal length which defines gap as a real gap
	void AddRegion(const Region& rgn, chrlen minGapLen);

	// Saves new instance
	void Write() const;

#if defined _READDENS || defined _BIOCC
	// Combines regions with a gap less than a specified value
	void Combine(chrlen minGapLen);
#endif
};

// 'FaFile' supports reading/writing chromosomes in FA format
class FaFile : public TxtInFile
{
	// 'DefRgnMaker' produced chrom defined regions while FaFile reading 
	class DefRgnMaker
	{
		chrlen	_minGapLen,			// minimal length which defines gap as a real gap
				_currPos;			// current readed nt position
		Region	_defRgn;			// current defined region
		ChromDefRegions& _defRgns;	// external chrom's defined regions

	public:
		inline DefRgnMaker(ChromDefRegions& rgns, chrlen minGapLen)
			: _minGapLen(minGapLen), _currPos(0), _defRgns(rgns) { _defRgns.Clear(); }

		// Adds last readed line length
		inline void AddLineLen(chrlen lineLen) { _currPos += lineLen; }

		// Remove last readed line length
		inline void RemoveLineLen(chrlen lineLen) { _currPos -= lineLen; }

		// Adds a gap to assign defined regions
		//	@start: gap's start position.
		//	@len: gap's length
		void AddGap(chrlen start, chrlen len);

		// Closes adding gaps, saves hrom's defined regions
		//	@cLen: chrom length
		void CloseAddGaps(chrlen cLen);
	};

	typedef const char* (FaFile::*tpGetLine)();
	static const char FaComment = '>';

	chrlen		_cLen;			// length of chromosome
	tpGetLine	_pGetLine;		// pointer to the 'GetNextLine' method
	DefRgnMaker	*_rgnMaker;

	// Search 'N' subsequence in the current line beginning from the start position
	//	@startPos: starting line reading position
	//	@NCnt: number of 'N' in the line beginning from the start position
	//	return: true if line trimming is complete, false for single 'N'
	void CountN(chrlen startPos, chrlen NCnt);

	// Reads line and set it as current with def filling regions
	// Since method scans the line, it takes no overhead expenses to do this.
	const char* GetLineWitnNControl();

public:
	// Opens existing file and reads first line.
	//	@rgns: def regions to fill, otherwise NULL to reading without 'N' control
	FaFile(const string& fName, ChromDefRegions* rgns=NULL);

	inline ~FaFile() { if(_rgnMaker) delete _rgnMaker; }

	// Gets chromosome's length
	inline chrlen ChromLength() const { return _cLen; }

	// Reads line and set it as current (with filling def regions if necessary)
	inline const char* NextGetLine()	{ return (this->*_pGetLine)(); }

	// Closes reading
	inline void CLoseReading()	{ if(_rgnMaker)	_rgnMaker->CloseAddGaps(_cLen); }
};

#endif	// _WIGREG
#endif	// no _FQSTATN
#if defined _CALLDIST || defined _FQSTATN

// 'FqFile' implements reading/writing file in FQ format.
class FqFile : public TxtInFile
{
	enum eLineLenIndex { HEADER1, READ, HEADER2, QUAL };

public:	
	// Creates new instance for reading by file name
	inline FqFile(const string& fileName)
		: TxtInFile(fileName, eAction::READ, 4, true, false) {}

	// Returns checked length of current readed Read.
	readlen ReadLength() const;
	
	// Gets checked Read from readed Sequence.
	const char* GetCurrRead() const;

	// Returns checked Sequence.
	const char*	GetSequence	();
	
	// Returns count of sequences.
	inline ULONG Count() const { return RecordCount(); }
};

#endif	// _CALLDIST || _FQSTATN

#if !defined _WIGREG && !defined _FQSTATN

// 'Read' represents Read (with name and score in case of _VALIGN) as item
class Read
{
public:
	static const readlen VarMinLen = 20;	// minimum Read length in variable Read mode
	static const readlen VarMaxLen = 3000;	// maximum Read length in variable Read mode
	//static const readlen sVarMaxLenLen;		// length of string representation of VarMaxLen
	static readlen	FixedLen;				// fixed length of Read
private:
#if defined _ISCHIP || defined _VALIGN
	static const char	Strands[2];			// strand markers: [0] - positive, [1] - negative
public:
	static const char	NmDelimiter = ':';		// delimiter between progTitle and chrom
	static const char	NmNumbDelimiter = DOT;	// delimiter between prog title and number
	static const char	NmPos1Delimiter = ':';	// delimiter before first recorded position
	static const char	NmPos2Delimiter = '-';	// delimiter between two recorded positions in pair
#endif

#ifdef _ISCHIP

private:
	static char		SeqQuality;			// the quality values for the sequence (ASCII)
	static readlen	LimitN;				// maximal permitted number of 'N' in Read or vUNDEF if all
	static bool		PosInName;			// true if Read name includes a position
	static const char Complements[];	// template for complementing Read

	typedef void (Read::* pCopyRead)(char*) const;

	static pCopyRead CopyRead[2];
	//static void (*spCopyRead[2])(const Read*, char*);

	const char* _seq;
	Region	_rgn;

	// Copies complemented Read into dst
	void CopyComplement(char* dst) const;

public:
	static const char* Title;
	static const char* title;

	// Initializes static members
	//	@len: length of Read
	//	@posInName: true if Read position is included in Read name
	//	@seqQual: quality values for the sequence
	//	@limN: maximal permitted number of 'N'
	static void Init(readlen len, bool posInName, char seqQual, short limN);

	inline static char StrandMark(bool reverse) { return Strands[int(reverse)]; }

	inline static bool IsPosInName() { return PosInName; }

	// Fills external buffer by quality values for the sequence
	inline static void FillBySeqQual(char* dst, readlen rlen) { memset(dst, SeqQuality, rlen); }

	// Constructor by sequence, start position and length
	Read(const char* seq, chrlen pos, readlen len) : _seq(seq), _rgn(Region(pos, pos + len)) {}

	// Gets Read's region
	//inline const Region& Rgn() const { return _rgn; }

	// Gets Read's length
	inline readlen Length() const { return _rgn.Length();  }

	// Gets Read's start position
	inline chrlen Start() const { return _rgn.Start; }

	// Gets Read's end position
	inline chrlen End() const { return _rgn.End; }

	// Gets Read's sequence
	inline const char* Seq() const { return _seq; }

	// Copies Read into dst
	inline void Copy(char* dst) const { memcpy(dst, _seq, Length()); }

	// Copies initial or complemented Read into dst
	inline void Copy(char* dst, bool reverse) const { (this->*CopyRead[reverse])(dst); }

	// Checks Read for number of 'N'
	//	return:	1: NULL Read; 0: success; -1: N limit is exceeded
	int CheckNLimit();

	// Prints quality values for the sequence
	inline static void PrintSeqQuality() { cout << '[' << SeqQuality << ']'; }

	// Prints Read values - parameters
	//	@signOut: output marker
	//	@isRVL: true if Read variable length is set
	static void PrintParams(const char* signOut, bool isRVL);

#else
public:
	chrlen	Pos;		// Read's actual start position
	readlen Len;		// Read's length
	ULLONG	Numb;		// read number keeped in name
	bool	Strand;		// true if strand is positive

	chrlen Centre() const { return Pos + (Len >> 1); }

#ifdef _VALIGN
	chrlen	RecPos;		// recorded (true) Read start position
	float	Score;		// Read's score
	chrid	InitCID;	// initial chrom - owner

	// Cobstructs extended Read
	//	@cid: chrom ID
	//	@num: uniq number within chrom or initial position
	Read(const Region& rgn, ULLONG numb, bool strand, chrlen recPos, float score, chrid cid) :
		Pos(rgn.Start),
		Len(readlen(rgn.Length())),
		Numb(numb),
		Strand(strand),
		RecPos(recPos),
		Score(score),
		InitCID(cid)
	{}
#elif defined _CALLDIST

	Read(const Region& rgn, ULLONG numb, bool strand)
		: Pos(rgn.Start), Len(readlen(rgn.Length())), Numb(numb), Strand(strand) {}
#else	// _BIOCC
	Read(const Region& rgn) 
		: Pos(rgn.Start), Len(readlen(rgn.Length())), Numb(0L), Strand(true) {}
#endif
	// Compares two Reads by position. For sorting a container.
	inline static bool CompareByStartPos(const Read& r1, const Read& r2) { return r1.Pos < r2.Pos; }

	//inline static bool CompareByNum(const Read& r1, const Read& r2) {	return r1.Num < r2.Num; }

	void inline Print() const {
		dout << Pos << TAB << Numb
#ifdef _VALIGN
			<< TAB << int(Score)
#endif
			<< LF;
	}
#endif
};
#endif