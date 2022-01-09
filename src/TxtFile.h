/**********************************************************
TxtFile.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07.01.2022
-------------------------
Provides read|write text file functionality
***********************************************************/

#pragma once
#ifndef _TXTFILE_H
#define _TXTFILE_H

#include "common.h"
//#include <variant>

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
	const BYTE	AvrLineLen;		// average line length; to reserve container size
	const char	Comment;		// char indicates that line is comment
	const char* LineSpec;		// substring on which each data line is beginning

	inline TabFilePar() : MinFieldCnt(0), MaxFieldCnt(0), AvrLineLen(0), Comment(cNULL), LineSpec(NULL) {}

	TabFilePar(BYTE minTabCnt, BYTE maxTabCnt, BYTE avrLineLen = 0, char comm = HASH, const char* lSpec = NULL) :
		MinFieldCnt(minTabCnt),
		MaxFieldCnt(maxTabCnt<minTabCnt ? minTabCnt : maxTabCnt),
		AvrLineLen(avrLineLen),
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
 * and less than size_t.
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
		ISCR		= 0x01,	// sign of Carriage Return ('\r')
		LFCHECKED	= 0x02,	// the presence of a symbol CR is checked; for Reading mode
		ZIPPED		= 0x04,	// file is zipped
		ABORTING	= 0x08,	// invalid file should be completed by throwing exception; for Reading mode
		ENDREAD		= 0x10,	// last call of GetNextRecord() has returned NULL; for Reading mode
		PRNAME		= 0x20,	// print file name in exception's message; for Reading mode
		MTHREAD		= 0x40,	// file in multithread mode: needs to be locked while writing
		CLONE		= 0x80,	// file is a clone
	};

private:
	const static int BlockSize = 1024 * 1024;	// 1 Mb
	const int CRLF = eFlag::ISCR + eFlag::LFCHECKED;	// combined 'flag'

	LLONG	_fSize;			// the length of uncompressed file; for zipped file more than
							// 4'294'967'295 its unzipped length is unpredictable
	string	_fName;			// file's name
	mutable short _flag;	// bitwise storage for signs included in eFlag

protected:
	void *	_stream;			// FILE* (for unzipped file) or gzFile (for zipped file)
	char *	_buff;				// basic I/O (read/write) buffer
	size_t	_buffLen;			// the length of the basic I/O buffer
	mutable size_t _currRecPos;	// start position of the last readed/writed record in the block
	mutable ULONG _recCnt;		// local counter of readed/writed records
	mutable Err::eCode	_errCode;
	//Stopwatch	_stopwatch;
	
	inline void SetFlag	(eFlag f, bool val)	const { val ? _flag |= f : _flag &= ~f; }
	inline void RaiseFlag(eFlag f)	const { _flag |= f; }
	inline bool IsFlag	(eFlag f)	const { return _flag & f; }
	inline bool IsZipped()			const { return _flag & ZIPPED; }
	inline bool IsClone()			const { return _flag & CLONE; }

	// *** 3 methods used by TxtInFile only

	// Gets the number of characters corresponded to LF
	// In Windows LF matches '\r\n', in Linux LF matches '\n',
	//	return: in Windows always 1, in Linux: 2 for file created in Windows, 1 for file created in Linux
	inline size_t LFSize() const	{ 	return 2 - (_flag & ISCR); }

	// Returns true if LF size is not defined
	inline bool IsLFundef() const { return !(_flag & LFCHECKED); }

	// Establishes the presence of CR symbol at the end of line.
	//	@c: if c is CR then the second bit is raised to 1, the first turn down to 0,
	//	so the value return by LFSZ mask is 2, otherwise remains in state 1
	//inline void SetLF(char c)	{ if(c==CR)	 RaiseFlag(ISCR); RaiseFlag(LFCHECKED); }
	inline void SetLF(char c) { if (c == CR)	_flag |= CRLF; }

private:
	// Initializes instance variables, opens a file, sets a proper error code.
	//	@fName: valid full name of file
	//	@mode: opening mode
	//	@fStream: clonable file stream or NULL
	//	return: true is success, otherwise false.
	bool SetBasic(const string& fName, eAction mode, void* fStream);
	
	// Allocates memory for the I/O buffer with checking.
	//	return: true if successful
	bool CreateIOBuff();

	//void operator = (const TxtFile&); //assignment prohibition

protected:
	// Constructs an TxtFile instance: allocates I/O buffer, opens an assigned file.
	//	@fName: valid full name of assigned file
	//	@mode: opening mode
	//	@msgFName: true if file name should be printed in the exception's message
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	TxtFile(const string& fName, eAction mode, bool msgFName =true, bool abortInvalid = true);

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

	// Returns true if instance is valid.
	inline bool IsGood() const	{ return _errCode == Err::NONE; }
	inline bool IsBad()  const	{ return _errCode != Err::NONE; }

	// Returns current error code.
	//inline Err::eCode ErrCode() const	{ return _errCode; }

//public:
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
	inline void ThrowExcept(const string& msg) const { Err(msg, CondFileName()).Throw(); }

	// Throws exception
//	@msg: exception message
	inline void ThrowExcept(Err::eCode code) const { Err(code, CondFileName().c_str()).Throw(); }
};

// 'TxtInFile' represents TxtFile for reading.
//	Ñorrectly reads records even if the latter does not end with LF (CR,LF) character (s)
class TxtInFile : public TxtFile
{
	size_t	_recLen;		// the length of record with LF marker
	size_t	_readedLen;		// number of actually readed chars in block
	size_t*	_linesLen;		// array of lengths of lines in a record
	BYTE	_recLineCnt;	// number of lines in a record

	// Raises ENDREAD sign an return NULL
	//inline char* ReadingEnded()		  { RaiseFlag(ENDREAD); return NULL; }

	// Returns record without control
	inline char* RealRecord() const { return _buff + _currRecPos - _recLen; }

	// Reads next block.
	//	@offset: shift of start reading position
	//	return: 0 if file is finished; -1 if unsuccess reading; otherwhise number of readed chars
	int ReadBlock(const size_t offset);

	// Returns true if block is complete, move remainder to the top
	bool CompleteBlock(size_t currLinePos, size_t blankLineCnt);

	// Fills I/O buffer with 0, beginning from @offset position
	//inline void ClearBuff(size_t offset = 0) { memset(_buff + offset, 0, _buffLen - offset); }

	// Gets string containing file name and current record number.
	//inline const string RecordNumbToStr() const { return LineNumbToStr(0); }

	// Gets number of line.
	//	@lineInd: index of line in a record
	ULONG LineNumber(BYTE lineInd) const { return (_recCnt-1) * _recLineCnt + lineInd + 1; }

protected:
	// Constructs an TxtInFile instance: allocates buffers, opens an assigned file.
	//	@fName: valid full name of assigned file
	//	@mode: opening mode: READ or READ_ANY
	//	@cntRecLines: number of lines in a record
	//	@msgFName: true if file name should be printed in an exception's message
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	TxtInFile(const string& fName, eAction mode, BYTE cntRecLines,
		bool msgFName = true, bool abortInvalid=true);

	inline ~TxtInFile() { if(_linesLen)	delete [] _linesLen; }

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
	//	@withoutLF: if true then without LF marker length, otherwise without it
	inline size_t	LineLengthByInd	(BYTE lineInd, bool withoutLF = true) const { 
		return _linesLen[lineInd] - withoutLF * LFSize();
	}

	// Throw exception if no record is readed.
	void CheckGettingRecord() const { 
		if (IsLFundef())
			ThrowExcept("attempt to get record's oinfo without reading record");
	}
#endif

	// Returns the read pointer to the beginning of the last read line and decreases line counter.
	//	Zeroes length of current reading record!
	//	@sep: field separator character
	void RollBackRecord(char sep);

public:
	// Gets length of current reading record including LF marker.
	//	Returns 0 after RollBackLastRecord() invoke.
	inline size_t	RecordLength() const { return _recLen; }

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

	// Throws exception with message included current reading line number
	//	@msg: exception message
	void ThrowExceptWithLineNumb(const string& msg) const { Err(msg, LineNumbToStr().c_str()).Throw(); }

	// Gets length of current line without LF marker: only for single-line record!
	inline chrlen LineLength()	const { return RecordLength() - LFSize(); }
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
	//typedef void(TxtOutFile::* fAddChar)();
	//// 'Add delimiter' methods: [0] - empty method, [1] - adds delimiter and increases current position
	//static fAddChar fLineAddChar[2];

	char	_delim;				// field delimiter in line
	// === line write buffer
	char*	_lineBuff = NULL;		// line write buffer; for writing mode only
	rowlen	_lineBuffLen = 0;		// length of line write buffer in writing mode, otherwise 0
	rowlen	_lineBuffOffset = 0;	// current shift from the _buffLine; replaced by #define!!!
#ifdef _MULTITHREAD
	// === total counter of writed records
	ULONG*	_totalRecCnt;	// pointer to total counter of writed records; for clone only
	Mutex::eType _mtype;
#endif

	//inline void AddCharEmpty() {}
	////	Adds delimiter to the current position in the line write buffer and increases current position.
	//inline void AddDelim() { LineAddChar(_delim); }

	//	Adds delimiter to the current position in the line write buffer and increases current position.
	//	@add: if true then adds delimiter
	//inline void LineAddDelim(bool add) { (this->*fLineAddChar[add])(); }
	inline void LineAddDelim(bool add) { if(add) LineAddChar(_delim); }

	// Allocates memory for write line buffer with checking.
	//	return: true if successful
	bool CreateLineBuff(rowlen len);

	// Closes adding record to the IO buffer: set current rec position and increases rec counter
	//	@len: length of added record
	void EndRecordToIOBuff(size_t len);

protected:
	// Returns current buffer write position
	inline rowlen CurrBuffPos() const { return _lineBuffOffset; }

	// Adds character to the current position in the line write buffer.
	//	@ch: char to be set
	inline void LineAddChar(char ch) { _lineBuff[_lineBuffOffset++] = ch; }

	// Adds character to the current position in the line write buffer with optional adding delimiter
	//	@ch: char to be set
	//	@addDelim: if true then adds delimiter and increases current position
	void LineAddChar(char ch, bool addDelim);

	// Copies block of chars to the current position in the line write buffer,
	//	adds delimiter after string	and increases current position.
	//	@src: pointer to the block of chars
	//	@num: number of chars
	//	@addDelim: if true then adds delimiter and increases current position
	//	return: new current position
	rowlen LineAddChars(const char* src, rowlen num, bool addDelim = true);

//public:
	// Constructs an TxtOutFile instance: allocates I/O buffer, opens an assigned file.
	//	@ftype: type of file
	//	@fName: valid file name without extention
	//	@printName: true if file name should be printed in the exception's message
	//	@abortInvalid: true if invalid instance shold be completed by throwing exception
	inline TxtOutFile(FT::eType ftype, const string& fName,
		char delim = TAB, bool printName = true, bool abortInvalid = true) :
		_delim(delim), _mtype(FT::MutexType(ftype)),
		TxtFile(fName + FT::Ext(ftype, Zipped), eAction::WRITE, printName, abortInvalid)
		{
#ifdef _MULTITHREAD
			_totalRecCnt = &_recCnt;	// for atomic increment
#endif
		}		// line buffer will be created in SetLineBuff()

	// Writes nonempty buffer, deletes line buffer and closes file
	~TxtOutFile();

//protected:
#ifdef _MULTITHREAD
	// Constructs a clone of an existing instance.
	// Clone is a copy of opened file with its own separate basic I/O and write line buffers.
	// Used for multithreading file recording
	//	@file: opened file which is cloned
	TxtOutFile	(const TxtOutFile& file);
#endif

	// Sets current position of the line write buffer.
	inline void LineSetOffset(rowlen offset = 0)	{ _lineBuffOffset = offset; }

	// Increases current position by the specified len
	inline void LineIncrOffset(rowlen len) { _lineBuffOffset += len; }

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
	void LineAddInt(LLONG val, bool addDelim = true);

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
	//void LineAddFloat(float val, BYTE ndigit, bool addDelim=true);

	// Adds line to the IO buffer (from 0 to the current position).
	//	@num: number of bytes from the beginning of line buffer to be added to file buffer
	//	@offset: start position for the next writing cycle
	void LineToIOBuff(rowlen offset=0);

	// Adds record to the IO buffer.
	// Generates exception if writing is fall.
	//	@src: record
	//	@len: length of the record
	void RecordToIOBuff	(const char *src, size_t len);

	// Adds string to IO buffer without checking buffer exceeding.
	void StrToIOBuff(const string&& str);

	// Adds commented string to the IO buffer.
	//	@str: recorded string
	inline void CommLineToIOBuff(const string& str) { StrToIOBuff("# " + str);	}

	// Allocates line write buffer.
	//	@len: length of buffer
	void SetLineBuff(rowlen len);
	
	// Adds to line 2 int values and adds line to the IO buff.
	//void WriteLine(ULONG val1, ULONG val2);

	// Adds to line C-string with delimiter, and adds line to the IO buff.
	//void WriteLine(const char* str);

	// Adds to line string with delimiter and int value, and adds line to the IO buff.
	//void WriteLine(const string& str, int val);

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
 * Number of fields can be equal 1, in that case ordinary plain text is processed.
 */
{
	FT::eType	_fType;
	USHORT	_lineSpecLen = 0;		// length of line specifier; for reading only
	short*	_fieldPos = nullptr;	// array of start positions in _currLine for each field
	char*	_currLine = nullptr;	// current readed line; for reading only
	ULONG	_estLineCnt = vUNDEF;	// estimated number of items

	// Checks if filed valid and throws exception if not.
	//	@fInd: field index
	//	return: true if field is valid
	bool IsFieldValid	(BYTE fInd) const;

	// Initializes new instance.
	//	@type: file bioinfo type
	//	@estLineCnt: if true then estimate count of lines
	void Init(FT::eType type, bool estLineCnt);

	// Frees allocated memory
	void Release() { if (_fieldPos) delete[] _fieldPos; _fieldPos = NULL; }

	// Returns the read pointer to the beginning of the last read line and decreases line counter
	inline void RollBackLine() { RollBackRecord(TAB); }

protected:
	// Sets estimation of number of lines by currently readed line, and rolls line back
	void SetEstLineCount();

	// Sets estimation of number of lines by type (for WIGGLE)
	void SetEstLineCount(FT::eType type) {
		_estLineCnt = ULONG(Length() / FT::FileParams(type).AvrLineLen);
	}

	// Initializes instance by a new type, correct estimated number of lines if it's predefined
	void ResetType(FT::eType type);

public:
	// Creates new instance for reading
	//	@fName: name of file
	//	@type: file bioinfo type
	//	@msgFName: true if file name should be printed in exception's message
	//	@abortInvalid: true if invalid instance should be completed by throwing exception
	TabFile(
		const string& fName,
		FT::eType type = FT::eType::UNDEF,
		eAction	mode = eAction::READ,
		bool estLineCnt = true,
		bool msgFName = true,
		bool abortInvalid = true
		) : TxtInFile(fName, mode, 1, msgFName, abortInvalid), _fType(type)
	{	
		if (mode != eAction::WRITE && IsGood())	
			Init(_fType, estLineCnt);	}

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

	// Returns a pointer to the substring defined by key.
	//	@str: null-terminated string to search the key
	//	@key: string to search for
	//	return: a pointer to the substring after key, or NULL if key does not appear in str
	static const char* KeyStr(const char* str, const string& key);

	// Checks definition or declaration line for key
	//	@str: null-terminated string to search the key
	//	@key: string to search for
	//	return: point to substring followed after the key
	const char* CheckSpec(const char* str, const string& key);

	// Returns required int value with check
//	@str: null-terminated string to search the key
//	@key: string to search for
//	return: key value, or throws an exception if key does not appear in str
	inline chrlen GetIntKey(const char* str, const string& key) { return atoi(CheckSpec(str, key) + 1); }

	// Gets file bioinfo type
	inline FT::eType Type() const { return _fType; }

	// Gets count of readed lines
	inline size_t Count() const { return RecordCount(); }

	// Returns estimated number of items
	inline size_t EstLineCount() const { return _estLineCnt; }

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

	// Gets current line.
	inline char* GetLine() const { return _currLine; }

	// Reads null-terminated string by field's index from current line without check up.
	inline const char* StrField	(BYTE fInd)	const {	return _currLine + _fieldPos[fInd]; }

	// Reads null-terminated string by field's index from current line with check up.
	const char* StrFieldValid	(BYTE fInd)	const {
		return IsFieldValid(fInd) ? StrField(fInd) : nullptr;
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
};


#ifndef _WIGREG

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

	const int DefCapacuty = 12;	// default container capacity

	string	_fName;
	chrlen	_gapLen;	// total length of gaps
	mutable bool _new;

public:
	static const string Ext;	// regions file extension

	// Creates new instance and initializes it from file if one exist
	//	@fName: chrom file name without extension
	//	@minGapLen: length, gaps less than which are ignored when reading; if 0 then read all regions 
	ChromDefRegions(const string& fName, chrlen minGapLen=2);
	
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

	const char FaComment = '>';
	chrlen		_cLen;						// length of chromosome
	const char* (FaFile::* _pGetLine)();	// pointer to the 'GetNextLine' method
	unique_ptr<DefRgnMaker> _rgnMaker;		// chrom defined regions store

	// Search 'N' subsequence in the current line beginning from the start position
	// and adds complete 'N' subsequence to _rgnMaker
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
		: TxtInFile(fileName, eAction::READ, 4, false) {}

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

#endif	//_TXTFILE_H