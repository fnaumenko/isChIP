/**********************************************************
TxtFile.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07.01.2022
-------------------------
Provides read|write text file functionality
***********************************************************/

#include "TxtFile.h"
#include <fstream>	// to write ChromDefRegions
#ifndef _FILE_WRITE
#include <fstream>
#endif

const BYTE TabFilePar::BGLnLen = Chrom::MaxAbbrNameLength + 2 * 9;	// 2*pos + correction
const BYTE TabFilePar::WvsLnLen = 9 + 3 + 2 + 25;	// pos + val + TAB + LF + correction
const BYTE TabFilePar::WfsLnLen = 5 + 1;			// val + LF

/************************ class FT ************************/

const char* FT::bedExt = "bed";
const char* FT::wigExt = "wig";
const string FT::Interval = "interval";
const string FT::Intervals = "intervals";
const string FT::Read = "read";
const string FT::Reads = "reads";

const FT::fTypeAttr FT::TypeAttrs[] = {
	{ "",		strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabFilePar(1, 1) },	// undefined type
	{ bedExt,	"feature",	"features",	Mutex::eType::WR_BED,TabFilePar(3, 6, 0, HASH, Chrom::Abbr) },	// ordinary bed
	{ bedExt,	Read,		Reads,		Mutex::eType::WR_BED,TabFilePar(6, 6, 0, HASH, Chrom::Abbr) },	// alignment bed
	{ "sam",	strEmpty,	strEmpty,	Mutex::eType::WR_SAM,TabFilePar(0, 0) },
	{ "bam",	Read,		Reads,		Mutex::eType::NONE,	 TabFilePar() },
	{ wigExt,	Interval,	Intervals,	Mutex::eType::NONE,	 TabFilePar(4, 4, TabFilePar::BGLnLen, HASH) },	// bedgraph: Chrom::Abbr isn't specified becuase of track definition line
	{ wigExt,	Interval,	Intervals,	Mutex::eType::NONE,	 TabFilePar(1, 1, TabFilePar::WvsLnLen, HASH) },// wiggle_0 variable step
	{ wigExt,	Interval,	Intervals,	Mutex::eType::NONE,	 TabFilePar(1, 1, TabFilePar::WfsLnLen, HASH) },// wiggle_0 fixed step
	{ "fq",		Read,		Reads,		Mutex::eType::WR_FQ, TabFilePar() },
	{ "fa",		strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabFilePar() },
	// for the chrom.size, do not specify TabFilePar::LineSpec to get an exception if the file is invalid
	{ "chrom.sizes",strEmpty,strEmpty,	Mutex::eType::NONE,	 TabFilePar(2, 2, 0, cNULL) },//, Chrom::Abbr) },
	{ "region",	strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabFilePar(2, 2) },
	{ "dist",	strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabFilePar(1, 2) },	// required 2 data fields, but set min=1 to skip optional text line
#ifdef _ISCHIP
	{ "ini",	strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabFilePar(4, 4) },	// isChIP ini file type
#endif
};
const BYTE FT::Count = sizeof(FT::TypeAttrs)/sizeof(FT::fTypeAttr);

// Returns file type
//	@fName: file name (with case insensitive extension)
//	@isABED: if true then returns ABED in case of .bed extention
FT::eType FT::GetType(const char* fName, bool isABED)
{
	const string s_ext = FS::GetExt(fName);
	if (s_ext == "fastq")	return eType::FQ;
	const char* ext = s_ext.c_str();
	eType type = eType::UNDEF;
	for(int i = int(eType::BED); i<Count; i++)	// start from ordinary bed
		if (!_stricmp(ext, TypeAttrs[i].Extens)) {
			type = eType(i); break;
		}
	if (type == eType::BED && isABED)	type = eType::ABED;
	return type;
}

// Gets file extension, beginning at DOT and adding .gz if needed
//	@t: file type
//	@isZip: true if add ".gz"
const string FT::Ext(eType t, bool isZip)
{
	string ext = string(1, DOT) + TypeAttrs[int(t)].Extens;
	if(isZip)	ext += ZipFileExt;
	return ext;
}

/************************ end of class FT ************************/

/************************ TxtFile ************************/
const char* modes[] = { "r", "w", "a+" };
const char* bmodes[] = { "rb", "wb" };

// Sets error code and throws exception if it is allowed.
void TxtFile::SetError(Err::eCode errCode, const string& senderSpec, const string& spec) const
{
	_errCode = errCode;
	if (IsFlag(ABORTING))
		Err(errCode, (CondFileName() + senderSpec).c_str(), spec.length() ? spec.c_str() : NULL).Throw();
}

// Initializes instance variables, opens a file, sets a proper error code.
//	@fName: valid full name of file
//	@mode: opening mode
//	@fStream: clonable file stream or NULL
//	return: true if success, otherwise false.
bool TxtFile::SetBasic(const string& fName, eAction mode, void* fStream)
{
	_buff = NULL;
	_stream = NULL;
	_errCode = Err::NONE;
	_fName = fName;
	_currRecPos = _recCnt = 0;
	_buffLen = BlockSize;	// by default; can be corrected
#ifdef _NO_ZLIB
	if(IsZipped()) { SetError(Err::FZ_BUILD); return false; }
#endif
	// set file stream
	if(fStream)	_stream = fStream;
	else 
#ifndef _NO_ZLIB
		if(IsZipped())
			if(mode == eAction::READ_ANY)	SetError(Err::FZ_OPEN);
			else {
				if( !(_stream = gzopen(fName.c_str(), bmodes[int(mode)])))
					SetError(Err::F_OPEN);
			}
		else
#endif
			if( _stream = fopen(fName.c_str(), modes[int(mode)]))
				setvbuf((FILE*)_stream, NULL, _IONBF, 0);
			else
				SetError(Err::F_OPEN);
	return IsGood();
}

// Allocates memory for the I/O buffer buffer with checking.
//  return: true if successful
bool TxtFile::CreateIOBuff()
{
	try { _buff = new char[_buffLen]; }
	catch(const bad_alloc)	{ SetError(Err::F_MEM); }
	return IsGood();
}

// Constructs an TxtFile instance: allocates I/O buffer, opens an assigned file.
//	@fName: valid full name of assigned file
//	@mode: opening mode
//	@msgFName: true if file name should be printed in exception's message
//	@abortInvalid: true if invalid instance shold be completed by throwing exception
TxtFile::TxtFile (const string& fName, eAction mode, bool msgFName, bool abortInvalid) :
	_flag(1)	// LF size set to 1
{
	SetFlag(ZIPPED, FS::HasGzipExt(fName));
	SetFlag(ABORTING, abortInvalid);
	SetFlag(PRNAME, msgFName);
	if( !SetBasic(fName, mode, NULL) )	return;
	_fSize = FS::Size(fName.c_str());
	if(_fSize == -1)	_fSize = 0;		// new file
#ifndef _NO_ZLIB
	else if(IsZipped()) {				// existed file
		LLONG size = FS::UncomressSize(fName.c_str());
		if( size > 0 ) {
			// wrong uncompressed size: increase initial size four times
			// since zip is too big to keep right size in archive
			if( size <= _fSize)	_fSize <<= 2;
			// right uncompressed size
			else				_fSize = size;
		}
		//_buffLen >>= 1;		// decrease block size twice because of allocating additional memory:
							// 2x_buffLen for writing or 3x_buffLen for reading by gzip
	}
#endif
	if( _fSize && _fSize < _buffLen )
		if( mode != eAction::WRITE )
			_buffLen = (ULONG)_fSize + 1;	// for reading
		else if( _fSize * 2 < _buffLen )
			_buffLen = (ULONG)_fSize * 2;	// for writing: increase small buffer for any case
	
	if( !CreateIOBuff() )	return;

#ifdef ZLIB_NEW
	if( IsZipped() && gzbuffer((gzFile)_stream, _buffLen) == -1 )
		SetError(Err::FZ_MEM);
#endif
}

#ifdef _MULTITHREAD
// Constructs a clone of an existing instance.
// Clone is a copy of opened file with its own separate write basic and write line buffers.
//	Used for multithreading file recording
//	@file: opened file which is cloned
//	Exception: file_error
TxtFile::TxtFile(const TxtFile& file) : _flag(file._flag)
{
	if( !SetBasic("Clone " + file._fName, eAction::WRITE, file._stream) )	return;
	RaiseFlag(CLONE);
	RaiseFlag(MTHREAD);
	file.RaiseFlag(MTHREAD);
	CreateIOBuff();
}
#endif

// Destructs the instance: close assigned file, releases I/O buffer.
TxtFile::~TxtFile()
{
	if(IsClone())	return;
	//_stopwatch.Stop(_fName);
	if(_buff)		delete [] _buff;
	if(_stream &&
#ifndef _NO_ZLIB
			IsZipped() ? gzclose( (gzFile)_stream) :
#endif
			fclose( (FILE*)_stream) )
		SetError(Err::F_CLOSE);
	//cout << "Free " << _fName << "\trecords = " << _recCnt << endl;
}

/************************ TxtFile: end ************************/

/************************ TxtInFile ************************/

// Constructs an TxtInFile instance: allocates buffers, opens an assigned file.
//	@fName: valid full name of assigned file
//	@mode: opening mode: READ or READ_ANY
//	@cntRecLines: number of lines in a record
//	@msgFName: true if file name should be printed in an exception's message
//	@abortInvalid: true if invalid instance should be completed by throwing exception
TxtInFile::TxtInFile(const string& fName, eAction mode, 
	BYTE cntRecLines, bool msgFName, bool abortInvalid) :
	_linesLen(NULL),
	_recLineCnt(cntRecLines),
	_recLen(0),
	_readedLen(0),
	TxtFile(fName, mode, msgFName, abortInvalid)
{
	//ClearBuff();
	if(Length() && ReadBlock(0) >= 0)		// read first block
		_linesLen = new size_t[cntRecLines];	// nonempty file: set lines buffer
	else
		RaiseFlag(ENDREAD);					// empty file
}

// Reads next block.
//	@offset: shift of start reading position
//	return: 0 if file is finished; -1 if unsuccess reading; otherwhise number of readed chars
int TxtInFile::ReadBlock(const size_t offset)
{
	size_t readLen;
	//_stopwatch.Start();
#ifndef _NO_ZLIB
	if( IsZipped() ) {
		int len = gzread((gzFile)_stream, _buff + offset, _buffLen - offset);
		if(len < 0) { SetError(Err::F_READ); return -1; }
		readLen = len;
	}
	else
#endif
	{
		readLen = fread(_buff + offset, sizeof(char), _buffLen - offset, (FILE*)_stream);
		if(readLen != _buffLen - offset && !feof((FILE*)_stream) )
		{ SetError(Err::F_READ); return -1; }
	}
	_readedLen = readLen + offset;
//#ifdef ZLIB_OLD
	//if( _readTotal + _readedLen > _fSize )
	//	_readedLen = _fSize - _readTotal;
	//_readTotal += _readedLen;
//#endif
	if( _readedLen != _buffLen && !IsZipped() && ferror((FILE*)_stream) )
	{ SetError(Err::F_READ); return -1; }
	_currRecPos = 0;
	//_stopwatch.Stop();
	return _readedLen;
}

// Returns true if block is complete, move remainder to the top
bool TxtInFile::CompleteBlock(size_t currLinePos, size_t blankLineCnt)
{
	if(_readedLen != _buffLen)	// final block, normal completion
	{ RaiseFlag(ENDREAD); return true; }
	if( !currLinePos )			// block is totally unreaded
	{	SetError(Err::F_BIGLINE); RaiseFlag(ENDREAD); return true; }

	blankLineCnt = _readedLen - _currRecPos - blankLineCnt;	// now length of unreaded remainder
	// move remainder to the beginning of buffer; if length of unreaded remain = 0, skip moving
	memmove(_buff, _buff + _currRecPos, blankLineCnt);
	//ClearBuff(blankLineCnt);
	_recLen = 0;
	if(!ReadBlock(blankLineCnt))
	{ RaiseFlag(ENDREAD); return true; };
	return false;
}

// Reads record
//	return: pointer to line or NULL if no more lines
const char* TxtInFile::GetNextRecord()
{
// Sets _currLinePos to the beginning of next non-empty line inread/write buffer
// GetNextRecord(), GetRecordN() and GetRecordTab() are quite similar,
// but they are separated for effectiveness
	if(IsFlag(ENDREAD))	return NULL;
	 
	size_t i, blanklCnt = 0,			// counter of empty lines
		 currPos = _currRecPos;		// start position of current readed record
	char* buf;

	_recLen = 0;
	for(BYTE rec=0; rec<_recLineCnt; rec++)
		for (buf = _buff + (i = currPos);;buf++, i++) {
			if(*buf == LF) {
				if( i==currPos ) {				// LF marker is first in line
					++currPos; ++blanklCnt;		// skip empty line
					continue;
				}
				if (IsLFundef())	SetLF(*(buf - 1));		// define LF size
lf:				_recLen += (_linesLen[rec] = ++i - currPos);
				currPos = i;
				break;							// mext line in a record
			}
			if(i >= _readedLen) {	// check for oversize buffer
				//cout << ">>> " << i << TAB << currPos << TAB << _readedLen << TAB << _buffLen << LF;
				if (_readedLen != _buffLen && i > currPos)	// last record does not end with LF
					goto lf;
				if(CompleteBlock(currPos, blanklCnt))	return NULL;
				currPos = blanklCnt = i = rec = 0;
				buf = _buff;
			}
		}
	_currRecPos = currPos;			// next record position
	_recCnt++;
	return RealRecord();
}

// Reads N-controlled record
//	@counterN: counter of 'N'
//	return: pointer to line or NULL if no more lines
const char* TxtInFile::GetNextRecord(chrlen& counterN)
{
	if(IsFlag(ENDREAD))	return NULL;

	size_t i, blanklCnt = 0,			// counter of empty lines
		 currPos = _currRecPos;		// start position of current readed record
	chrlen cntN = 0;				// local counter of 'N'
	char* buf;

	_recLen = 0;
	for(BYTE rec=0; rec<_recLineCnt; rec++)
		for(buf = _buff + (i=currPos);; buf++,i++)	{
			if(*buf == cN)		cntN++;
			else if(*buf == LF) {
				if( i==currPos ) {				// LF marker is first in line
					currPos++; blanklCnt++;		// skip empty line
					continue;
				}
				if (IsLFundef())	SetLF(*(buf - 1));		// define LF size
lf:				_recLen += (_linesLen[rec] = ++i - currPos);
				currPos = i;
				break;
			}
			if (i >= _readedLen) {	// check for oversize buffer
				if (_readedLen != _buffLen && i > currPos)	// last record does not end with LF
					goto lf;
				if(CompleteBlock(currPos, blanklCnt))	return NULL;
				cntN = currPos = blanklCnt = i = rec = 0;
				buf = _buff;
			}
		}
	_currRecPos = currPos;			// next record position
	_recCnt++;
	counterN += cntN;
	return RealRecord();
}

// Reads tab-controlled record
//	@tabPos: TAB's positions array that should be filled
//	@cntTabs: maximum number of TABS in TAB's positions array
//	return: pointer to line or NULL if no more lines
char* TxtInFile::GetNextRecord(short* const tabPos, const BYTE tabCnt)
{
	if(IsFlag(ENDREAD))	return NULL;
	 
	size_t i, blanklCnt = 0,			// counter of empty lines
		 currPos = _currRecPos;		// start position of current readed record
	BYTE tabInd = 1;				// index of TAB position in tabPos
	char* buf;

	_recLen = 0;
	for(BYTE rec=0; rec<_recLineCnt; rec++)
		for (buf = _buff + (i = currPos);; buf++, i++) {
			if(*buf == TAB) {
				if(tabInd < tabCnt)
					tabPos[tabInd++] = short(i + 1 - currPos);
			}
			else if(*buf == LF) {
				if( i==currPos ) {				// LF marker is first in line
					currPos++; blanklCnt++;		// skip empty line
					continue;
				}
				if (IsLFundef())	SetLF(*(buf - 1));		// define LF size
lf:				_recLen += (_linesLen[rec] = ++i - currPos);
				currPos = i;
				break;
			}
			if (i >= _readedLen) {	// check for oversize block
				if (_readedLen != _buffLen && i > currPos)	// last record does not end with LF
					goto lf;
				if(CompleteBlock(currPos, blanklCnt))	return NULL;
				currPos = blanklCnt = i = rec = 0;
				buf = _buff;
				tabInd = 1;
			}
		}
	_currRecPos = currPos;			// next record position
	_recCnt++;
	return RealRecord();
}

// Returns the read pointer to the beginning of the last read line and decreases line counter. 
//	Zeroes length of current reading record! 
//	For this reason it cannot be called more then once after each GetNextRecord() invoke
//	@sep: field separator character
void TxtInFile::RollBackRecord(char sep)
{
	for (size_t i = _recLen; i > 1 ; --i)
		if (!_buff[_currRecPos - i])
			_buff[_currRecPos - i] = sep;	// return back TAB in the last record instead of inserted '0'
	_buff[_currRecPos - 1] = LF;			// return back LineFeed in the end of last record instead of inserted '0'
	_currRecPos -= _recLen;
	_recLen = 0;
	_recCnt--;
}

// Gets string containing file name and current line number.
//	@code: code of error occurs
//	@lineInd: index of line in a record; if 0, then first line
const string TxtInFile::LineNumbToStr(Err::eCode code, BYTE lineInd) const
{
	_errCode = code;
	ostringstream ss;
	if(IsFlag(PRNAME))	ss << FileName();
	ss << ": line " << LineNumber(lineInd);
	return ss.str();
}

/************************ TxtInFile: end ************************/

#ifdef _FILE_WRITE

/************************ TxtOutFile ************************/

bool TxtOutFile::Zipped;												// true if filed should be zippped

//TxtOutFile::fAddChar TxtOutFile::fLineAddChar[] = {	// 'Add delimiter' methods
//	&TxtOutFile::AddCharEmpty,	// empty method
//	&TxtOutFile::AddDelim		// adds delimiter and increases current position
//};

// Allocates memory for write line buffer with checking.
//  return: true if successful
bool TxtOutFile::CreateLineBuff(rowlen len)
{
	try { _lineBuff = new char[_lineBuffLen=len]; }
	catch(const bad_alloc)	{ SetError(Err::F_MEM); };
	return IsGood();
}

// Closes adding record to the IO buffer: set current rec position and increases rec counter
//	@len: length of added record
void TxtOutFile::EndRecordToIOBuff(size_t len)
{
	_currRecPos += len;
	_buff[_currRecPos++] = LF;
#ifdef _MULTITHREAD
	if (IsFlag(MTHREAD) && !IsClone())
		InterlockedIncrement(_totalRecCnt);		// atomic increment _recCnt
	else
#endif
		_recCnt++;
}


// Writes nonempty buffer, deletes line buffer and closes file
TxtOutFile::~TxtOutFile()
{
	if(_currRecPos)	Write();
	delete [] _lineBuff;
}

#ifdef _MULTITHREAD
// Constructs a clone of an existing instance.
// Clone is a copy of opened file with its own separate basic I/O and write line buffers.
//	Used for multithreading file recording
//	@file: opened file which is cloned
TxtOutFile::TxtOutFile(const TxtOutFile& file):
	_lineBuffOffset(file._lineBuffOffset),
	_delim(file._delim),
	TxtFile(file)
{
	if( !CreateLineBuff(file._lineBuffLen) )	return;
	memcpy(_lineBuff, file._lineBuff, _lineBuffLen);
	_totalRecCnt = &file._recCnt;
}
#endif

// Adds character to the current position in the line write buffer with optional adding delimiter
//	@ch: char to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddChar(char ch, bool addDelim)
{
	LineAddChar(ch);
	LineAddDelim(addDelim);
}

// Copies block of chars to the current position in the line write buffer,
//	adds delimiter after string	and increases current position.
//	@src: pointer to the block of chars
//	@num: number of chars
//	@addDelim: if true then adds delimiter and increases current position
//	return: new current position
rowlen TxtOutFile::LineAddChars(const char* src, rowlen num, bool addDelim)
{
	memcpy(_lineBuff + _lineBuffOffset, src, num);
	_lineBuffOffset += num;
	LineAddDelim(addDelim);
	return _lineBuffOffset;
}

// Adds integral value to the current position of the line write buffer,
//	adds delimiter after value and increases current position.
//	@v: value to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddInt(LLONG v, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%d", v);
	LineAddDelim(addDelim);
	//LineAddStr(to_string(val), addDelim);
}

// Adds two integral values separated by default delimiter to the current position 
// of the line write buffer, adds delimiter after value and increases current position.
//	@v1: first value to be set
//	@v2: second value to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddInts(ULONG v1, ULONG v2, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%u%c%u", v1, _delim, v2);
	LineAddDelim(addDelim);
	//ostringstream ss;
	//ss << v1 << _delim << v2;
	//LineAddStr(ss.str(), addDelim);
}

// Adds three integral values separated by default delimiter to the current position 
// of the line write buffer, adds delimiter after value and increases current position.
//	@v1: first value to be set
//	@v2: second value to be set
//	@v3: third value to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddInts(ULONG v1, ULONG v2, ULONG v3, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%u%c%u%c%u", v1, _delim, v2, _delim, v3);
	LineAddDelim(addDelim);
	//ostringstream ss;
	//ss << v1 << _delim << v2 << _delim << v3;
	//LineAddStr(ss.str(), addDelim);
}

// Adds floating point value to the current position of the line write buffer,
//	adds delimiter after value and increases current position.
//	@val: value to be set
//	@ndigit: number of digits to generate
//	@addDelim: if true then adds delimiter and increases current position
//void TxtOutFile::LineAddFloat(float val, BYTE ndigit, bool addDelim)
//{
//	// double val, because casting int to float is not precise, f.e. (float)61342430 = 61342432
//	_gcvt(val, ndigit, _lineBuff + _lineBuffOffset);
//	_lineBuffOffset += rowlen(strlen(_lineBuff + _lineBuffOffset));
//	LineAddDelim(addDelim);
//}

// Adds line to IO buffer (from 0 to the current position).
//	@offset: start position for the next writing cycle
void TxtOutFile::LineToIOBuff(rowlen offset)
{
	RecordToIOBuff(_lineBuff, _lineBuffOffset);
	_lineBuffOffset = offset;
}

// Adds record to the the IO buffer.
// Generate exception if writing is fall.
//	@src: record
//	@len: length of record
void TxtOutFile::RecordToIOBuff(const char *src, size_t len)
{
	if( _currRecPos + len + 1 > _buffLen )	// write buffer to file if it's full
		Write();
	memcpy(_buff + _currRecPos, src, len);
	EndRecordToIOBuff(len);
}

// Adds string to IO buffer without checking buffer exceeding.
void TxtOutFile::StrToIOBuff(const string&& str)
{
	memmove(_buff + _currRecPos, str.c_str(), str.length());
	EndRecordToIOBuff(str.length());
}

// Allocates line write buffer.
//	@len: length of buffer
void TxtOutFile::SetLineBuff(rowlen len)
{
	CreateLineBuff(len);
	memset(_lineBuff, _delim, len);
}

// Adds to line 2 int values and adds line to the IO buff.
//void TxtOutFile::WriteLine(ULONG val1, ULONG val2)
//{
//	LineAddInts(val1, val2, false);
//	LineToIOBuff();
//}

// Adds to line C-string with delimiter, and adds line to the IO buff.
//void TxtOutFile::WriteLine(const char* str)
//{
//	LineAddStr(str, strlen(str));
//	LineToIOBuff();
//}

// Adds to line string with delimiter and int value, and adds line to the IO buff.
//void TxtOutFile::WriteLine(const string& str, int val)
//{
//	LineAddStr(str);
//	LineAddInt(val, false);
//	LineToIOBuff();
//}

// Writes thread-safely current block to file.
void TxtOutFile::Write() const
{
	//_stopwatch.Start();
#ifdef _MULTITHREAD
	const bool lock = IsFlag(MTHREAD) && Mutex::IsReal(_mtype);
	if(lock)	Mutex::Lock(_mtype);
#endif
	int res = 
#ifndef _NO_ZLIB
		IsZipped() ?
		gzwrite((gzFile)_stream, _buff, _currRecPos) :
#endif
		fwrite(_buff, 1, _currRecPos, (FILE*)_stream);
	if(res == _currRecPos)	_currRecPos = 0;
	else { SetError(Err::F_WRITE); /*cout << "ERROR!\n";*/ }
#ifdef _MULTITHREAD
	if(lock) {
		if(IsClone())
			//InterlockedExchangeAdd(_totalRecCnt, _recCnt);
			*_totalRecCnt += _recCnt,
			_recCnt = 0;
		Mutex::Unlock(_mtype);
	}
#endif
	//_stopwatch.Stop();
}

//bool	TxtOutFile::AddFile(const string fName)
//{
//	TxtFile file(fName, *this);
//	while( file.ReadBlock(0) ) {
//		_currRecPos = file._readedLen;
//		if( !Write() )	return false;
//	}
//	return true;
//}

/************************ class TxtOutFile: end ************************/

#endif	// _FILE_WRITE

#ifndef _FQSTATN

/************************ class TabFile ************************/

// Checks if field valid and throws exception if not.
//	@ind: field index
//	return: true if field is valid
bool TabFile::IsFieldValid	(BYTE ind) const
{
	if(_fieldPos[ind] == vUNDEF) {		// vUNDEF was set in buff if there was no such field in the line
	//|| !StrField(ind)[0]) {			// empty field
		if(ind < FT::FileParams(_fType).MinFieldCnt)
			Err(Err::TF_FIELD, LineNumbToStr(Err::TF_FIELD).c_str()).Throw();
		return false;
	}
	return true;
}

// Initializes new instance.
//	@type: file bioinfo type
//	@isEstLineCnt: if true then estimate count of lines
void TabFile::Init(FT::eType type, bool isEstLineCnt)
{
	if(!_fieldPos)	
		_fieldPos = new short[FT::FileParams(type).MaxFieldCnt];
	_lineSpecLen = FT::FileParams(type).LineSpecLen();

	// estimate line cnt
	if(isEstLineCnt)
		if (FT::FileParams(type).AvrLineLen)
			SetEstLineCount(type);
		else if (GetNextLine(true))
			SetEstLineCount();
		else _estLineCnt = 0;
	//cout << "  file estLen: " << _estLineCnt << LF;
}

// Sets estimation of number of lines by currently readed line, and rolls line back
void TabFile::SetEstLineCount()
{
	_estLineCnt = ULONG(Length() / RecordLength());
	RollBackLine();
}

// Initializes instance by a new type, correct estimated number of lines if it's predefined
void TabFile::ResetType(FT::eType type)
{
	if (FT::FileParams(_fType).MaxFieldCnt < FT::FileParams(type).MaxFieldCnt)
		Release();
	Init(_fType = type, true);
}

// Returns a pointer to the substring defined by key.
//	@str: null-terminated string to search the key
//	@key: string to search for
//	return: a pointer to the substring after key, or NULL if key does not appear in str
const char* TabFile::KeyStr(const char* str, const string& key)
{
	if (str) {
		const char* strKey = strstr(str, key.c_str());
		if (strKey)	return strKey + key.length();
	}
	return NULL;
}

// Checks definition or declaration line for key
//	@str: null-terminated string to search the key
//	@key: string to search for
//	return: point to substring followed after the key
const char* TabFile::CheckSpec(const char* str, const string& key)
{
	const char* strKey = KeyStr(str, key);
	if (!strKey)
		ThrowExcept(LineNumbToStr() + ": absent or wrong '" + key + "' key");
	return strKey;
}

// Skip commented lines and returns estimated number of uncommented lines
//ULONG TabFile::GetUncommLineCount()
//{
//	const char*	line;
//	const char comment = FT::FileParams(_fType).Comment;
//	ULONG	cnt = 0;
//	
//	for(USHORT pos=0; line = GetNextRecord(); pos=0) {
//		while( *(line+pos)==SPACE )	pos++;		// skip blanks at the beginning of line
//		if (*(line+pos) != comment)	break;		// skip comment line
//	}
//	if(line) {
//		cnt = ULONG(Length() / RecordLength());
//		RollBackLastRecord();
//	}
//	return cnt;
//}

// Reads first line and set it as current.
// Throws exception if invalid
//	@cntLines: returned estimated count of lines.
//	It works properly only if lines are sorted by ascending, f.e. in sorted bed-files.
//	return: current line
//const char* TabFile::GetFirstLine(ULONG& cntLines)
//{
//	cntLines = 0;
//	if(GetNextLine() )
//		cntLines = (ULONG)(Length() / RecordLength());	
//	return _currLine;
//}

// Reads next line and set it as current.
//	@checkTabs: it true then check the number of incoming fields (tabs)
//		An invoke with the false is performed only for WIG files, to avoid checking the declaration line
//	return: current line.
const char*	TabFile::GetNextLine(bool checkTab)
{
	const TabFilePar par = FT::FileParams(_fType);
	
	// fill _fieldPos 0 to check if all fields will be initialize
	// _fieldPos keeps start positions (next after TAB) for each field in line.
	// First field's position usually is 0, but if there are blanks in the beginning of the line,
	// it will be in first unblank position.
	// GetNextRecord(..) fills _fieldPos beginning from second field.
	memset(_fieldPos, vUNDEF, sizeof(short) * par.MaxFieldCnt);

	char* line = GetNextRecord(_fieldPos, par.MaxFieldCnt);	// a bit faster than using _currLine
	if(line) {
		USHORT	currPos = 0;	// a bit faster than using _currPos in heep
	
		// skip blanks at the beginning of line
		while (line[currPos] == SPACE)	currPos++;
		// skip comment line or line without specifier
		if (*(line + currPos) == par.Comment
		|| (par.LineSpec && memcmp(line + currPos, par.LineSpec, _lineSpecLen)))
			return GetNextLine(checkTab);

		_fieldPos[0] = currPos;		// set start position of first field
#ifdef _WIG
		if (checkTab || isdigit(line[0]))
#endif // _WIG
			for (BYTE i = 1; i < par.MaxFieldCnt; line[_fieldPos[i++] - 1] = cNULL)	// replace TABs by 0
				if (_fieldPos[i] == vUNDEF)		// numbers of TABs is less than number of fields
					if (i >= par.MinFieldCnt)				// less than number of optional fields
						break;
					else {									// less than number of required fields
						SetError(Err::TF_FIELD, LineNumbToStr(),
							": " + to_string(i) + " against " + to_string(par.MinFieldCnt) + "; wrong format?");
						return _currLine = NULL;
					}
		line[RecordLength() - 1] = cNULL;	// replace '\n' by 0
		
		//cout << this->FileName() << "\trecLen = " << RecordLength() << LF;

		// this version when _fieldPos is not initialized in TxtFile::GetNextRecord(): a bit more slower
		//for(BYTE i=1; i<_cntFields; i++) {
		//	// search for the next TAB
		//	for(; currPos<_recLen; currPos++)
		//		if(line[currPos]==TAB)
		//			break;
		//	line[currPos] = cNULL;		// replace TABs by 0
		//	_fieldPos[i] = ++currPos;
		//}
	}
	return _currLine = line;
}

/************************ end of TabFile ************************/

#ifndef _WIGREG
/************************ ChromDefRegions ************************/

const string ChromDefRegions::Ext = ".region";	// regions file extension

// Returns false if the next region is separated from the current one
// by a value less than the established minimum gap.
// Otherwise returns true and combined region as parameter.
bool ChromDefRegions::Combiner::ExceptRegion(Region& rgn)
{
	if(_rgn.Empty())				 { _rgn = rgn;			return false; }	// first income region
	if(rgn.Start-_rgn.End < _gapLen) { _rgn.End = rgn.End;	return false; }
	Region	next = rgn;
	rgn = _rgn;		// returned region
	_rgn = next;
	return true;
}

// Creates new instance and initializes it from file if one exist
//	@fName: chrom file name without extension
//	@minGapLen: length, gaps less than which are ignored when reading; if 0 then read all regions 
ChromDefRegions::ChromDefRegions(const string& fName, chrlen minGapLen) : _gapLen(0), _new(true)
{
	_fName = fName + Ext;
	if(FS::IsFileExist(_fName.c_str())) {
		TabFile file(_fName, FT::eType::RGN);
		Combiner comb(minGapLen);
		Region rgn;

		Reserve(file.EstLineCount());
		if (file.GetNextLine()) {
			_gapLen = file.LongField(1);	// first line contains total length of gaps
			while (file.GetNextLine()) {
				rgn.Set(file.LongField(0), file.LongField(1));
				if (!minGapLen || comb.ExceptRegion(rgn))
					Add(rgn);
			}
			if (minGapLen)	Add(comb.LastRegion());
			_new = false;
		}
	}
	else	Reserve(DefCapacuty);
}

// Adds def region beginning of current gap position.
//	@rgn: def region
//	@minGapLen: minimal length which defines gap as a real gap
void ChromDefRegions::AddRegion(const Region& rgn, chrlen minGapLen)
{
	if( rgn.End && minGapLen && rgn.Length())	// current gap is closed
		if( _regions.size() && rgn.Start-_regions.back().End < minGapLen )
			_regions.back().End = rgn.End;		// pass minimal allowed gap
		else 
			_regions.push_back(rgn);			// add new def-region
}

void ChromDefRegions::Write() const
{
	if(!_new || FS::IsShortFileName(_fName))	return;
	ofstream file;

	file.open (_fName.c_str(), ios_base::out);
	file << "SumGapLen:" << TAB << _gapLen << LF;
	for(Iter it=Begin(); it!=End(); it++)
		file << it->Start << TAB << it->End << LF;
	file.close();
	_new = false;
}

#if defined _READDENS || defined _BIOCC
// Combines regions with a gap less than minGapLen
void ChromDefRegions::Combine(chrlen minGapLen)
{
	Region rgn;
	Regions rgns;	// new combined regions
	Combiner comb(minGapLen);

	rgns.Reserve(Count());
	for(Iter it=Begin(); it!=End(); it++) {
		rgn.Set(it->Start, it->End);
		if(comb.ExceptRegion(rgn))	rgns.Add(rgn);
	}
	rgns.Add(comb.LastRegion());
	Copy(rgns);
}
#endif

/************************ end of ChromDefRegions ************************/

/************************ class FaFile ************************/

// Adds a gap to assign defined regions
//	@start: gap's start position.
//	@len: gap's length
void FaFile::DefRgnMaker::AddGap(chrlen start, chrlen len)
{
	_defRgn.End = (start += _currPos);
	_defRgns.AddRegion(_defRgn, _minGapLen);
	_defRgn.Start = start + len;
	_defRgns.IncrGapLen(len);
}

// Closes adding gaps, saves hrom's defined regions
//	@cLen: chrom length
void FaFile::DefRgnMaker::CloseAddGaps(chrlen cLen)
{
	_defRgn.End = cLen;
	_defRgns.AddRegion(_defRgn, _minGapLen);
	_defRgns.Write(); 
}

// Search 'N' subsequence in the current line beginning from the start position
// and adds complete 'N' subsequence to _rgnMaker
//	@startPos: starting line reading position
//	@NCnt: number of 'N' in the line beginning from the start position
//	return: true if line trimming is complete, false for single 'N'
void FaFile::CountN(chrlen startPos, chrlen NCnt)
{
	if(NCnt == 1)	return;

	const char* line = Line();
	const chrlen lineLen = LineLength();
	chrlen nCnt = 0, iN;				// local count N, first N index
	bool firstN = true;

	for(chrlen i=startPos; i<lineLen; i++)
		if(line[i] == cN) {					// some 'N' subsequence begins?
			nCnt++;
			if(firstN)	{ iN = i; firstN = false; }
		}		
		else if(nCnt)						// is current subsequence finished?
			if(nCnt == 1) {					// was it a single 'N'?
				nCnt = 0; NCnt--; firstN = true;	
			}
			else {							// continuous 'N' sequence?
				_rgnMaker->AddGap(iN, nCnt);	// add complete 'N' subsequence
				CountN(i++, NCnt - nCnt);	// search in the rest of line
				break;
			}
}

// Reads line and set it as current, and generates regions
// Since method scans the line, it takes no overhead expenses to do this.
//	@pocket: external temporary variables
//	return: current line.
const char* FaFile::GetLineWitnNControl() 
{
	chrlen nCnt = 0;	// number of 'N' in line
	const char* line = GetNextRecord(nCnt);
	if(line) {
		const chrlen len = LineLength();
		if(nCnt)
			if(nCnt == len)	_rgnMaker->AddGap(0, nCnt);	// the whole line is filled by 'N'
			else			CountN(0, nCnt);
		_rgnMaker->AddLineLen(len);
	}
	return line;
}

// Opens existing file and reads first line.
//	@rgns: def regions to fill, otherwise NULL to reading without 'N' control
FaFile::FaFile(const string& fName, ChromDefRegions* rgns) : TxtInFile(fName, eAction::READ, 1)
{
	if (rgns) {
		_rgnMaker.reset(new DefRgnMaker(*rgns, 2));
		_pGetLine = &FaFile::GetLineWitnNControl;
	}
	else _pGetLine = &FaFile::GetNextRecord;

	chrlen len = chrlen(Length());
	if(NextGetLine()[0] == FaComment) {		// is first line a header?
		len -= RecordLength();
		if(_rgnMaker)	_rgnMaker->RemoveLineLen(LineLength());
		NextGetLine();
	}
	// set chrom length
	_cLen = len - LFSize() *
		(len/RecordLength() +		// amount of LF markers for whole lines
		bool(len%RecordLength()));	// LF for part line
}

/************************ end of class FaFile ************************/

#endif	// _ISCHIP, _READDENS, _BIOCC
#endif	// no _FQSTATN
#if defined _CALLDIST || defined _FQSTATN

/************************ FqFile ************************/

// Returns checked length of current readed Read.
readlen FqFile::ReadLength() const
{
	//CheckGettingRecord();
	return LineLengthByInd(READ); 
}

// Gets checked Read from readed Sequence.
const char* FqFile::GetCurrRead() const
{ 
	//CheckGettingRecord();
	return NextRecord() - RecordLength() + LineLengthByInd(HEADER1); 
}

// Returns checked Sequence
const char* FqFile::GetSequence()
{
	const char* record = GetNextRecord();
	if(record != NULL) {
		if(*record != AT)
			Err("non '@' marker; missed header line", LineNumbToStr()).Throw();
		if( *(record + LineLengthByInd(HEADER1, false) + LineLengthByInd(READ, false)) != PLUS )
			Err("non '+' marker; missed second header line", LineNumbToStr()).Throw();
	}
	return record;
}
/************************ FqFile: end ************************/

#endif	// _FQSTATN

// Creates new instance with read buffer belonges to aggregated file: constructor for concatenating.
//	For reading only.
//	@fName: valid full name of file
//	@file: file that is aggregated
//TxtFile::TxtFile(const string & fName, const TxtFile& file) :
//	_flag(file._flag),
//	_recLineCnt(file._recLineCnt),
//	_buffLineLen(0)
//{
//	if( !SetBasic(fName, READ, NULL) )	return;
//	_buffLen = file._buffLen;
//	_buff = file._buff;
//	RaiseFlag(CONSTIT);	
//}