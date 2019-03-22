/**********************************************************
TxtFile.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
Last modified: 21.03.2019
Provides read|write text file functionality
***********************************************************/

#include "TxtFile.h"
#ifndef _FILE_WRITE
#include <fstream>
#endif

/************************ class FT ************************/

const FT::fType FT::Types[] = {
	{ "", strEmpty, strEmpty,			TabFilePar( 0, 0, '\0', NULL) },		// undefined type
	{ "bed", "read",	"reads",		TabFilePar( 6, 6, HASH, Chrom::Abbr) },	// alignment bed
	{ "bed", "feature", "features",		TabFilePar( 3, 6, HASH, Chrom::Abbr) },	// ordinary bed
	{ "wig", "interval", "intervals",	TabFilePar( 2, 2, HASH, NULL) },		// wiggle
	{ "sam", strEmpty, strEmpty,		TabFilePar( 0, 0, HASH, NULL) },		// sam
	{ "fq", strEmpty, strEmpty,			TabFilePar( 0, 0, '\0', NULL) },		// fastQ
	{ "fa", strEmpty, strEmpty,			TabFilePar( 0, 0, '\0', NULL) }			// fasta
};
const BYTE FT::Count = sizeof(FT::Types)/sizeof(FT::fType);

// Returns file format
//	@fName: file name (with case insensitive extension)
FT::eTypes FT::GetType(const char* fName)
{
	const string ext = FS::GetExt(fName);
	const char* c_ext = ext.c_str();
	for(int i = 2; i<Count; i++)	// start from ordinary bed
		if(!_stricmp(c_ext, Types[i].Extens))	return eTypes(i);
	return UNDEF;
}

// Validates file extension
//	@fName: file name (with case insensitive extension and [.gz])
//	@t: file type
//	@printfName: true if file name should be ptinted
//	@throwExc: true if throw exception, otherwise throw warning
//	return: true if file extension correspondes to file type
bool FT::CheckType(const char* fName, eTypes t, bool printfName, bool throwExc) { 
	if( GetType(fName) != (t == ABED ? BED : t) ) {
		Err("wrong extension", printfName ? fName : NULL).Throw(throwExc);
		return false;
	}
	return true;
}

#ifdef _ISCHIP
// Gets file extension, beginning at DOT and adding .gz if needed
//	@t: file type
//	@isZip: true if add ".gz"
const string FT::RealExt(eTypes t, bool isZip)
{
	string ext = string(".") + Types[t].Extens;
	if(isZip)	ext += ZipFileExt;
	return ext;
}
#endif

/************************ end of class FT ************************/

/************************ TxtFile ************************/
const char* modes[] = { "r", "w", "a+" };
const char* bmodes[] = { "rb", "wb" };

//int TxtFile::NUMB_BLK = 4;	// size of basic block in Mb

// Sets error code and throws exception if it is allowed.
void TxtFile::SetError(Err::eCode errCode) const
{
	_errCode = errCode;
	if( IsFlag(ABORTING) )
		Err(errCode, FileNameToExcept().c_str()).Throw();
}

// Initializes instance variables, opens a file, sets a proper error code.
//	@fName: valid full name of file
//	@mode: opening mode
//	@flStream: clonable file stream or NULL
//	return: true if success, otherwise false.
bool TxtFile::SetBasic(const string& fName, eAction mode, void* flStream)
{
	_buff = NULL;
	_stream = NULL;
	_errCode = Err::NONE;
	_fName = fName;
	_currRecPos = _recCnt = 0;
	_buffLen = NUMB_BLK * BasicBlockSize;	// by default; can be corrected
#ifdef _NO_ZLIB
	if(IsZipped()) { SetError(Err::FZ_BUILD); return false; }
#endif
	// set file stream
	if(flStream)	_stream = flStream;
	else {
#ifndef _NO_ZLIB
		if(IsZipped())
			if(mode == READ_ANY)	SetError(Err::FZ_OPEN);
			else	_stream = gzopen(fName.c_str(), bmodes[mode]);
		else
#endif
			_stream = fopen(fName.c_str(), modes[mode]);
		if( _stream == NULL )	SetError(Err::F_OPEN);
	}
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
//	@abortInvalid: true if invalid instance shold be completed by throwing exception
//	@rintName: true if file name should be printed in exception's message
TxtFile::TxtFile (const string& fName, eAction mode, bool abortInvalid, bool printName) :
	_flag(1)	// EOL size set to 1
{
	SetFlag(ZIPPED, FS::HasGzipExt(fName));
	SetFlag(ABORTING, abortInvalid);
	SetFlag(PRNAME, printName);
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
		if( mode != WRITE )	
			_buffLen = (ULONG)_fSize + 1;	// for reading
		else if( _fSize * 2 < _buffLen )
			_buffLen = (ULONG)_fSize * 2;	// for writing: increase small buffer for any case
	
	if( !CreateIOBuff() )	return;

#ifdef ZLIB_NEW
	if( IsZipped() && gzbuffer( (gzFile)_stream, _buffLen) == -1 )
	{ SetError(Err::FZ_MEM); return; }
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
	if( !SetBasic("Clone " + file._fName, WRITE, file._stream) )	return;
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
//	@abortInvalid: true if invalid instance should be completed by throwing exception
//	@rintName: true if file name should be printed in an exception's message
TxtInFile::TxtInFile(const string& fName, eAction mode, 
	BYTE cntRecLines, bool abortInvalid, bool printName) : 
	_linesLen(NULL),
	_linesInRecCnt(cntRecLines),
	_readedLen(0),
	TxtFile(fName, mode, abortInvalid, printName)
{
	if(ReadBlock(0) >= 0)					// read first block
		_linesLen = new UINT[cntRecLines];	// set lines buffer
}

// Reads next block.
//	@offset: shift of start reading position
//	return: 1 if file is not finished; 0 if it is finished; -1 if unsuccess reading
int TxtInFile::ReadBlock(const UINT offset)
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
	return _readedLen == 0 ? 0 : 1;
}

// Sets _currLinePos to the beginning of next non-empty line in _buff
//	@counterN: if not NULL, adds to counterN the number of 'N' in a record. Used in Fa() only.
//	@posTab: if not NULL, sets TABs positions in line to this array
//	@cntTabs: if @posTab is not NULL, the length of @posTab array
//	return: point to line or NULL if no more lines
char* TxtInFile::GetRecord(chrlen* const counterN, short* const posTab, BYTE cntTabs)
{
	if(Length() == 0)	return ReadingEnded();
	UINT i, 
		 currLinePos = _currRecPos,	// start position of current readed record
		 cntEmpty = 0;				// counter of empty lines
	chrlen cntN = 0;				// local counter of 'N'
	BYTE indTab = 1;				// index of TAB position in @posTab
	char c;

	_recLen = 0;
	for(BYTE r=0; r<_linesInRecCnt; r++)
		for(i=currLinePos;; i++)	{
			if( i >= _readedLen ) {				// check for the next block
				if( _readedLen != _buffLen )	// final block
					if( i == _buffLen			// current line is not ended by EOL
					|| _buff[i-1] == EOL		// current line which is ended by EOL
					|| _buff[i-1] == '\0')		// EOL has been replaced by previous call TabFile::GetLine()
						return ReadingEnded();	// normal completion
					else {						// current line is not ended by EOL
						if(i<_buffLen)
							_buff[i] = '\0';	// close current line
						goto A;
					}
				// jump to the next block
				if( !currLinePos ) {			// this block is totally unreaded
					SetError(Err::F_BIGLINE);
					return ReadingEnded();
				}
				// length of untreated rest of record
				UINT restLen = _readedLen - _currRecPos - cntEmpty;
				// move untreated remainder to the beginning of buffer; if restLen=0, skip moving
				memmove(_buff, _buff + _currRecPos, restLen);
				if( !ReadBlock(restLen) )
					return ReadingEnded();	// record was finished exactly by previous ReadBlock()
				// restart at the beginning of the record
				indTab = 1;
				cntN = currLinePos = cntEmpty = i = _recLen = r = 0;
			}
			c = _buff[i];
			if( c==cN )
				cntN++;
			else if( c==TAB ) {
				if(indTab < cntTabs)
					posTab[indTab++] = i - currLinePos + 1;
			}
			else if( c==EOL ) {
				if( i==currLinePos ) {			// EOL marker is first in line
					currLinePos++; cntEmpty++;	// skip empty line
					continue;
				}
				if( IsEOLnotDef() )		SetEOL(_buff[i-1]);		// define EOL size
A:				_recLen += (_linesLen[r] = i - currLinePos + 1);
				currLinePos = i + 1;
				break;
			}
		}
	_currRecPos = currLinePos;			// next record position
	 if(counterN)		(*counterN) += cntN;
	 _recCnt++;
	return _buff + _currRecPos - _recLen;
}

// Gets string containing file name and current line number.
//	@lineInd: index of line in a record; if 0, then first line
const string TxtInFile::LineNumbToStr(BYTE lineInd) const
{
	string str;
	if(IsFlag(PRNAME))	str = FileName();
	//str += SepCl;
	return str + ": line " + NSTR(LineNumber(lineInd));
}


/************************ TxtInFile: end ************************/

#ifdef _FILE_WRITE

/************************ TxtOutFile ************************/

// Allocates memory for write line buffer with checking.
//  return: true if successful
bool TxtOutFile::CreateLineBuff(rowlen len)
{
	try { _lineBuff = new char[_lineBuffLen=len]; }
	catch(const bad_alloc)	{ SetError(Err::F_MEM); };
	return IsGood();
}

// Writes nonempty buffer, deletes line buffer and closes file
TxtOutFile::~TxtOutFile()
{
	if(_currRecPos)	Write();
	if( _lineBuff )	delete [] _lineBuff;
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

// Copies block of chars to the current position in the line write buffer,
//	adds delimiter after string	and increases current position.
//	@src: pointer to the block of chars
//	@num: number of chars
//	@addDelim: if true then adds delimiter and increases current position
//	return: new current position
rowlen TxtOutFile::LineAddChars(const char* src, size_t num, bool addDelim)
{
	memcpy(_lineBuff+_lineBuffOffset, src, num);
	_lineBuffOffset += num;
	if(addDelim)	LineAddDelim();
	return _lineBuffOffset;
}

// Adds integral value to the current position of the line write buffer,
//	adds delimiter after value and increases current position.
//	@val: value to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddInt(LLONG val, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%d", val);
	if(addDelim)	LineAddDelim();
}

// Adds two integral values separated by default delimiter to the current position 
// of the line write buffer, adds delimiter after value and increases current position.
//	@val1: first value to be set
//	@val2: second value to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddInts(ULONG val1, ULONG val2, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%d%c%d", val1, _delim, val2);
	if(addDelim)	LineAddDelim();
}

// Adds three integral values separated by default delimiter to the current position 
// of the line write buffer, adds delimiter after value and increases current position.
//	@val1: first value to be set
//	@val2: second value to be set
//	@val3: third value to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddInts(ULONG val1, ULONG val2, ULONG val3, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%d%c%d%c%d", val1, _delim, val2, _delim, val3);
	if(addDelim)	LineAddDelim();
}

// Adds floating point value to the current position of the line write buffer,
//	adds delimiter after value and increases current position.
//	@val: value to be set
//	@ndigit: number of digits to generate
//	@addDelim: if true then adds delimiter and increases current position
void TxtOutFile::LineAddFloat(float val, BYTE ndigit, bool addDelim)
{
	// double val, because casting int to float is not precise, f.e. (float)61342430 = 61342432
	_gcvt(val, ndigit, _lineBuff + _lineBuffOffset);
	_lineBuffOffset += strlen(_lineBuff + _lineBuffOffset);
	if(addDelim)	LineAddDelim();
}

// Adds line to IO buffer (from 0 to the current position).
//	@offset: start position for the next writing cycle
//	@closeLine: if true then close line by EOL
void TxtOutFile::LineToIOBuff(rowlen offset, bool closeLine)
{
	RecordToIOBuff(_lineBuff, _lineBuffOffset, closeLine);
	_lineBuffOffset = offset;
}

// Adds record to the the IO buffer.
// Generate exception if writing is fall.
//	@src: record
//	@len: length of record
//	@closeLine: if true then close line by EOL
void TxtOutFile::RecordToIOBuff	(const char *src, UINT len, bool closeLine)
{
	if( _currRecPos + len + 1 > _buffLen )	// write buffer to file if it's full
		Write();
	memcpy(_buff + _currRecPos, src, len );
	_currRecPos += len;
	if(closeLine)	_buff[_currRecPos++] = EOL;
#ifdef _MULTITHREAD
	if(IsFlag(MTHREAD) && !IsClone())
		InterlockedIncrement(_totalRecCnt);		// atomic increment _recCnt
	else	
#endif
		_recCnt++;
}

// Allocates line write buffer.
//	@len: length of buffer
//	@delim: delimiter between line fields
void TxtOutFile::SetLineBuff(rowlen len)
{
	CreateLineBuff(len);
	memset(_lineBuff, _delim, len);
}

// Adds to line int and float values and adds line to the IO buff.
//	@ndigit: number if float digits to write
//void TxtOutFile::WriteLine(int val1, float val2, BYTE ndigit)
//{
//	LineAddInt(val1);
//	LineAddFloat(val2, ndigit, false);
//	LineToIOBuff();
//}

// Adds to line int and float values and adds line to the IO buff.
//	@ndigit: number if float digits to write
void TxtOutFile::WriteLine(int val1, float val2, float val3, BYTE ndigit, const char*str)
{
	LineAddInt(val1);
	LineAddFloat(val2, ndigit);
	LineAddFloat(val3, ndigit);
	LineAddStr(str);
	LineToIOBuff();
}

// Adds to line 2 int values and adds line to the IO buff.
void TxtOutFile::WriteLine(ULONG val1, ULONG val2)
{
	LineAddInts(val1, val2, false);
	LineToIOBuff();
}

// Adds to line C string and adds line to the IO buff.
void TxtOutFile::WriteLine(const char* str)
{
	LineAddStr(str, strlen(str));
	LineToIOBuff();
}

// Adds to line string and int value and adds line to the IO buff.
void TxtOutFile::WriteLine(const string& str, int val)
{
	LineAddStr(str);
	LineAddInt(val, false);
	LineToIOBuff();
}

// Writes thread-safely current block to file.
void TxtOutFile::Write() const
{
	//_stopwatch.Start();
#ifdef _MULTITHREAD
	if(IsFlag(MTHREAD))	Mutex::Lock(Mutex::WR_FILE);
#endif
	//cout << "Write " << FileName() << "\trecords = " << _recCnt;// << endl;
	int res = 
#ifndef _NO_ZLIB
		IsZipped() ?
		gzwrite((gzFile)_stream, _buff, _currRecPos) :
#endif
		fwrite(_buff, 1, _currRecPos, (FILE*)_stream);
		if(res != _currRecPos)	{ SetError(Err::F_WRITE); cout << "ERROR!\n"; }
	else					_currRecPos = 0;
#ifdef _MULTITHREAD
	if(IsFlag(MTHREAD)) {
		//cout << "\tUNLOCK\n";
		Mutex::Unlock(Mutex::WR_FILE);
		if(IsClone()) {
			InterlockedExchangeAdd(_totalRecCnt, _recCnt);
			_recCnt = 0;
		}
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

#ifndef _WIGREG

/************************ struct Region ************************/

// Extends Region with chrom length control.
// If extended Region starts from negative, or ends after chrom length, it is fitted.
//	@extLen: extension length in both directions
//	@cLen: chrom length; if 0 then no check
void Region::Extend(chrlen extLen, chrlen cLen)
{
	Start -= extLen > Start ? Start : extLen;
	End += extLen;
	if(cLen && End > cLen)	End = cLen;
}

/************************ end of struct Region ************************/

/************************ class Regions ************************/

// Geta total length of regions.
chrlen Regions::Length () const
{
	chrlen len = 0;
	for(Iter it=_regions.begin(); it<_regions.end(); it++)
		len += it->Length();
	return len;
}

#if defined _DENPRO || defined _BIOCC

// Returns an iterator referring to the past-the-end element, where end is external
//	@curr_it: region's const iterator, from which the search is started
//	@end: external pre-defined end coordinate
Regions::Iter Regions::ExtEnd(Iter curr_it, chrlen end) const
{
	Iter it = curr_it;
	for(; it != _regions.end(); it++)
		if(it->End > end)	break;
	return it;
}

	// Initializes this instance by intersection of two Regions.
//	Typically for that purpose is used Interval Tree,
//	but this implementation uses Regions ordering and is simpler.
void Regions::FillOverlap(const Regions &regn1, const Regions &regn2)
{
	chrlen start=0, end=0, start1, start2, end1, end2;
	Iter it1 = regn1._regions.begin();
	Iter it2 = regn2._regions.begin();
	Reserve( max(regn1.Count(), regn2.Count()) );
	for(; it1 != regn1._regions.end(); it1++) {
		start1 = it1->Start;	end1 = it1->End;
		start2 = it2->Start;	end2 = it2->End;
		if(start1 < start2) {
			if(end1 > start2) {
				start = start2;
				if(end1 > end2) {	end = end2;	it2++; it1--; }
				else				end = end1;
			}
		}
		else
			if(start1 >= end2)	{	it2++; it1--;	}
			else {
				start = max(start1, start2);
				if(end1 > end2) {	end = end2;	it2++; it1--; }
				else				end = end1;
			}
		if( end ) {
			AddRegion(start, end);
			if( it2 == regn2._regions.end() )		return;
			end = 0;
		}
	}
}

// Initializes this instance by inverted external Regions,
// so the regions turns to the gaps and vice versa.
// Each new region is less than proper old gap by 1 from each side.
//	@regn: external Regions
//	@maxEnd: the maximum possible end-coordinate of region:
//	the chromosome length in case of nucleotides sequance.
void Regions::FillInvert(const Regions &regn, chrlen maxEnd)
{
	chrlen start=0, end=0;
	Reserve(regn.Count() + 1);
	Iter it = regn._regions.begin();
	for(; it != regn._regions.end(); it++) {
		end = it->Start - 1;
		AddRegion(start, end);
		start = it->End + 1;
	}
	if(start < (maxEnd + 1))
		AddRegion(start, maxEnd);
}

// Reads data from file @fname
//	return: written minimal gap length
short Regions::Read(const string & fName)
{
	TabFile file(fName, TxtFile::READ, 2);
	ULONG evalLineLen = 0;
	// first line is definition line with gap length;
	// no need to check since aborting invalid file is set
	const char* currLine = file.GetFirstLine(&evalLineLen);
	short minGapLen = file.IntField(1);
	Reserve( static_cast<chrlen>( file.Length() / evalLineLen ) );
	//do		AddRegion( file.IntField(0), file.IntField(1) );
	while( (currLine = file.GetLine()) != NULL )
		AddRegion( file.IntField(0), file.IntField(1) );
	return minGapLen;
}

// Saves data to file
//	@fname: name of file
//	@minGapLen: minimal length which defines gap as a real gap
void Regions::Write(const string & fName, short minGapLen) const
{
#ifdef _FILE_WRITE
	TxtOutFile file(fName);

	file.SetLineBuff(2*INT_CAPACITY+2);
	file.WriteLine("minGap", minGapLen);
	for(Iter it = Begin(); it != End(); it++)
		file.WriteLine(it->Start, it->End);
	file.Write();
#else
	ofstream file;

	file.open (fName);
	file << "minGap\t" << minGapLen << EOL;
	for(Iter it = Begin(); it != End(); it++)
		file << it->Start << TAB << it->End << EOL;
	file.close();
#endif
}

#endif	// _DENPRO, _BIOCC

// Adds gap, beginning of current gap position.
//	@gapStart: gap's start position.
//	@currGapStart: current gap start position holder
//	@minGapLen: minimal length which defines gap as a real gap
void Regions::AddGap(chrlen gapStart, chrlen currGapStart, chrlen minGapLen)
{
	if( gapStart && gapStart != currGapStart )		// current N-region is closed
		if( currGapStart-_regions.back().End/*-1*/ < minGapLen )
			_regions.back().End = gapStart/*-1*/;		// pass minimal allowed undefined nt
		else 
			_regions.push_back(Region(currGapStart, gapStart));	// add new def-region
}

#ifdef DEBUG
void Regions::Print () const
{
	int i=0;
	cout << "Regions:\n";
	//for(auto it=_regions.begin(); it<_regions.end(); it++)
	for(Iter it=_regions.begin(); it<_regions.end(); it++)
		cout <<
		//setw(2) << 
		++i << COLON << TAB <<
		//setw(9) << 
		it->Start << TAB << it->End << endl;
}
#endif
/************************ end of class Regions ************************/
#endif	// _WIGREG

/************************ class TabFile ************************/

// Checks if field valid and throws exception if not.
//	@ind: field index
//	return: true if field is valid
bool TabFile::IsFieldValid	(BYTE ind) const
{
	if(_fieldPos[ind] == vUNDEF) {		// vUNDEF was set in buff if there was no such field in the line
	//|| !StrField(ind)[0]) {			// empty field
		if(ind < _params.MinFieldCnt)
			ThrowLineExcept(Err::TF_FIELD);
		return false;
	}
	return true;
}

// Initializes new instance.
//	@mode: action mode (read, write, all)
void TabFile::Init(eAction mode)
{
	_fieldPos = NULL;
	_lineSpecLen = 0;
	if( mode != WRITE && IsGood() )	{
		_fieldPos = new short[_params.MaxFieldCnt];
		if(_params.LineSpec)	_lineSpecLen = strlen(_params.LineSpec);
	}
}

// Reads first line and set it as current.
// Throws exception if invalid and aborting file is set
//	@cntLines: returned estimated count of lines.
//	It works properly only if lines are sorted by ascending, f.e. in sorted bed-files.
//	return: current line
const char* TabFile::GetFirstLine(ULONG *cntLines)
{
	*cntLines = 0;
	if( GetLine() )
		// Add 1 for case of absence EOL marker in the last line.
		*cntLines = (ULONG)(Length() / RecordLength() + 1);	
	return _currLine;
}

const char*	TabFile::GetLine()
{
	char*	currLine;	// a bit faster than using _currLine in heep

	// fill _fieldPos 0 to check if all fields will be initialize
	// _fieldPos keeps start positions (next after TAB) for each field in line.
	// First field's position usually is 0, but if there are blanks in the beginning of the line,
	// it will be in first unblank position.
	// GetRecord(..) fills _fieldPos beginning from second field.
	memset(_fieldPos, vUNDEF, sizeof(short) * _params.MaxFieldCnt);

	currLine = GetRecord(NULL, _fieldPos, _params.MaxFieldCnt);
	if( currLine != NULL ) {
		USHORT	currPos = 0;	// a bit faster than using _currPos in heep
	
		// skip blanks at the beginning of line
		while( *(currLine+currPos)==BLANK )	currPos++;
		// skip comment line
		if( *(currLine+currPos)==_params.Comment )
			return GetLine();
		// check line for specifier
		if(_params.LineSpec && memcmp(currLine+currPos, _params.LineSpec, _lineSpecLen) )
			return GetLine();

		_fieldPos[0] = currPos;		// set start position of first field
		BYTE i = 1;
		// check obligatory field's positions, excluding first & last one
		for(; i<_params.MinFieldCnt-1; i++) {
			if( _fieldPos[i] == vUNDEF ) {		// numbers of TABs is less than number of required fields
				SetError(Err::TF_FIELD);
				return _currLine = NULL;
			}
			currLine[_fieldPos[i]-1] = '\0';	// replace TABs by 0
		}
		
		for(; i<_params.MaxFieldCnt; i++) {
			if( _fieldPos[i] == vUNDEF )		// numbers of TABs is less than number of optional fields
				break;
			currLine[_fieldPos[i]-1] = '\0';	// replace TABs by 0
		}
		
		currLine[RecordLength()-1] = '\0';	// close last position by 0

		// this version when _fieldPos is not initialized in TxtFile::GetRecord(): a bit more slower
		//for(BYTE i=1; i<_cntFields; i++) {
		//	// search for the next TAB
		//	for(; currPos<_recLen; currPos++)
		//		if(currLine[currPos]==TAB)
		//			break;
		//	currLine[currPos] = '\0';		// replace TABs by 0
		//	_fieldPos[i] = ++currPos;
		//}
	}
	return _currLine = currLine;
}

/************************ end of class TabFile ************************/

#if !defined _WIGREG

/************************ class FaFile ************************/
const string FaFile::Ext = ".fa";
#define FA_COMMENT	'>'

// Sets the initial state.
void FaFile::Pocket::Clear()
{
	_currPos = _countN = _currGapStart = 0;
	_defRgns.Clear();
}

// Adds some amount of 'N' as a gap to assign defined regions
//	@gapStart: gap's start position.
//	@gapLen: gap's length; if 0, complete adding gaps
void FaFile::Pocket::AddN(chrlen shiftGapStart, chrlen gapLen)
{
	shiftGapStart += _currPos;
	_defRgns.AddGap(shiftGapStart, _currGapStart, _minGapLen);
	_currGapStart = shiftGapStart + gapLen;
	_countN += gapLen;
}

// Opens an existing .fa fil and reads first line.
//	@fName: full .fa file name
//	@pocket: external temporary variables
FaFile::FaFile(const string & fName, Pocket& pocket) : TxtInFile(fName, READ, 1)
{
	// set header length
	chrlen len = 0;
	if( GetLine(pocket)[0] == FA_COMMENT ) {	// header line?
		len = RecordLength();
		pocket.Clear();
		GetLine(pocket);				// first data line
	}
	// set genome length
	len = chrlen(Length() - len);
	pocket._cLen = len - EOLSize() * 
		(len/RecordLength() +			// amount of EOL markers for whole lines
		(len%RecordLength() ? 1 : 0) );	// EOL for part line
}

// Reads line and set it as current, and generates regions
// Since method scans the line, it takes no overhead expenses to do this.
//	@pocket: external temporary variables
//	return: current line.
const char* FaFile::GetLine(Pocket& pocket) 
{
	chrlen countN = 0;
	const char* line = GetRecord(&countN);
	if( line ) {
		if( countN )			// some nt in line are undefined
			for(chrlen i=0; i<LineLength(); i++)
				if( line[i]==cN ) {
					pocket.AddN(i, countN);
					break;
				}
		pocket._currPos += LineLength();
	}
	return line;
}

//#if defined _FILE_WRITE  && defined DEBUG
//FaFile::FaFile(const string & fName, const char *chrName) : TxtFile(fName, WRITE, 1)
//{
//	const char *title = ">Chrom::Abbr";
//	char *header = new char[strlen(title) + strlen(chrName) + 1];
//	strcpy(header, title);
//	strcpy(header+strlen(title), chrName);
//	RecordToIOBuff(header, strlen(title) + strlen(chrName));
//	delete [] header;
//}
//#endif

/************************ end of class FaFile ************************/

#endif	// _ISCHIP, _DENPRO, _BIOCC
#endif	// _FQSTATN

#ifdef _FQSTATN

/************************ FqFile ************************/

// Returns checked length of current readed Read.
readlen FqFile::ReadLength() const
{
	CheckGettingRecord();
	return LineLengthByInd(READ) - EOLSize(); 
}

// Gets checked Read from readed Sequence.
const char* FqFile::GetCurrRead() const
{ 
	CheckGettingRecord();
	return NextRecord() - RecordLength() + LineLengthByInd(HEADER1); 
}

// Returns checked Sequence
const char* FqFile::GetSequence()
{
	const char* record = GetRecord();
	if(record != NULL) {
		if(*record != AT)
			Err(Err::FQ_HEADER, LineNumbToStr(0).c_str()).Throw();
		if( *(record + LineLengthByInd(HEADER1, true) + LineLengthByInd(READ, true)) != PLUS )
			Err(_errCode = Err::FQ_HEADER2, LineNumbToStr(2).c_str()).Throw();
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
//	_linesInRecCnt(file._linesInRecCnt),
//	_buffLineLen(0)
//{
//	if( !SetBasic(fName, READ, NULL) )	return;
//	_buffLen = file._buffLen;
//	_buff = file._buff;
//	RaiseFlag(CONSTIT);	
//}
