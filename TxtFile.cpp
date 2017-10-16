#include "TxtFile.h"

/************************ class TxtFile ************************/
const char* modes[] = { "r", "w", "a+" };
const char* bmodes[] = { "rb", "wb" };

// Sets error code and throws exception if it is allowed.
void TxtFile::SetError(Err::eCode errCode) const
{
	_errCode = errCode;
	if( IsFlagSet(ABORTING) )
		Err(errCode, FileNameToExcept().c_str()).Throw();
}

// Initializes instance variables and opens a file with setting a proper error code in case of fault.
//	@fName: valid full name of file
//	@mode: opening mode
//	@file: clonable file or NULL
//	return: true if success, otherwise false.
bool TxtFile::SetBasic(const string& fName, eAction mode, void* file)
{
	_buff = _buffLine = NULL;
	_linesLen = NULL;
	_stream = NULL;
	_errCode = Err::NONE;
	_fName = fName;
	_currRecPos = _recLen = _cntRecords = _readingLen = 0;
#ifdef _NO_ZLIB
	if( IsZipped() ) { SetError(Err::FZ_BUILD); return false; }
#endif
	// set file stream
	if( file )	_stream = file;
	else {
#ifndef _NO_ZLIB
		if(IsZipped())
			if( mode == ALL )	
				SetError(Err::FZ_OPEN);
			else				
				_stream = gzopen(fName.c_str(), bmodes[mode]);
		else
#endif
			_stream = fopen(fName.c_str(), modes[mode]);
		if( _stream == NULL )	SetError(Err::F_OPEN);
	}
	return IsGood();
}

// Allocates memory for file read/write buffer or write line buffer with checking.
//  return: true if successful
bool TxtFile::CreateBuffer(eBuff buffType)
{
	try {
		if( buffType == BUFF_BASIC )	_buff = new char[_buffLen];
		else if(!_buffLine)	{
			_buffLine = new char[_buffLineLen];
			memset(_buffLine, _delim, _buffLineLen);
			_buffLineOffset = 0;
		}
	}
	catch(const bad_alloc)	{ SetError(Err::F_MEM); };
	return IsGood();
}

// General constructor, always is called first.
//	@fName: valid full name of file
//	@mode: opening mode
//	@cntRecLines: number of lines in a record
//	@abortInvalid: true if invalid instance shold be completed by throwing exception
//	@rintName: true if file name should be printed in exception's message
TxtFile::TxtFile (const string& fName, eAction mode, BYTE cntRecLines, bool abortInvalid, bool printName) :
	_flag(0),
	_cntRecLines(cntRecLines),
	_buffLineLen(0)
{
	SetFlag(ZIPPED, FS::HasGzipExt(fName));
	SetFlag(ABORTING, abortInvalid);
	SetFlag(PRNAME, printName);
	if( !SetBasic(fName, mode, NULL) )	return;
	// set file's and buffer's sizes
	_buffLen = NUMB_BLK * BASE_BLK_SIZE;
	_fSize = FS::Size(fName.c_str());
	if( _fSize == -1 )	_fSize = 0;		// new file
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
		_buffLen >>= 1;		// decrease block size twice because of allocating additional memory:
							// 2*_buffLen for writing or 3*_buffLen for reading by gzip
	}
#endif
	if( _fSize && _fSize < _buffLen )
		if( mode != WRITE )	
			_buffLen = (ULONG)_fSize + 1;	// for reading
		else if( _fSize * 2 < _buffLen )
			_buffLen = (ULONG)_fSize * 2;	// for writing: increase small buffer for any case
	
	if( !CreateBuffer(BUFF_BASIC) )	return;

#ifdef ZLIB_NEW
	if( IsZipped() && gzbuffer( (gzFile)_stream, _buffLen) == -1 )
	{ SetError(Err::FZ_MEM); return; }
#endif

	if(mode != WRITE) {
		if(ReadBlock(0) < 0)	return;		// read first block
		_linesLen = new UINT[cntRecLines];	// set lines buffer
	}
	else {
		memset(_buff, 0, _buffLen);
		_linesLen = NULL;
	}
	//cout << "constructor " << _fName << " done" << endl;
}

#ifdef _MULTITHREAD
// Creates new instance with basic buffer belonges to aggregated file: constructor for concatenating.
//	For reading only.
//	@fName: valid full name of file
//	@file: file that is aggregated
TxtFile::TxtFile(const string & fName, const TxtFile& file) :
	_flag(file._flag),
	_cntRecLines(file._cntRecLines),
	_buffLineLen(0)
{
	if( !SetBasic(fName, READ, NULL) )	return;
	_buffLen = file._buffLen;
	_buff = file._buff;
	RaiseFlag(CONSTIT);	
}

// Creates a clone of existing instance.
// Clone is a copy of opened file with its own buffer for writing only.
//	Used for multithreading file recording
//	@file: opened file which is cloned
//	@threadNumb: number of thread
//	Exception: file_error
TxtFile::TxtFile(const TxtFile& file, threadnumb threadNumb) :
	_flag(file._flag),
	_cntRecLines(file._cntRecLines),
	_buffLineLen(file._buffLineLen),
	_delim(file._delim)
{
	if( !SetBasic("Clone " + file._fName, WRITE, file._stream) )	return;
	// block size is slightly increased by increasing thread number
	// for mistiming block's writing to file as possible
	_buffLen = BASE_BLK_SIZE * (threadNumb + NUMB_BLK);
	RaiseFlag(CLONE);
	CreateBuffer(BUFF_BASIC);
	CreateBuffer(BUFF_LINE);
}
#endif

TxtFile::~TxtFile()
{
	if( _linesLen )						delete [] _linesLen;
	if( _buff && !IsFlagSet(CONSTIT) )	delete [] _buff;
	if( _buffLine )	{
		//cout << "delete _buffLine\n";
		delete [] _buffLine;
	}
	if( _stream && !IsClone() )	{
		int res = 
#ifndef _NO_ZLIB
			IsZipped() ?
				gzclose( (gzFile)_stream) :
#endif
				fclose( (FILE*)_stream);
		if( res )	SetError(Err::F_CLOSE);
	}
}

// Reads next block.
//	@offset: shift of start reading position
//	return: 1 if file is not finished; 0 if it is finished; -1 if unsuccess reading
int TxtFile::ReadBlock(const UINT offset)
{
	size_t readLen;
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
	
	_readingLen = readLen + offset;
//#ifdef ZLIB_OLD
	//if( _readTotal + _readingLen > _fSize )
	//	_readingLen = _fSize - _readTotal;
	//_readTotal += _readingLen;
//#endif
	if( _readingLen != _buffLen && !IsZipped() && ferror((FILE*)_stream) )
	{ SetError(Err::F_READ); return -1; }
	_currRecPos = 0;
	return _readingLen == 0 ? 0 : 1;
}

// Sets _currLinePos to the beginning of next non-empty line in _buff
//	@counterN: if not NULL, adds to counterN the number of 'N' in a record. Used in Fa() only.
//	@posTab: if not NULL, sets TABs positions in line to this array
//	@cntTabs: if @posTab is not NULL, the length of @posTab array
//	return: point to line or NULL if no more lines
char* TxtFile::GetRecord(chrlen* const counterN, short* const posTab, BYTE cntTabs)
{
	if( _fSize == 0 )	return ReadingEnded();
	UINT i, 
		 currLinePos = _currRecPos,	// start position of current readed record
		 cntEmpty = 0;				// counter of empty lines
	chrlen cntN = 0;				// local counter of 'N'
	BYTE indTab = 1;				// index of TAB position in @posTab
	char c;

	_recLen = 0;
	for(BYTE r=0; r<_cntRecLines; r++)
		for(i=currLinePos;; i++)	{
			if( i >= _readingLen ) {			// check for the next block
				if( _readingLen != _buffLen )	// final block
					if( i == _buffLen			// current line is not ended by EOL
					|| _buff[i-1] == EOL		// current line which is ended by EOL
					|| _buff[i-1] == '\0'	)	// EOL in the current line has been replaced by 0 in previous call TabFile::GetLine()
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
				UINT restLen = _readingLen - _currRecPos - cntEmpty;
				// move untreated rest to the beginning of buffer; if restLen=0, skip moving
				memmove(_buff, _buff + _currRecPos, restLen);
				if( !ReadBlock(restLen) )
					return ReadingEnded();	// record was finished exactly by previous ReadBlock()
				indTab = 1;
				cntN = currLinePos = cntEmpty = i = _recLen = r = 0;
			}
			c = _buff[i];
			if( counterN && c==cN )
				cntN++;
			else if( indTab < cntTabs && c==TAB )
				posTab[indTab++] = i - currLinePos + 1;
			else if( c==EOL ) {
				if( i==currLinePos ) {			// EOL marker is first in line
					currLinePos++; cntEmpty++;	// skip empty line
					continue;
				}
				if( !IsFlagSet(EOLSZSET) ) {
					SetFlag(EOLSZ, _buff[i-1]==CR);
					RaiseFlag(EOLSZSET);
				}
A:				_recLen += (_linesLen[r] = UINT(i-currLinePos) + 1);
				currLinePos = i+1;
				break;
			}
		}
	_currRecPos = currLinePos;			// next record position
	 if(counterN)		(*counterN) += cntN;
	 _cntRecords++;
	return _buff + _currRecPos - _recLen;
}

#ifdef _FILE_WRITE

//	Adds delimiter on the given shift and increases current position.
//	return: new current position
rowlen TxtFile::LineAddDelim()
{
	_buffLine[_buffLineOffset] = _delim;
	return ++_buffLineOffset;
}

// Copies block of chars to the current position in the line write buffer,
//	adds delimiter after string	and increases current position.
//	@src: pointer to the block of chars
//	@num: number of chars
//	return: new current position
rowlen TxtFile::LineAddChars(const char* src, size_t num)
{
	memcpy(_buffLine+_buffLineOffset, src, num);
	_buffLineOffset+=num;
	return LineAddDelim();
}

// Copies block of chars before the current position in the line write buffer,
//	adds delimiter before string and decreases current position.
//	@src:  pointer to the block of chars
//	@num: number of chars
void TxtFile::LineAddCharsBack(const char* src, size_t num)
{
	_buffLine[--_buffLineOffset] = _delim;
	_buffLineOffset -= num;
	memcpy(_buffLine + _buffLineOffset, src, num);
}

// Initializes line write buffer.
//	@len: length of buffer
//	@delim: delimiter between line fields
void TxtFile::SetWriteBuffer(rowlen len, char delim)
{
	_delim = delim;
	_buffLineLen = len;
	CreateBuffer(BUFF_LINE);
}

// Sets some bytes in the line write buffer to the specified value.
//	@offset: shift of buffer start position of filling bytes
//	@val: value to be set
//	@len: number of bytes to be set to the value, or the rest of buffer by default
void TxtFile::LineFill(rowlen offset, char val, rowlen len)
{
	memset(_buffLine + offset, val, len ? len : (_buffLineLen-offset));
}

// Adds byte to the current position in the line write buffer,
//	adds delimiter after byte and increases current position.
//	@chr: value to be set
//	@addDelim: if true then adds delimiter and increases current position
void TxtFile::LineAddChar(char chr, bool addDelim)
{ 
	_buffLine[_buffLineOffset++] = chr;
	if(addDelim)	_buffLine[_buffLineOffset++] = _delim;
}
	
// Adds byte before the current position of the line write buffer,
//	adds delimiter before byte and decreases current position.
//	@chr: value to be set
void TxtFile::LineAddCharBack(char chr, bool addDelim)
{ 
	_buffLine[--_buffLineOffset] = chr;
	if(addDelim)	_buffLine[--_buffLineOffset] = _delim;
}

// Adds floating point value to the current position of the line write buffer,
//	adds delimiter after value and increases current position.
//	@val: value to be set
//	@ndigit: number of digits to generate
//	@addDelim: if true then adds delimiter and increases current position
void TxtFile::LineAddFloat(double val, BYTE ndigit, bool addDelim)
{
	// double val, because casting int to float is not precise, f.e. (float)61342430 = 61342432
	_gcvt(val, ndigit, _buffLine+_buffLineOffset);
	_buffLineOffset+=strlen(_buffLine+_buffLineOffset);
	if(addDelim)
		LineAddDelim();
}

// Adds first part of the line write buffer (from 0 to the current position)
//	to the file write buffer.
//	@offset: start shift for the next writing cycle
//	@closeLine: if true then close line by EOL
void TxtFile::LineToBuffer(rowlen offset, bool closeLine)
{
	AddRecord(_buffLine, _buffLineOffset, closeLine);
	_buffLineOffset = offset;
}

// Adds record to the file write buffer.
// Generate exception if writing is fall.
//	@src: record
//	@len: length of record
//	@closeLine: if true then close line by EOL
//	@mate: first, second file-mate or single file
void TxtFile::AddRecord(const char *src, UINT len, bool closeLine, eMate mate)
{
	if( _currRecPos + len + 1 > _buffLen )	// write buffer to file if it's full
		Write(mate);
	memcpy( _buff+_currRecPos, src, len );
	_currRecPos += len;
	if( closeLine )	_buff[_currRecPos++] = EOL;
	_cntRecords++;
}

// Writes current block to file.
//	@mate: SINGLE for single file or MATE_FIRST | MATE_SECOND for pair of files.
// Set up mutex for writing to synchronous files while multithreading.
void TxtFile::Write(eMate mate) const
{
	int res;
#ifdef _MULTITHREAD
	//if( mate != MATE_SECOND )	// rouse mutex for first mate and always by single file
	if( mate != MATE_SINGLE )	// rouse mutex for first mate and always by single file
		Mutex::Lock(Mutex::WR_FILE);
#endif
	res = 
#ifndef _NO_ZLIB
		IsZipped() ?
		gzwrite((gzFile)_stream, _buff, _currRecPos) :
#endif
		fwrite(_buff, 1, _currRecPos, (FILE*)_stream);
#ifdef _MULTITHREAD
	//if( mate != MATE_FIRST )	// drop mutex by second mate and always by single file
	if( mate != MATE_SINGLE )	// drop mutex by second mate and always by single file
		Mutex::Unlock(Mutex::WR_FILE);
#endif
	if( res != _currRecPos )	SetError(Err::F_WRITE);
	else						_currRecPos = 0;
}

//bool	TxtFile::AddFile(const string fName)
//{
//	TxtFile file(fName, *this);
//	while( file.ReadBlock(0) ) {
//		_currRecPos = file._readingLen;
//		if( !Write() )	return false;
//	}
//	return true;
//}

#endif	// _FILE_WRITE

/************************ end of class TxtFile ************************/

#ifdef _FILE_WRITE

/************************ class LineFile ************************/

// Adds to line int and float values and writes line.
void LineFile::WriteLine(int val1, float val2, BYTE ndigit)
{
	LineAddInt(val1);
	LineAddFloat(val2, ndigit);
	LineToBuffer();
}

// Adds to line 2 int values and writes line.
void LineFile::WriteLine(int val1, int val2)
{
	LineAddStr(NNSTR(val1, _delim, val2));
	LineToBuffer();
}

// Adds to line string and int value and writes line.
void LineFile::WriteLine(const string& str, int val)
{
	LineAddStr(str);
	LineAddInt(val, false);
	LineToBuffer();
}

/************************ end of class LineFile ************************/

#endif	// _FILE_WRITE

#ifndef _FQSTATN

#ifndef _WIGREG

/************************ struct Region ************************/

// Extends Region with chrom length control.
// If extended Region starts from negative, or ends after chrom length, it is fitted.
//	@extLen: extension length in both directions
//	@cLen: chrom length
void Region::Extend(chrlen extLen, chrlen cLen)
{
	Start -= extLen > Start ? Start : extLen;
	if(cLen && (End += extLen) > cLen)	End = cLen;
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
	LineFile file(fName, TAB);

	file.BeginWrite(2*INT_CAPACITY+2);
	file.WriteLine("minGap", minGapLen);
	for(Iter it = Begin(); it != End(); it++)
		file.WriteLine(it->Start, it->End);
	file.Write();
}

#endif	// _DENPRO, _BIOCC

// Adds gap, beginning of current gap position.
//	@gapStart: gap's start position.
//	@currGapStart: current gap start position holder
//	@minGapLen: minimal length which defines gap as a real gap
void Regions::AddGap(chrlen gapStart, chrlen currGapStart, chrlen minGapLen)
{
	if( gapStart && gapStart != currGapStart )		// current N-region is closed
		if( currGapStart-_regions.back().End-1 < minGapLen )
			_regions.back().End = gapStart-1;		// pass minimal allowed undefined nt
		else 
			_regions.push_back(Region(currGapStart, gapStart-1));	// add new def-region
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
const char TabFile::Comment = '#';

// Checks if field valid and throws exception if not.
//	@ind: field index
//	return: true if field is valid
bool TabFile::IsFieldValid	(BYTE ind) const
{
	if(_fieldPos[ind] == vUNDEF		// vUNDEF was set in buff if there was no such field in the line
	|| !SField(ind)[0]) {			// empty string was set if the line was ended by TAB (empty field)
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
		*cntLines = (ULONG)(Length() / RecordLength() + 1);	// Add 1 for case no EOL marker in the last line.
	// check _fieldPos
	//for(BYTE i=0; i<_cntFields; i++)
	//	cout << _fieldPos[i] << EOL;
	//if(!*cntLines)	SetError(Err::EMPTY);
	return _currLine;
}

const char*	TabFile::GetLine()
{
	USHORT	currPos;	// a bit faster than using _currPos in heep
	char*	currLine;	// a bit faster than using _currLine in heep

	// fill _fieldPos 0 to check if all fields will be initialize
	// _fieldPos keeps start positions for each field in line.
	// First field's position usually is 0, but in case of blanks at the beginning of line
	// it will be in first unblank position.
	// So GetRecord(..) fiils _fieldPos beginning from second field.
	memset(_fieldPos, vUNDEF, sizeof(short)*_params.MaxFieldCnt);

	currLine = GetRecord(NULL, _fieldPos, _params.MaxFieldCnt);
	if( currLine != NULL ) {
		// skip blanks on the beginning of line
		for(currPos=0; *(currLine+currPos)==BLANK; currPos++);
		// skip comment line
		if( *(currLine+currPos) == _params.Comment )
			return GetLine();
		// check line for specifier
		if(_params.LineSpec && memcmp(currLine+currPos, _params.LineSpec, _lineSpecLen) ) {
			//SetError(Err::TF_SPEC);
			//return _currLine = NULL;
			return GetLine();
		}

		// check readed field's positions, excluding first & last one
		if( _checkFieldCnt )
			for(BYTE i=1; i<_params.MaxFieldCnt-1; i++)
				if( _fieldPos[i] == vUNDEF ) {
					SetError(Err::TF_FIELD);
					return _currLine = NULL;
				}

		_fieldPos[0] = currPos;		// set start position of first field
		// replace TABs by 0
		for(BYTE i=1; i<_params.MaxFieldCnt; i++) {
			//if( (currPos=_fieldPos[i]) == vUNDEF ) {	// last position may be 0 if numbers of TABs in line is less than number of fields
			if( _fieldPos[i] == vUNDEF ) {
				currLine[RecordLength() -1] = '\0';		// close record by 0: it is necessery for last position
				break;
			}
			//currPos=_fieldPos[i];
			//currLine[currPos-1] = '\0';
			currLine[_fieldPos[i]-1] = '\0';
		}
		// this variant when _fieldPos is not initialized in GetRecord(): a bit more slower
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
FaFile::FaFile(const string & fName, Pocket& pocket) : TxtFile(fName, READ, 1)
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

#if defined _FILE_WRITE  && defined DEBUG
FaFile::FaFile(const string & fName, const char *chrName) : TxtFile(fName, WRITE, 1)
{
	const char *title = ">Chrom::Abbr";
	char *header = new char[strlen(title) + strlen(chrName) + 1];
	strcpy(header, title);
	strcpy(header+strlen(title), chrName);
	AddRecord(header, strlen(title) + strlen(chrName));
	delete [] header;
}
#endif

/************************ end of class FaFile ************************/

#endif	// _ISCHIP, _DENPRO, _BIOCC
#endif	// _FQSTATN
