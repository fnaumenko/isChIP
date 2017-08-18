#include "bioTxtFile.h"

/************************ class FqFile ************************/

#ifdef _FILE_WRITE

#define AT '@'
#define PLUS '+'

string FqFile::Ext = ".fq";	// Extention of output FQ file (.fq or .fq.gz)
rowlen FqFile::ReadStartPos = 0;	// Read field constant start position 

#endif	// _FILE_WRITE

#ifdef _FQSTATN

// Returns checked length of current readed Read.
readlen FqFile::ReadLength() const
{
	CheckGettingRecord(FileName());
	return LineLength(READ) - EOLSize(); 
}

// Gets checked Read from readed Sequence.
char* FqFile::GetCurrRead() const
{ 
	CheckGettingRecord(FileName());
	return NextRecord() - RecordLength() + LineLength(HEADER1); 
}

// Returns checked Sequence
char* FqFile::GetSequence()
{
	char* record = GetRecord();
	if( record != NULL ) {
		if( *record != AT )
			Err(_errCode = Err::FQ_HEADER, LineNumbToStr(0)).Throw();
		if( *(record + LineLength(HEADER1) + LineLength(READ)) != PLUS )
			Err(_errCode = Err::FQ_HEADER2, LineNumbToStr(2)).Throw();
	}
	return record;
}

#endif	// _FQSTATN

#ifdef _FILE_WRITE

// Presets line write buffer
void FqFile::InitBuffer()
{
	LineSetOffset(ReadStartPos + Read::Len + 1);
	LineAddChar(PLUS, true);
	LineFill(ReadStartPos + Read::Len + 3, Read::SeqQuality, Read::Len);	// quality line
}

// Initializes line write buffer
void FqFile::InitToWrite()
{
	if(!ReadStartPos)
		ReadStartPos = Read::OutNameLength + 2;		// Read name + AT + EOL
	
	SetWriteBuffer(ReadStartPos + 2*Read::Len + 3, EOL);
	InitBuffer();
}

// Forms Read from fragment and adds it to the file.
//	@rName: Read's name
//	@read: valid read
//	@reverse: if true then complement added read 
//	@mate: SINGLE for one-side sequensing, or MATE_FIRST | MATE_SECOND for paired-end
void FqFile::AddRead(const string& rName, const char* read, bool reverse, eMate mate)
{
	// save Read
	LineSetOffset(ReadStartPos);
	if(reverse)
		Read::CopyComplement(LineCurrPosBuf(), read);
	else
		LineCopyChars(read, Read::Len);
	LineAddStrBack(rName);			// Read name
	LineAddCharBack(AT);
	LineBackToBuffer();

	//if( IsBadReadPtr(_record + SeqHeaderLen, Read::Len) )	cout << "BAD PTR\n";
}

#endif	// _FILE_WRITE

/************************ end of class FqFile ************************/

/************************ class BedRFile ************************/

string BedRFile::Ext = ".bed";			// Extention of output BED file .bed or .bed.gz)
const string BedRFile::Score = "60";

// Initializes line write buffer; only for master, clones are initialized by master
inline void BedRFile::InitToWrite() {
	SetWriteBuffer( 
		Chrom::MaxAbbrNameLength + // length of chrom name
		Read::OutNameLength + 		// length of Read name
		2 * CHRLEN_CAPAC +			// start + stop positions
		Score.length() +			// score
		1 + 6 + 1,					// strand + 5 TABs + EOL  + one for safety
		TAB
	);
}

// Sets chrom's name to line write buffer
void BedRFile::BeginWriteChrom(chrid cID)
{
	LineSetOffset(0);
	_offset = LineAddStr(Chrom::AbbrName(cID));
	LineSetOffset(_offset);
}

// Adds Read to the line's write buffer.
//	@rName: Read's name
//	@pos: valid Read's start position
//	@reverse: if true then set reverse strand, otherwise set forward
void BedRFile::AddRead(const string& rName, chrlen pos, bool reverse)
{
	LineAddStr(NNSTR(pos, TAB, pos+Read::Len));	// start, end positions
	LineAddStr(rName);							// Read name
	LineAddStr(Score);							// score
	LineAddChar(Read::Strand[int(reverse)]);	// strand
	LineToBuffer(_offset);
}

// Adds two mate Reads to the line's write buffer.
//	@rName: Read's name
//	@pos1: valid first mate Read's start position
//	@pos2: valid second mate Read's start position
void BedRFile::AddTwoReads(const string& rName, chrlen pos1, chrlen pos2)
{
	AddRead(rName + Read::NmSuffMate1, pos1, false);
	AddRead(rName + Read::NmSuffMate2, pos2, true);
}

/************************ end of class BedRFile ************************/

/************************ class SamFile ************************/

/*
SAM format: field FLAG
Bit		Description
-------------------
1 0x1	template having multiple segments in sequencing
2 0x2	each segment properly aligned according to the aligner
3 0x4	segment unmapped
4 0x8	next segment in the template unmapped
5 0x10	SEQ being reverse complemented
6 0x20	SEQ of the next segment in the template being reversed
7 0x40	the rst segment in the template
8 0x80	the last segment in the template
9 0x100	secondary alignment
0 0x200	not passing quality controls
1 0x400	PCR or optical duplicate
2 0x800	supplementary alignment

 0	   00000	SE +
16	   10000	SE -
99	01100011	PE +
147	10010011	PE -
*/

//#define RNEXT_SE '*'	// Ref. name of the mate/next SE read
#define RNEXT_PE '='	// Ref. name of the mate/next PE read

string SamFile::Ext = ".sam";			// Extention of output SAM file (.sam or .sam.gz)
string SamFile::Comb5_6;				// combined value from 5 to 6 field: defined in constructor
const string SamFile::Comb7_9 = "*\t0\t0\t";	// combined value from 7 to 9 field for SE mode
string SamFile::Flag[2];
rowlen SamFile::Offset5_9 = 0;
rowlen SamFile::ReadStartPos =	0;

// Creates and initialize line write buffer.
void SamFile::InitBuffer()
{
	LineFill(0, TAB);
	LineFill(ReadStartPos + Read::Len + 1, Read::SeqQuality, Read::Len);	// Read quality
	if( !OutFile::PairedEnd() ) {						// preset 5-9 fields for SE mode
		Offset5_9 = Comb7_9.length() + Comb5_6.length() + 1;
		LineSetOffset(ReadStartPos - Offset5_9);
		LineAddStr(Comb5_6);
		LineAddStr(Comb7_9);
	}
}

// Writes header line to line write buffer.
//	@tag0: line tag
//	@tag1: first subtag
//	@val1: first value
//	@tag2: second subtag
//	@val2: second value
//	@closeLine: true if line is ended
void SamFile::SetHeaderLine(
	const char* tag0,
	const char* tag1, const string& val1,
	const char* tag2, const string& val2,
	bool closeLine
	)
{
	ostringstream oss;

	if( tag0 ) {
		LineAddChar(AT);
		LineAddStr(tag0);
	}
	oss << tag1 << COLON << val1;
	LineAddStr(oss.str());
	if( tag2 ) {
		oss.rdbuf()->str(sBLANK);		// clear too long stream by any one-letter string
		oss.seekp(0);
		oss << tag2 << COLON << val2;
		LineAddStr(oss.str());
	}
	if(closeLine)	
		LineDecreaseOffset();
	LineToBuffer(0, closeLine);
	_headLineCnt++;
}

// Generates and writes SAM header
//	@cFiles: chrom files
//	@commandLine: command line
void SamFile::CreateHeader(const ChromFiles& cFiles, const string& commandLine)
{
	const string& path = cFiles.Path();
	vector<string> szFiles;
	FS::GetFiles(szFiles, path, ChromSizes::Ext);
	if(szFiles.size() > 1)
		Err("more than one chromosome sizes files", path).Throw();
	ChromSizes* cSizes;
	if(szFiles.size())	// chrom.sizes exists
		cSizes = new ChromSizes(path + szFiles[0]);
	else {				// chrom.sizes doesn't exist
		Timer tm;
		cout << "Chromosome sizes file" << MsgFileAbsent;
		fflush(stdout);
		tm.Start();
		cSizes = new ChromSizes(cFiles);
		cout << MsgDone;
		tm.Stop();
	}

	SetHeaderLine("HD", "VN", "1.0", "SO", "unsorted");
	for(ChromFiles::cIter it=cFiles.cBegin(); it!=cFiles.cEnd(); it++)
		if( cFiles.IsTreated(it) )
			SetHeaderLine("SQ", "SN", Chrom::AbbrName(CID(it)), "LN", NSTR((*cSizes)[CID(it)]));
	SetHeaderLine("PG", "ID", Product::Title, "PN", Product::Title, false);
	SetHeaderLine(NULL, "VN", Product::Version, NULL, StrEmpty, false);
	SetHeaderLine(NULL, "CL", commandLine);
	Write();	// write header in case of multithread
	delete cSizes;
}

#define CL_LEN 3	// length of the header line tag 'CL:'

// Initializes line write buffer and header and makes ready for writing
//	@cFiles: chrom files
//	@commandLine: command line
void SamFile::InitToWrite(const ChromFiles& cFiles, const string& commandLine)
{
	if(!ReadStartPos) {		// first constructor call?
		if( OutFile::PairedEnd() )	{ Flag[0] = "99";	Flag[1] = "147"; }
		else						{ Flag[0] = "0";	Flag[1] = "16";	 }
		Comb5_6 = NNSTR(Read::MapQuality, TAB, Read::Len) + "M"; // MAPping Quality + CIGAR: Read length

		ReadStartPos =				// maximal length of write line buffer
									// without Read & Quality fields, with delimiters
		Read::OutNameLength +		// QNAME: Read name
		3 +							// FLAG: bitwise FLAG
		Chrom::MaxAbbrNameLength +	// RNAME: AbbrChromName
		CHRLEN_CAPAC +				// POS: 1-based start pos
		2 +							// MAPQ: MAPping Quality
		4 +							// CIGAR: Read length in 3 digits + letter 'M' or '='
		1 +							// RNEXT: Ref. name of the mate/next read: SE '*', PE '='
		CHRLEN_CAPAC +				// PNEXT: Position of the mate/next read
		3 + 						// TLEN: observed Template LENgth
		10 + 1 + 1;					// number of TABs + EOL  + one for safety

		rowlen buffLen = max(
			rowlen(ReadStartPos+2*Read::Len), rowlen(commandLine.length()+CL_LEN));
		SetWriteBuffer(buffLen, TAB);
		ReadStartPos = buffLen - 2*Read::Len - 1;
		CreateHeader(cFiles, commandLine);
	}
	InitBuffer();
}

// Adds full-defined Read to the line's write buffer.
//	@rName: Read's name
//	@read: valid Read
//	@flag: FLAG field value
//	@pos1: valid start position of mate Read (PE) or Read (SE)
//	@pos2: valid start position of mate Read (PE) or -1 (SE)
//	@fLen: +/- fragment's length (PE) or 0 (SE)
void SamFile::AddStrongRead(const string& rName, const char* read, const string& flag,
	chrlen pos1, chrlen pos2, int fLen)
{
	ostringstream oss;

	LineSetOffset(ReadStartPos);
	LineCopyChars(read, Read::Len);	// SEQ: Read
	if(fLen) {						// PE
		oss	<< RNEXT_PE << TAB		// RNEXT: Ref. name of the mate/next PE read
			<< ++pos2 << TAB		// PNEXT: 1-based second mate Read's start position
			<< fLen;				// TLEN: observed Template LENgth
		LineAddStrBack(oss.str());
		oss.rdbuf()->str(sBLANK);	// clear too long stream by any one-letter string
		oss.seekp(0);
		LineAddStrBack(Comb5_6);	// MAPQ + CIGAR
	}
	else							// SE
		LineSetOffset(ReadStartPos - Offset5_9);
	oss << ++pos1;	
	LineAddStrBack(oss.str());		// POS: 1-based start position
	LineAddStrBack(_cName);			// RNAME: chrom's name
	LineAddStrBack(flag);			// FLAG
	LineAddStrBack(rName);			// QNAME: Read name

	LineBackToBuffer();
}

// Adds Read to the line's write buffer.
//	@rName: Read's name
//	@read: valid Read
//	@pos: valid Read's start position
//	@reverse: if true then set reverse strand, otherwise set forward
inline void SamFile::AddRead(const string& rName, const char* read, chrlen pos, bool reverse)
{
	AddStrongRead(rName, read, Flag[int(reverse)], pos);
}

// Adds two mate Reads to the line's write buffer.
//	@rName: name of Read
//	@read1: valid first mate Read
//	@read2: valid second mate Read
//	@pos1: valid first mate Read's start position
//	@pos2: valid second mate Read's start position
//	@fLen: fragment's length
void SamFile::AddTwoReads(const string& rName,
	const char* read1, const char* read2, chrlen pos1, chrlen pos2, int fLen)
{
	AddStrongRead(rName, read1, Flag[0], pos1, pos2, fLen);
	AddStrongRead(rName, read2, Flag[1], pos2, pos1, -fLen);
}

/************************ end of class SamFile ************************/

/************************ class OutFile ************************/

BYTE OutFile::Mode = 0;	// 0 in case of one-side sequencing, 1 in case paired-end

OutFile::AddReads OutFile::callAddRead[] =
	{ &OutFile::AddReadSE, &OutFile::AddReadPE, &OutFile::NoAddRead };

// Creates new instance for writing.
//	@fName: common file name without extention
//	@mask: combination of 0x1, 0x2, 0x4
//	@pairedEnd: 0 in case of one-side sequencing, 1 in case of paired-end
//	@isZipped: true if output files should be zipped
OutFile::OutFile(const string& fName, int mask, UINT pairedEnd,	bool isZipped)
{
	const string& zipExt = isZipped ? ZipFileExt : StrEmpty;
	_mode = Mode = BYTE(pairedEnd);

	_fqFile1 = _fqFile2 = NULL;
	if( mask & 0x1 )
		if( PairedEnd() ) {
			_fqFile1 = new FqFile(fName, zipExt, 1);
			_fqFile2 = new FqFile(fName, zipExt, 2);
		}
		else
			_fqFile1 = new FqFile(fName, zipExt);
	_bedFile = mask & 0x2 ? new BedRFile(fName, zipExt) : NULL;
	_samFile = mask & 0x4 ? new SamFile	(fName, zipExt) : NULL;
}

#ifdef _MULTITHREAD
// Creates a clone of existed instance for writing.
//	@file: original instance
//	@threadNumb: number of thread
OutFile::OutFile(const OutFile& file, threadnumb threadNumb)
{
	_mode = Mode;
	_fqFile1 = file._fqFile1 ?	new FqFile	(*file._fqFile1, threadNumb) : NULL;
	_fqFile2 = file._fqFile2 ?	new FqFile	(*file._fqFile2, threadNumb) : NULL;
	_bedFile = file._bedFile ?	new BedRFile(*file._bedFile, threadNumb) : NULL;
	_samFile = file._samFile ?	new SamFile	(*file._samFile, threadNumb) : NULL;
}
#endif

OutFile::~OutFile()
{
	if(_fqFile1)	delete _fqFile1;
	if(_fqFile2)	delete _fqFile2;
	if(_bedFile)	delete _bedFile;
	if(_samFile)	delete _samFile;
}

// Initializes buffers and makes ready for writing
//	@cFiles: chrom files
//	@commandLine: command line
void OutFile::Init(const ChromFiles& cFiles, const string& commandLine)
{
	if(_fqFile1)	_fqFile1->InitToWrite();
	if(_fqFile2)	_fqFile2->InitToWrite();
	if(_bedFile)	_bedFile->InitToWrite();
	if(_samFile)	_samFile->InitToWrite(cFiles, commandLine);
}

ULONG OutFile::Count() const
{
	if(_fqFile1)	return _fqFile1->RecordCount() << PairedEnd();
	if(_bedFile)	return _bedFile->RecordCount();
	if(_samFile)	return _samFile->Count();
	return 0;
}

// Adds one SE Read
int OutFile::AddReadSE(string& rName, const Nts& nts,
	ULONG rNumb, chrlen pos, fraglen fragLen, bool reverse)
{
	int ret;
	if(reverse)	pos += (fragLen - Read::Len);
	const char* read = nts.Read(pos);
	if( (ret = Read::CheckNLimit(read)) <= 0 )	return ret;

	rName += NSTR(rNumb ? rNumb : pos);
	if(_fqFile1)	_fqFile1->AddRead(rName, read, reverse, TxtFile::MATE_SINGLE);
	if(_bedFile)	_bedFile->AddRead(rName, pos,  reverse);
	if(_samFile)	_samFile->AddRead(rName, read, pos, reverse);
	return 1;
}

// Adds two PE Reads
//	@reverse: not used
int OutFile::AddReadPE(string& rName, const Nts& nts,
	ULONG rNumb, chrlen pos, fraglen fragLen, bool reverse)
{
	int ret;
	chrlen pos2 = pos + fragLen - Read::Len;
	const char* read2 = nts.Read(pos2);
	if( (ret = Read::CheckNLimit(read2)) <= 0 )	return ret;
	const char* read1 = nts.Read(pos);
	if( (ret = Read::CheckNLimit(read1)) <= 0 )	return ret;

	rName += (rNumb ? NSTR(rNumb) : NNSTR( pos, Read::NmPosDelimiter, pos2 ));
	if(_fqFile1) {
		_fqFile1->AddRead(rName + Read::NmSuffMate1, read1, false, TxtFile::MATE_FIRST);
		_fqFile2->AddRead(rName + Read::NmSuffMate2, read2, true, TxtFile::MATE_SECOND);
	}
	if(_bedFile)	_bedFile->AddTwoReads(rName, pos, pos2);
	if(_samFile)	_samFile->AddTwoReads(rName, read1, read2, pos, pos2, fragLen);
	return 1;
}

void OutFile::Write() const
{
	if(_fqFile1)
		if( PairedEnd() ) {
			_fqFile1->Write(TxtFile::MATE_FIRST);
			_fqFile2->Write(TxtFile::MATE_SECOND); 
		}
		else
			_fqFile1->Write(); 
	if(_bedFile)	_bedFile->Write();
	if(_samFile)	_samFile->Write();
}

// Prints output files as head info
//	@signOut: output marker
void OutFile::Print(const char* signOut) const
{
	if(_fqFile1) {
		cout << signOut << "Output sequence: " << _fqFile1->FileName();
		if(_fqFile2)
			cout << ", " << _fqFile2->FileName(); 
		cout << endl;
	}
	if(_bedFile || _samFile) {
		cout << signOut << "Output alignment: ";
		if(_bedFile)	cout << _bedFile->FileName();
		cout << " ";
		if(_samFile)	cout << _samFile->FileName();
		cout << endl;
	}
	cout << signOut << "Sequencing: " << (PairedEnd() ? "paired" : "single") << "-end\n";
}

/************************ end of class OutFile ************************/

//#define MIN_NUMB	100

SamFragDistr::SamFragDistr(const string& fname)
{
	//ULONG cntLines;
	TabFile file(fname, 9, true, TxtFile::READ, NULL, '@', true);
	 
	cout << "Fragments distribution in " << fname << endl;

	int numbLen = 1000;
	Array<int> numbers(numbLen);
	int v, i, insCnt = 0;
	bool checkPair = false;

	const char* currLine;	// = file.GetFirstLine(&cntLines);
	for( i=1; currLine=file.GetLine(); i++ )
		if( i%2 && (v=file.IntField(8)) ) {
			if(v<0)	v *= -1;
			if(v < numbLen) {	numbers[v]++;	insCnt++; }
		}
	// print
	//cout << EOL;
	for(int k=0; k<numbLen; k++)
		if(numbers[k])
			cout << k << TAB << numbers[k] << EOL;
	cout << "counted " << insCnt << " from " << i << BLANK << sPercent(insCnt, i, 3, 0, true);
	cout << endl;
}
