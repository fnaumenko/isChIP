#include "OutTxtFile.h"

/************************ class BedRFile ************************/

// Initializes line write buffer; only for master, clones are initialized by master
//	@commandLine: command line to add as comment on the first line
void BedRFile::InitToWrite(const string& commandLine) {
	rowlen buffLen =
		Chrom::MaxAbbrNameLength +		// length of chrom name
		Read::OutNameLength + 			// length of Read name
		2 * CHRLEN_CAPAC +				// start + stop positions
		OutFile::MapQual.length() +		// score
		1 + 6;							// strand + 5 TABs + EOL

	ostringstream oss;
	oss << HASH << BLANK << commandLine;
	string firstLine = oss.str();
	if(firstLine.length() > int(buffLen))
		buffLen = firstLine.length();
	SetWriteBuffer(buffLen + 1, TAB);	// + 1 for safety
	AddRecord(firstLine.c_str(), commandLine.length());
	Write();	// write header in case of multithread
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
	LineAddStr(OutFile::MapQual);				// score
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

string SamFile::Comb5_6;				// combined value from 5 to 6 field: defined in constructor
const string SamFile::Comb7_9 = "*\t0\t0\t";	// combined value from 7 to 9 field for SE mode
string SamFile::Flag[2];
rowlen SamFile::Offset5_9 = 0;
rowlen SamFile::ReadStartPos =	0;

// Creates and initialize line write buffer.
//	@rQualPatt: Read quality pattern, or NULL
void SamFile::InitBuffer(const char* rQualPatt)
{
	LineFill(0, TAB);
	rowlen startPos = ReadStartPos + Read::Len + 1;
	
	// Read quality line
	if(rQualPatt)	LineFill(startPos, rQualPatt, Read::Len);
	else			LineFill(startPos, Read::SeqQuality, Read::Len);	
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
//	@cSizes: chrom sizes
//	@commandLine: command line
void SamFile::CreateHeader(const ChromSizes& cSizes, const string& commandLine)
{
	SetHeaderLine("HD", "VN", "1.0", "SO", "unsorted");
	for(ChromSizes::cIter it=cSizes.cBegin(); it!=cSizes.cEnd(); it++)
		SetHeaderLine("SQ", "SN", Chrom::AbbrName(CID(it)), "LN", NSTR(it->second));
	SetHeaderLine("PG", "ID", Product::Title, "PN", Product::Title, false);
	SetHeaderLine(NULL, "VN", Product::Version, NULL, strEmpty, false);
	SetHeaderLine(NULL, "CL", commandLine);
	Write();	// write header in case of multithread
}

#define CL_LEN 3	// length of the header line tag 'CL:'

// Initializes line write buffer and header and makes ready for writing
//	@commandLine: command line
//	@cSizes: chrom sizes
//	@rQualPatt: Read quality pattern, or NULL
void SamFile::InitToWrite(const string& commandLine, const ChromSizes&cSizes, const char* rQualPatt)
{
	if(!ReadStartPos) {		// first constructor call?
		if( OutFile::PairedEnd() )	{ Flag[0] = "99";	Flag[1] = "147"; }
		else						{ Flag[0] = "0";	Flag[1] = "16";	 }
		Comb5_6 = OutFile::MapQual + "\t" + NSTR(Read::Len) + "M"; // MAPping Quality + CIGAR: Read length

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
		CreateHeader(cSizes, commandLine);
	}
	InitBuffer(rQualPatt);
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

string OutFile::MapQual;	// the mapping quality
OutFile::eMode OutFile::Mode = mSE;

OutFile::AddReads OutFile::callAddRead[] =
	{ &OutFile::AddReadSE, &OutFile::AddReadPE, &OutFile::NoAddRead };

// Creates new instance for writing.
//	@fName: common file name without extention
//	@outType: types of output files
//	@mode: SE or PE
//	@fqQualPattFName: name of valid file with FQ quality pattern, or NULL
//	@mapQual: the mapping quality
//	@isZipped: true if output files should be zipped
OutFile::OutFile(const string& fName, eFormat outType, eMode mode,
	const char* fqQualPattFName, BYTE mapQual, bool isZipped) :
	_rQualPatt(NULL)
{
	_mode = Mode = mode;
	_fqFile1 = _fqFile2 = NULL;
	if( fqQualPattFName && (outType & (ofFQ | ofSAM)) ) {	// fill Read quality pattern
		TabFile file(fqQualPattFName);
		const char* line = file.GetLine();
		if(!line)	Err(Err::TF_EMPTY, fName.c_str(), "lines").Throw();
		_rQualPatt = new char[Read::Len + 1];
		_rQualPatt[Read::Len] = '\0';	// needed only to print in parameters info
		readlen lineLen = (readlen)file.LineLength();
		if(lineLen < Read::Len)
			memset(_rQualPatt, Read::SeqQuality, Read::Len);
		else	lineLen = Read::Len;
		memcpy(_rQualPatt, line, lineLen);
	}
	MapQual = BSTR(mapQual);

	if( outType & ofFQ ) {
		if( PairedEnd() ) {
			_fqFile1 = new FqFile(fName, isZipped, 1);
			_fqFile2 = new FqFile(fName, isZipped, 2);
		}
		else
			_fqFile1 = new FqFile(fName, isZipped);
	}
	_bedFile = outType & ofBED ? new BedRFile(fName, isZipped) : NULL;
	_samFile = outType & ofSAM ? new SamFile (fName, isZipped) : NULL;
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
	if(_rQualPatt)	delete[] _rQualPatt;
	if(_fqFile1)	delete _fqFile1;
	if(_fqFile2)	delete _fqFile2;
	if(_bedFile)	delete _bedFile;
	if(_samFile)	delete _samFile;
}

// Initializes buffers and makes ready for writing
//	@cSizes: chrom sizes or NULL
//	@commandLine: command line
void OutFile::Init(const ChromSizes* cSizes, const string& commandLine)
{
	if(_fqFile1)	_fqFile1->InitToWrite(_rQualPatt);
	if(_fqFile2)	_fqFile2->InitToWrite(*_fqFile1);
	if(_bedFile)	_bedFile->InitToWrite(commandLine);
	if(_samFile)	_samFile->InitToWrite(commandLine, *cSizes, _rQualPatt);
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
	if(_fqFile1)	_fqFile1->AddRead(rName, read, reverse);
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
		_fqFile1->AddRead(rName + Read::NmSuffMate1, read1, false);
		_fqFile2->AddRead(rName + Read::NmSuffMate2, read2, true);
	}
	if(_bedFile)	_bedFile->AddTwoReads(rName, pos, pos2);
	if(_samFile)	_samFile->AddTwoReads(rName, read1, read2, pos, pos2, fragLen);
	return 1;
}

void OutFile::Write() const
{
	if(_fqFile1)	_fqFile1->Write();
	if(_fqFile2)	_fqFile2->Write();
	if(_bedFile)	_bedFile->Write();
	if(_samFile)	_samFile->Write();
}

// Prints output file formats and sequencing mode
//	@signOut: output marker
void OutFile::PrintFormat(const char* signOut) const
{
	if(_fqFile1) {
		cout << signOut << "Output sequence: " << _fqFile1->FileName();
		if(_fqFile2)
			cout << SepCm << _fqFile2->FileName(); 
		cout << endl;
	}
	if(_bedFile || _samFile) {
		cout << signOut << "Output alignment: ";
		if(_bedFile) {
			cout << _bedFile->FileName();
			if(_samFile)	cout << SepCm;
		}
		if(_samFile)	cout << _samFile->FileName();
		cout << endl;
	}
	cout << signOut << "Sequencing: " << (PairedEnd() ? "paired" : "single") << "-end\n";
}

// Prints Read quality settins
//	@signOut: output marker
void OutFile::PrintReadQual(const char* signOut) const
{
	bool prMapQual = _samFile || _bedFile;
	cout << signOut << "Read quality: ";
	if(_fqFile1) {
		cout << "initial ";
		if(_rQualPatt)	cout << "pattern" << Equel << _rQualPatt;
		else			cout << Equel << Read::SeqQuality;
		if(prMapQual)	cout << SepSCl;
	}
	if(prMapQual)		cout << "mapping" << Equel << MapQual;
	cout << EOL;
}

/************************ end of class OutFile ************************/

//#define MIN_NUMB	100

//SamFragDistr::SamFragDistr(const string& fname)
//{
//	//ULONG cntLines;
//	TabFile file(fname, 9, true, TxtFile::READ, NULL, '@', true);
//	 
//	cout << "Fragments distribution in " << fname << endl;
//
//	int numbLen = 1000;
//	Array<int> numbers(numbLen);
//	int v, i, insCnt = 0;
//	bool checkPair = false;
//
//	const char* currLine;	// = file.GetFirstLine(&cntLines);
//	for( i=1; currLine=file.GetLine(); i++ )
//		if( i%2 && (v=file.IntField(8)) ) {
//			if(v<0)	v *= -1;
//			if(v < numbLen) {	numbers[v]++;	insCnt++; }
//		}
//	// print
//	//cout << EOL;
//	for(int k=0; k<numbLen; k++)
//		if(numbers[k])
//			cout << k << TAB << numbers[k] << EOL;
//	cout << "counted " << insCnt << " from " << i << BLANK << sPercent(insCnt, i, 3, 0, true);
//	cout << endl;
//}
