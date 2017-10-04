#include "common.h"
#include <sstream>
#ifdef OS_Windows
	#include <algorithm>
	#define SLASH '\\'		// standard Windows path separator
	#define REAL_SLASH '/'	// is permitted in Windows too
#else
	#define SLASH '/'	// standard Linux path separator
#endif

/************************ common Functions ************************/

// Gets number of digist in a integral value
//	@val: integral value
//	return: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCount (LLONG number)
{
	//BYTE res = BYTE(number<0);
	BYTE res = 0;
	for(; number; number/=10, res++);
	return res;
}

// Returns string represents the percent of part relatively total
//	@percent: value of percent
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default
//	@parentheses: if true parenthesize the value (default)
string	sPercent(float percent, BYTE precision, BYTE fieldWith, bool parentheses)
{
	float threshold = (float)pow(10.0, -precision);
	stringstream sout;
	if(parentheses)		sout << " (";
	sout << /*setfill(BLANK) <<*/ setw(fieldWith);
	if( percent && percent < threshold )
		sout << '<' << threshold;
	else {
		if(precision)	sout << setprecision(precision);
		sout << percent;
	}
	sout << PERS;
	if(parentheses)		sout << ')';
	return sout.str();
}

#if defined _WIGREG || defined _BIOCC

// Align position to the up or down resoluation level
// f.e. by resoluation==5 pos 102 -> 100+relative, pos 104 -> 105++relative
//	@pos: chromosome's position
//	@res: resoluation
//	@relative: 0 or 1
//	1 used for 1-relative position (the first base is 1, WIG)
//	0 used for 0-relative position (BED)
//	return: aligned position
chrlen AlignPos(chrlen pos, BYTE res, BYTE relative)
{
	short rest = pos % res;
	return rest - relative ?
		pos + relative - rest + (rest > res<<1 ? res : 0) :
		pos;
}

#endif

//size_t getAvailSystemMemory()
//{
//#ifdef OS_Windows
//	MEMORYSTATUSEX status;
//	status.dwLength = sizeof(status);
//	GlobalMemoryStatusEx(&status);
//	return (size_t)status.ullAvailPhys;
//#else
//	long pages = sysconf(_SC_PHYS_PAGES);
//	long page_size = sysconf(_SC_PAGE_SIZE);
//	return pages * page_size;
//#endif
//}

//#ifdef OS_Windows
//string Wchar_tToString(const wchar_t* wchar)
//{
//	string str = strEmpty;
//	while(*wchar)	str += char(*wchar++);
//	return str;
//}
//
//wchar_t* StringToWchar_t(const string &str, wchar_t* wchar)
//{
//	int i = (int)str.size();
//	wchar[i] = 0;
//	for(--i; i >= 0; i--)
//		wchar[i] = (wchar_t)str[i];
//	return wchar;
//}
//#endif

/************************ end of common Functions ************************/

/************************ class Options ************************/

// prints value in parentheses
#define PRINT_IN_PRTHS(v)	cout<<" ["<<(v)<<']'

const char*	OptTitle = "Option ";
const char*	optTitle = " option ";
const char* Ambiguous = "Ambiguous";
const char*	Default = " Default: ";
const char*	Missing = "missing ";


const char* Options::Booleans [] = {"OFF", "ON"};

const char* Options::_TypeNames [] = {
	StrEmpty.c_str(), "<name>", "<chars>", "<int>", "<float>", "<long>",
	StrEmpty.c_str(), StrEmpty.c_str(), StrEmpty.c_str(), StrEmpty.c_str()
};

const char Options::Option::EnumDelims [] = {'|', ','};

#define DESCR_SHIFT	3		// number of TABs before decsriptions
#define ENUM_REPLACE '?'	// symbol in description that is replaced by enum value

const string sValue = "value";
//const string sWrongValue = "wrong " + sValue;

// Checks is string represents digital value
bool isValidDigit(const char *str)
{
	if (*str == HPH)	++str;		// handle negative numbers
	if (!*str)	return false;	// handle empty string or just "-"
	// check for non-digit chars in the rest of the stirng.
	BYTE dotCnt = 0;
	for(; *str; ++str)
		if(*str == DOT)	
			if(dotCnt)	return false;
			else		dotCnt++;
		else if(!isdigit(*str))	return false;
	return true;
}

// Recursively prints string with replaced ENUM_REPLACE symbol by enum/combi value.
//	@buff: external buffer to copy and output temporary string
//	@vals: enum/combi values or NULL for other types
//	@cnt: external counter of enum/combi values
void PrintTransformDescr(char* buff, const char** vals, short* cnt)
{
	if(vals) {			// enum/combi?
		const char* subStr = strchr(buff, ENUM_REPLACE);
		if(subStr) {	// something to replace by enum value
			size_t strLen = subStr - buff;
			buff[strLen] = 0;
			cout << buff << vals[(*cnt)++];	// output substring and enum value
			buff[strLen] = ENUM_REPLACE;
			PrintTransformDescr(buff+strLen+1, vals, cnt);
		}
		else	cout << buff;
	}
	else	cout << buff;
}

// Recursively prints string with EOL inside as a set of left-margin strings
// Used to output aligned option descriptions
// First string is printed from current stdout position.
// Last substring doesn't include EOL.
//	@buff: external buffer to copy and output temporary string
//	@str: input string with possible EOLs
//	@subStr: substring of input string to the first EOL, or NULL if input string is not ended by EOL
//	@vals: enum/combi values or NULL for other types
//	@cnt: external counter of enum/combi values
void PrintSubLine(char* buff, const char* str, const char* subStr, const char** vals, short* cnt)
{
	if(subStr) {	// is substring ended by EOL exist?
		// form substring
		size_t strLen = subStr - str;
		strncpy(buff, str, strLen);
		buff[strLen] = 0;
		PrintTransformDescr(buff, vals, cnt);	// output enum values
		cout << EOL;
		for(BYTE t=0; t<DESCR_SHIFT; t++)	cout << TAB;
		str = subStr + 1;		// skip EOL
		subStr = strchr(str, EOL);
		PrintSubLine(buff, str, subStr, vals, cnt);
	}
	else {			// output rest of initial string without EOL
		strcpy(buff, str);
		PrintTransformDescr(buff, vals, cnt);
	}
}

// Sets option value.
//	@isWord: true if option is a word, false if option is a char
//	@opt: option without HYPHENs
//	@val: value of option
//	@isNextOpt: true if next parameter is option
//	Return: 0 if success, -1 if not found, 1 if option or value is wrong
int Options::Option::SetVal(bool isWord, char* opt, char* val, bool isNextOpt)
{
	if( isWord ) { if(Str==NULL || strcmp(Str, opt))	return -1; }
	else if(Char != *opt)	return -1;

	if(Sign.Is(Signs::Recogn))	return PrintAmbigOpt(isWord, opt, "Duplicated", NULL);
	// check value existence
	bool noVal = val==NULL || val[0]==HPH;	// true if no value
	if(ValRequired())	{ if(noVal)	return PrintWrongOpt(isWord, opt, NULL, sValue + " required"); }
	else if(!noVal&&isNextOpt)	return PrintWrongOpt(isWord, opt, NULL, sValue + " prohibited");

	Sign.MarkAs(Signs::Recogn);
	switch(ValType) {
		case tENUM:	return SetEnum(val) ? 0 : PrintWrongOpt(isWord, opt, val);
		case tCOMB:	return SetComb(val) ? 0 : PrintWrongOpt(isWord, opt, val);
		case tCHAR:	
			if(NVal != vUNDEF)
				return SetTriedDigit(*val, isWord);	// value is treated as int
													// and tCHAR where value is treated as string
		case tNAME: SVal = val;			return 0;
		case tHELP:	return PrintUsage(true);
		case tVERS:	return PrintVersion();
		default:
			if(isValidDigit(val))
				return SetTriedDigit(atof(val), isWord);	// numerical value
	}
	return PrintWrongOpt(isWord, opt, val);
}

// Check option for obligatory.
//	return: -1 if option is obligatory but not stated, otherwise 1
int Options::Option::CheckOblig()
{
	if( Sign.Is(Signs::Oblig) && ValRequired()
	&& ( (ValType == tNAME && SVal == NULL)
		|| (ValType != tNAME && NVal == vUNDEF) ) ) {
			cerr << Missing << "required option " 
					<< OptToStr(Char==HPH, Char==HPH ? Str : &Char) << EOL;
			return -1;
	}
	return 1;
}

// Checks limits and set numerical value
//	@val: numerical value
//	@isWord: true if option enters by long name
//	return: 1 if limits are exceeded, otherwise 0
int Options::Option::SetTriedDigit(double val, bool isWord)
{
	int outOfLimit = -1;
	if( val < MinNVal )		outOfLimit = 0;
	else if( val > MaxNVal )	outOfLimit = 1;
	if( outOfLimit >= 0 ) {
		const char* sign[] = {"less", "more"};
		cerr << OptTitle << OptToStr(isWord, isWord ? Str : &Char);
		cerr << MSGSEP_BLANK << sValue << setprecision(4) << BLANK << val
				<< " is " << sign[outOfLimit] << " than permissible " 
				<< (!outOfLimit ? MinNVal: MaxNVal ) << EOL;
		return 1;
	}
	NVal = val;
	return 0;
}

// Checks and sets enum option value.
//	@val: input value as string
//	return: true if success
bool Options::Option::SetEnum(char* val)
{
	if( ValRequired() ) {				// non-boolean template
		int ind = GetEnumInd(val);
		if( ind < 0 )	return false;
		NVal = ind + MinNVal;	// add possible minimum value as a shift
	}
	else
		NVal = !NVal;			// invers boolean value
	return true;
}

// Checks and sets enum option value.
//	@val: input value as C string
//	return: true if success
bool Options::Option::SetComb(char* vals)
{
	char* pdelim;	// a pointer to the first occurrence of delimiter COMMA in vals
	int	ind;
	bool ret = true;
	char delim = EnumDelims[1];

	NVal = 0;	// reset default value
	// run through given values
	for(char* val = vals; true; val = pdelim + 1) {
		pdelim = strchr(val, delim);
		if(pdelim)	*pdelim = '\0';		// temporary turn val into a pointer to C string
		ind = GetEnumInd(val);
		if(ind < 0)		ret = false;
		else// set bitwise val, in which running number of each bit (from right end)
			// corresponds to number of finded value in enum
			NVal = int(NVal) ^ (1<<ind);
		if(pdelim)	*pdelim = delim;	// give delimiter back
		else	break;					// no delimiter: last or single value
	}
	return ret;
}

// Prints option with double blank before if descr==true, single blank otherwise.
void Options::Option::Print(bool descr)
{
	BYTE cntChars = 0;
	if(descr)	{ cout << BLANK;	cntChars++;	}
	//Option& sopt = _Options[i];
	// PRINT OPTION
	cout << BLANK << HPH << Char;	cntChars += 3;
	if(Str)	{
		if(Char != HPH) {	cout << "|--";	cntChars += 3; }
		cout << Str;
		cntChars += strlen(Str);
	}
	// PRINT OPTION'S VALUE 
	if( ValType == tENUM ||  ValType == tCOMB )
		cntChars += PrintEnumVals();
	else {
		const char* tname = _TypeNames[ValType];
		cout << BLANK << tname;
		cntChars += 1 + strlen(tname);
	}
	// PRINT DESCSIPTION
	if( descr ) {
		// align description 
		short cnt = DESCR_SHIFT - cntChars / 8;	// 3*8: right boundary of descriptions
		if( !cnt )						cout << TAB;
		else if( cnt > 0 )
			for(BYTE i=0; i<cnt; i++)	cout << TAB;

		// print description
		cnt = 0;	// use as external enum counter
		char* buffer = new char[strlen(Descr) + 1];
		PrintSubLine(
			buffer,
			//new char[strlen(Descr) + 1],
			Descr,
			strchr(Descr, EOL), 
			ValType == tENUM || ValType == tCOMB ? (const char**)SVal : NULL, 
			&cnt
		);
		delete [] buffer;

		// print DEFAULT
		if( ValRequired() && NVal != vUNDEF )
			if( ValType == tENUM || ValType == tCOMB ) {
				if( NVal >= MinNVal )	// do not print default if it is never set by user
					PRINT_IN_PRTHS(	((char**)SVal)
						[int(NVal)-int(MinNVal)] );	// offset by min enum value
			}
			else if(ValType == tCHAR)
				PRINT_IN_PRTHS(char(NVal));
			else
				PRINT_IN_PRTHS(NVal);
		else if( SVal != NULL )
			PRINT_IN_PRTHS(ValType == tENUM ? "NONE" : SVal);
		cout << EOL;
	}
}

// Prints enum or combi values
//	return: number of printed symbols
BYTE Options::Option::PrintEnumVals()
{
	if( ValRequired() ) {
		char** vals = (char**)SVal;
		BYTE wCnt = BYTE(MaxNVal),
			 cCnt = BYTE(strlen(vals[0]));
		cout << " <" << vals[0];
		for(BYTE i=1; i<wCnt; i++) {
			cout << EnumDelims[ValType-tENUM] << vals[i];
			cCnt += strlen(vals[i]);
		}
		cout << '>';
		return cCnt + wCnt + 2;
	}
	return 0;
}

// Performs a case-insensitive search of given string value among enum values.
//	@val: input value as string
//	return: index of finded value in enum, or -1 if the value is not present in enum
int Options::Option::GetEnumInd (char* val)
{
	char** templ = (char**)SVal;
	BYTE i, cnt = BYTE(MaxNVal);
	for(i=0; i<cnt; i++)
		if( !_stricmp(val, templ[i]) )
			return i;
	return -1;
}

// Set value of option 
//	@isWord: true if option is a word, false if option is a char
//	@opt: option with HYPHENs
//	@val: value of option
//	@isNextOpt: true if next parameter is option
//	@argIndex: the current index in argc; is increased by 1 in case of required value
//	return: 0 if success, 1 otherwise
int Options::SetOption (bool isWord, char* opt, char* val, bool isNextOpt, int *argIndex)
{
	int res;

	opt += int(isWord)+1;
	if(isWord) { if(strlen(opt) == 1)
		return PrintAmbigOpt(isWord, opt, Ambiguous, "excess '-'?"); }
	else if(strlen(opt) > 1)
		return PrintAmbigOpt(isWord, opt, Ambiguous, "forgot '-'?");

	for(int i=0; i<_OptCount; i++) {
		res = _Options[i].SetVal(isWord, opt, val, isNextOpt);
		if(!res) {
			if( _Options[i].ValRequired() )	(*argIndex)++;
			return 0;
		}
		if(res>0)		return 1;
	}
	cerr << "wrong option: " << OptToStr(isWord, opt) << EOL;
	return 1;
}

// Prints version
//	return: always 1
int	Options::PrintVersion()
{
	cout<< Product::Version
#ifndef _NO_ZLIB
		<< "\tzlib "ZLIB_VERSION
#endif
		<< endl;
	return 1;
}

// Ouptuts option with error message to cerr
//	@isWord: true if option is long
//	@opt: option
//	@val: value or NULL
//	@msg: error message about value
int Options::PrintWrongOpt(bool isWord, const char* opt, const char* val, const string msg)
{
	cerr << OptTitle << OptToStr(isWord, opt)
		 << MSGSEP_BLANK << (msg == StrEmpty ? "wrong " + sValue : msg);
	if( val ) cerr << BLANK << val;
	cerr << EOL;
	return 1;
}

// Ouptuts ambiguous option with error message to cerr
//	@isWord: true if option is long
//	@opt: option
//	@headMsg: message at the beginning
//	@tailMsg: message at the end or NULL
int Options::PrintAmbigOpt(bool isWord, const char* opt, const char* headMsg, const char* tailMsg)
{
	cerr << headMsg << optTitle << OptToStr(isWord, opt);
	if( tailMsg )	cerr << MSGSEP_BLANK << tailMsg;
	cerr << EOL;
	return 1;
}


// Check obligatory options and output message about first absent obligatory option.
//	return: -1 if some of obligatory options does not exists, otherwise 1
int Options::CheckObligs()
{
	for(int i=0; i<_OptCount; i++)
		if( _Options[i].CheckOblig() < 0 )	return -1;
	return 1;
}

// Prints 'usage' information
//	@title: if true prints title before information
//	return: 1 if title is settinf to true, 0 otherwise
int Options::PrintUsage (bool title)
{
	BYTE i, k;
	if( title )		cout << Product::Descr << endl << endl;
	
	// output 'Usage' section
	cout << "Usage:";
	for(k=0; k<_UsageCount; k++) {
		cout << TAB << Product::Title << " [options]";
		const Usage& usage = _Usages[k];
		// output required options
		for(i=0; i<_OptCount; i++)
			if( _Options[i].Sign.Is(Signs::Oblig) )
				_Options[i].Print(false);
		// output parameters
		if( usage.OptVal != vUNDEF )
			_Options[usage.OptVal].Print(false);	// output option value
		else
			cout << usage.Text;						// output text
		cout << endl;
	}
	cout << endl;

	// output options section
	for(k=0; k<_GroupCount; k++) {
		cout << _OptGroups[k] << ":\n";
		for(i=0; i<_OptCount; i++)
			if( _Options[i].OptGroup == k )
				_Options[i].Print(true);
	}
	return int(title);
}

// Returns command line.
//	@argc: count of main() parameters
//	@argv: array of main() parameters
string const Options::CommandLine(int argc, char* argv[])
{
	ostringstream oss;
	int i;
	for (i = 0; i < argc-1; i++)
		oss << argv[i] << BLANK;
	oss << argv[i];
	return oss.str();
}

// Parses and checks main() parameters and their values. Output message if some of them is wrong.
//	@argc: count of main() parameters
//	@argv: array of main() parameters
//	@obligPar: name of required application parameter or NULL if not required
//	return: index of first parameter (not option) in argv[], argc if it is absent, negative if  tokenize complets wrong
int Options::Tokenize(int argc, char* argv[], const char* obligPar)
{
	int i, res = 1;
	for (i = 1; i < argc; i++) {	// Remember argv[0] is the path to the program
		if( argv[i][0] != HPH )	{
			if( i < argc-1 				// not a last option or parameter
			&& argv[i+1][0] == HPH ) {	// next word is an option
				cerr << argv[i] << ": neither option nor parameter"  << EOL;
				res = -1;
			}		
			break;
		}
		if( SetOption(
			argv[i][1] == HPH, argv[i], 
			i+1 < argc ? argv[i+1] : NULL,
			i+2 < argc ? argv[i+2][0]==HPH : false,
			&i) )
		{	res = -1; break; }
	}
	if( res > 0 )	res = CheckObligs();	// check required options
	if( res > 0 && obligPar && i == argc ) {
		cerr << Missing << obligPar << EOL;
		res = -1;
	}
	return i * res;
}

//void Options::GetOpt(int i)
//{
//	_Options[i].Print(true);
//}
/************************ end of class Options ************************/

/************************ class Err ************************/

const char* Err::_msgs [] = {
/* NONE */		"WARNING: ",//WARNING,
/* P_MISSED */	"missing",
/* F_NON */		"no such file",
/* FD_NON */	"no such file or directory",
/* F_MEM */		"memory exceeded",
/* F_OPEN */	"could not open",
/* F_CLOSE */	"could not close",
/* F_READ */	"could not read",
/* F_WRITE */	"could not write",
/* F_BIGLINE */	"buffer is less than length of line",
/* F_NOREADREC*/"attempt to get record's info without reading record",
/* FZ_MEM */	"wrong setting internal gzip buffer",
/* FZ_OPEN */	"wrong reading mode ALL for zipped file",
/* FZ_BUILD */	"this build does not support zipped files",
#ifndef _FQSTATN
/* TF_FIELD */	"wrong format: number of fields is less than expected",
/* TF_SPEC */	"wrong line format",
/* TF_EMPTY */	"no",
/* BP_BADEND */	"'start' position is equal or more than 'end'",
/* BP_NEGPOS */	"negative position is not allowed",
/* BP_EXCEED */	"position exceeds chromosome length",
#ifdef _BEDR_EXT
/* BR_RNAME */	"wrong read name format:",
#endif
///* BR_SIZE */	"different size of read",
/* FA_LONGLEN */"length of chromosome is more than ULONG_MAX",
#endif
#if defined(_ISCHIP) || defined(_FQSTATN)
/* FQ_HEADER */	"no '@' marker; miss header line",
/* FQ_HEADER2 */"no '+' marker; miss second header line",
#elif defined _DENPRO || defined _BIOCC
/* ARR_OUTRANGE */	" is out of range",
/* SUM_EXCEED */	"exceeded digital limit while S calculated. Anormous density. Calculate P only",
#endif
/* EMPTY */		""
};

// Initializes _outText by cstring contained message kind of "<@sender>: <@text> <@specifyText>".
//void Err::set_message(const string& sender, const char *text, const char *specifyText)
//{
//	size_t size = sender.size() + strlen(text) + 1;
//	//size_t size0 = size;
//	if( sender.size() )	size += strlen(MSGSEP_BLANK);
//	if( specifyText )	size += strlen(specifyText) + 1;
//	_outText = new char[size];
//	strcpy(_outText, sender.c_str());
//	if( sender.size() )
//		strcat(_outText, MSGSEP_BLANK);
//	strcat(_outText, text);
//	if( specifyText ) {
//		strcat(_outText, sBLANK);
//		strcat(_outText, specifyText);
//	}
//}

void Err::set_message(const char* sender, const char *text, const char *specifyText)
{
	size_t senderLen = sender ? strlen(sender) : 0;
	size_t outLen = senderLen + strlen(text) + 1;
	if(senderLen)	outLen += strlen(MSGSEP_BLANK);
	if(specifyText)	outLen += strlen(specifyText) + 1;
	_outText = new char[outLen];
	memset(_outText, '\0', outLen);
	if(senderLen) {
		strcpy(_outText, sender);
		strcat(_outText, MSGSEP_BLANK);
	}
	strcat(_outText, text);
	if(specifyText) {
		strcat(_outText, sBLANK);
		strcat(_outText, specifyText);
	}
}

Err::Err(const Err & src)
{
	_code = src._code;
	int size = strlen(src._outText) + 1;
	_outText = new char[size];
	strcpy(_outText, src._outText);
	//_specifyText = src._specifyText;
}

void Err::Throw(bool throwException, bool endOfLine) {
	if( throwException )	
		throw *this;
	else {		
		dout << _outText;
		if( endOfLine )		dout << EOL;
	}
}

// Outputs warning with prefix "WARNING" and additional text, if it is setting.
void Err::Warning(string const& addText) {
	dout << _msgs[ErrWARNING] << _outText;
	if( !addText.empty() )	
		dout << addText;
	dout << EOL;
}

/************************ end of class Err ************************/

/************************ class FileSystem ************************/

// Returns true if file system's object exists
//	@name: object's name
//	@st_mode: object's system mode
bool FS::IsExist(const char* name, int st_mode)
{
	struct_stat64 st;
	return ( !_stat64(name, &st) && st.st_mode & st_mode );
}

// Checks if file system's object doesn't exist
//	@name: object's name
//	@st_mode: object's system mode
//	@throwExcept: if true throws excwption,
//	otherwise outputs Err message as warning without EOL
//	@ecode: error's code
//	return: true if file or directory doesn't exist
bool FS::CheckExist	(const char* name,  int st_mode, bool throwExcept, Err::eCode ecode)
{
	if( IsExist(name, st_mode) )	
		return false;
	Err(ecode, name).Throw(throwExcept);
	return true;
}

// Searches through a file name for the any extention (slash|back-slash insensible).
//	@fname: file name
//	return: the index of the DOT mathed extention; otherwise npos
size_t FS::GetExtPos(const string &fname) {
	size_t pos = fname.find_last_of(DOT);
	if( pos != string::npos ) {
		if( pos==0 )					return string::npos;	// ./name
		if( pos==1 && fname[0]==DOT )	return string::npos;	// ../name
	}
	return pos;
}

// Gets size of file or -1 if file doesn't exist
LLONG FS::Size (const char* fname)
{
	struct_stat64 st;
	return _stat64(fname, &st) == -1 ? -1 : st.st_size;
}

// Gets real size of zipped file  or -1 if file cannot open; limited by UINT
LLONG FS::UncomressSize	(const char* fname)
{
	FILE *file = fopen(fname, "rb");	// "read+binary"
	if( file == NULL )		return -1;
	BYTE sz[4] = {0,0,0,0};
	_fseeki64(file, -4, SEEK_END);
	fread(sz, 1, 4, file);
	fclose(file);
	return (sz[3] << 3*8) + (sz[2] << 2*8) + (sz[1] << 8) + sz[0];
}

// Returns true if file has a specified  extension.
//	@fname: file name
//	@ext: extension includes dot symbol and can be composite
bool FS::HasExt	(const string& fname, const string& ext)
{ 
	size_t pos = fname.find(ext);
	return pos == string::npos ? false : fname.size() - pos == ext.size();
}

// Returns string containing real file extension (without zip extention).
//	@fname: pointer to the file name
//	return: string containing real file extension or empty string if no real extention
string const FS::GetExt(const char* fname) {
	const char * pdot = strrchr(fname, DOT);
	if( !pdot )		return StrEmpty;
	if( strcmp(pdot, ZipFileExt.c_str()) )
		return string(pdot+1);				// no zip extention
	const char * pprevdot = pdot - 1;
	for(; pprevdot >= fname; pprevdot--)	// find previous DOT
		if( *pprevdot == DOT )	
			break;
	return pprevdot+1 == fname ? StrEmpty : string(pprevdot+1, pdot-pprevdot-1);
}

// Returns file name without extentiom (slash|back-slash insensible)
string const FS::FileNameWithoutExt (const string& fname)
{
	size_t pos = GetExtPos(fname);
	return pos == string::npos ? fname : fname.substr(0, pos);
}

// Returns short file name by long one
//	@fname: long file name
string const FS::ShortFileName (const string& fname)
{
#ifdef OS_Windows
	if(fname.find(REAL_SLASH) != string::npos) {
		string tmp(fname);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp.substr(tmp.find_last_of(SLASH) + 1);
	}
#endif
	return fname.substr(fname.find_last_of(SLASH) + 1);
}

// Returns directory name by long file name
//	@fname: long file name
//	@addSlash: true if slash sould be added at the end
string const FS::DirName (const string& fname, bool addSlash)
{
#ifdef OS_Windows
	if(fname.find(REAL_SLASH) != string::npos) {
		string tmp(fname);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp.substr(0, tmp.find_last_of(SLASH) + int(addSlash));
	}
#endif
	return fname.substr(0, fname.find_last_of(SLASH) + int(addSlash));
}

// Returns the name of last subdirectory by long file name
//	@fname: long file name
string const FS::LastSubDirName (const string& fname)
{
	const string& dir = FS::DirName(fname, false);
	size_t pos = dir.find_last_of(SLASH);
	return pos == string::npos ? dir :
		dir.substr(dir.substr(0, pos).length() + 1);
}

// Returns the name ended by slash without checking the name
string const FS::MakePath(const string& name)
{
#ifdef OS_Windows
	if(name.find(REAL_SLASH) != string::npos) {
		string tmp(name + SLASH);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp;
	}
#endif
	return name + SLASH;
}

#ifndef _WIGREG
// Fills external vector of strings by file's names found in given directory
// Implementation depends of OS.
//	@files: external vector of strings that should be filled by file's names
//	@dirName: name of directory
//	@fileExt: file's extention as a choosing filter
//	@all: true if all files with given extention should be placed into external vector,
//	otherwise only one (any)
//	return: true if files with given extention are found
bool FS::GetFiles	(vector<string>& files, const string& dirName,
	const string& fileExt, bool all)
{
	BYTE count = 0;
#ifdef OS_Windows
	string fileTempl = FS::MakePath(dirName) + '*' + fileExt;
	WIN32_FIND_DATA ffd;

	HANDLE hFind = FindFirstFile( fileTempl.c_str(), &ffd );
	if( hFind == INVALID_HANDLE_VALUE )
		return false;
	if( all ) {
		// count files to reserve files capacity
		do	count++;
		while (FindNextFile(hFind, &ffd));
		files.reserve(count);
		// fill files
		hFind = FindFirstFile( fileTempl.c_str(), &ffd );
		do	files.push_back( string(ffd.cFileName) );
		while (FindNextFile(hFind, &ffd));
	}
	else
		files.push_back( string(ffd.cFileName) );
	FindClose(hFind);
	return true;
#else
	struct dirent *entry;
	DIR *dir = opendir(dirName.c_str());	// doesn't need to check because of programmes options inspection
	if( all ) {
		// count all files to reserve files capacity
		while( readdir(dir) )	count++;
		closedir(dir);
		if( !count )	return false;
		files.reserve(count);
		dir = opendir(dirName.c_str());
	}
	// fill files
	string name;
	while( entry = readdir(dir) ) {
		name = string(entry->d_name);
		if( HasExt(name, fileExt) ) {
			files.push_back(name);
			if( !all )	break;
		}
	}
	closedir (dir);
	return files.size() > 0;
#endif	// OS_Windows
}
#endif	// _WIGREG
/************************ end of class FileSystem ************************/

/************************  class Timer ************************/
bool	Timer::Enabled = false;
clock_t	Timer::_StartCPUClock;

// Prints elapsed time interval
//	@title: string printed before time output
//	@elapsed: elapsed time interval
//	@parentheses: if true then output time in parentheses
//	@isCarrgReturn: if true then ended output by EOL
void Timer::PrintElapsed(const char *title, long elapsed, bool parentheses, bool isCarrgReturn)
{
	if( title )			dout << title;
	if( parentheses )	dout << '(';
	long hours = elapsed/60/60;
	if( hours )
		dout << setfill('0') << setw(2) << hours << COLON;
	dout<< setfill('0') << setw(2) << (elapsed/60%60) << COLON
	    << setfill('0') << setw(2) << (elapsed%60);
	if( parentheses )	dout << ')';
	if( isCarrgReturn ) {
		dout << EOL;
		fflush(stdout);
	}
}

// Stops enabled CPU timer and print elapsed time
//	@isCarrgReturn: if true then ended output by EOL
void Timer::StopCPU(bool isCarrgReturn)
{
	if( Enabled )
		PrintElapsed("CPU: ", (clock()-_StartCPUClock)/CLOCKS_PER_SEC, false, isCarrgReturn);
}

// Stops enabled timer and print elapsed time with title
//	@title: string printed before time output
//	@parentheses: if true then output time in parentheses
//	@isCarrgReturn: if true then ended output by EOL
void Timer::Stop(const char *title, bool parentheses, bool isCarrgReturn)
{
	if( Enabled ) {
		time_t stopTime;
		time( &stopTime );
		PrintElapsed(title, (long)difftime(stopTime, _startTime), parentheses, isCarrgReturn);
	}
}
/************************  end of class Timer ************************/

/************************  class CPU_Timer ************************/
//bool CPU_Timer::Enabled = false;
//clock_t CPU_Timer::_startTime;
//
//void CPU_Timer::Stop(const char *title, bool isCarrgReturn)
//{
//	if( Enabled )
//		Timer::PrintElapsed(title, (clock()-_startTime)/CLOCKS_PER_SEC, false, isCarrgReturn);
//}
/************************  end of class CPU_Timer ************************/

#ifdef _MULTITHREAD
/************************  class Mutex ************************/
pthread_mutex_t	Mutex::_mutexes[Mutex::Count];

void Mutex::Init() {
	for(BYTE i=0; i<Count; i++)
#ifdef OS_Windows
		InitializeCriticalSection(&_mutexes[i]);
#else
		pthread_mutex_init(&_mutexes[i], NULL);
#endif
}

void Mutex::Finalize() {
	for(BYTE i=0; i<Count; i++)
#ifdef OS_Windows
		DeleteCriticalSection(&_mutexes[i]);
#else
		pthread_mutex_destroy(&_mutexes[i]);
#endif
}

void Mutex::Lock(const eType type) {
#ifdef OS_Windows
	EnterCriticalSection(&_mutexes[type]);
#else
	pthread_mutex_lock(&_mutexes[type]);
#endif
}

void Mutex::Unlock(const eType type) {
#ifdef OS_Windows
	LeaveCriticalSection(&_mutexes[type]);
#else
	pthread_mutex_unlock(&_mutexes[type]);
#endif
}

/************************  end of class Mutex ************************/

/************************  class Thread ************************/
Thread::Thread(retThreadValType(
	#ifdef OS_Windows
	__stdcall 
	#endif
	*proc)(void*), void *arglist)
{
#ifdef OS_Windows
	//_thread = (HANDLE) _beginthread(proc, 0, arglist);	// gets unstable call WaitForSingleObject
	_thread = (HANDLE) _beginthreadex(NULL, 0, proc, arglist, 0, NULL);
#else
	pthread_create(&_thread, NULL, proc, arglist);
	//int code = pthread_create(&_thread, NULL, proc, arglist);
	//if (code) { 
	//	char buf[256]; 
	//	strerror_r(code, buf, sizeof buf);
	//	cout << "Thread constructor: " << buf << endl;
	//}
#endif
}
/********************  end of class Thread *********************/
#endif	// _MULTITHREAD

#ifdef _BIOCC
/********************  class CCaggr *********************/
// Increases P values in consideration of means
//	@x: first value to increase
//	@y: second value to increase
void CCaggr::IncreasePVars(double x, double y) {
	x -= Mean1;
	y -= Mean2;
	VarP1 += x * x;
	VarP2 += y * y;
	CovP += x * y;
}

// Increases S value with check-up for exceeding
//	@i: index of value on VarS
//	@val: value to increase
//	@cID: chrom's ID, needed for exception only
//	@return: false if all right 
bool CCaggr::IncreaseSVar(BYTE i, ULLONG val, chrid cID)
{
	if(val) {
		check = VarS[i];
		if((VarS[i] += val) <= check) {
			Err(Err::SUM_EXCEED, Chrom::AbbrName(cID)).Warning();
			VarS[0] = Undef;
			return true;
		}
	}
	return false;
}

// Increases S values with check-up for exceeding
//	@x: first value to increase
//	@y: second value to increase
//	@cID: chrom's ID, needed for exception only
//	@return: false if all right 
bool CCaggr::IncreaseSVars(ULLONG x, ULLONG y, chrid cID)
{
	if( IncreaseSVar(0, x*x, cID) )		return true;
	if( IncreaseSVar(1, y*y, cID) )		return true;
	if( IncreaseSVar(2, x*y, cID) )		return true;
	return false;
}
/********************  end of class CCaggr *********************/
#endif	// _BIOCC

/************************ class Chrom ************************/

const char*		Chrom::Abbr = "chr";
const char*		Chrom::UndefName = "UNDEF";
const string	Chrom::Title = "chromosome";
const BYTE		Chrom::MaxAbbrNameLength = BYTE(strlen(Chrom::Abbr)) + MaxShortNameLength;
const BYTE		Chrom::MaxNamedPosLength = 
	BYTE(strlen(Chrom::Abbr)) + MaxShortNameLength + CHRLEN_CAPAC + 1;

chrid Chrom::_cID = UnID;	// user-defined chrom ID

// Sets chromosome's ID stated by user with validation.
//	@cID: chromosome's ID
//	return: true if the check has passed, otherwise do not set and print message to cerr
bool Chrom::SetStatedID(chrid cID)
{
	if( cID>64 && cID!=M && cID!=X && cID!=Y ) {
		cerr << Name(cID) << ": wrong " << Title << "'s name\n";
		return false;
	}
	_cID = cID;
	return true;
}

// Returns a pointer to the first occurrence of C sunstring. Recurcive.
//	@str: C string to find in
//	@templ: C string to be located
//	@templLen: length of templ (extern because of avoiding recursive recalculate)
const char* SubStr(const char* str, const char* templ, int templLen)
{
	str = strchr(str, *templ);
	if( str )
		for(short i=1; i<templLen; i++)
			if( *++str != templ[i] )
				return SubStr(str, templ, templLen);
	return str;
}

chrid Chrom::ID(const char* cName, size_t prefixLen)
{
	if( !cName )		return UnID;
	cName += prefixLen;									// skip prefix
	//if(*cName == M || strchr(cName, USCORE))	return UnID;
	if( strchr(cName, USCORE))	return UnID;			// exclude chroms with '_'
	if(*cName <= '9')			return atoi(cName);		// numeric chromosome
	return islower(*cName) ? toupper(*cName) : *cName;	// letter's chromosome
}

// Locate chrom number in string.
//	@str: string checked for chrom number
//	return: pointer to the chrom number in str, or a null pointer if Chrom::Abbr is not part of str.
const char* Chrom::FindNumb(const char* str) {
	const char* substr = strstr(str, Abbr);
	return substr ? (substr+strlen(Abbr)) : NULL;
}

// Returns the length of prefix (substring before short chromosome's name)
//	@cLongName: long chromosome's name
//	return: length of substring before short chromosome's name, or -1 if short name is not finded
short Chrom::PrefixLength(const char* cLongName)
{
	// search from the beginning of the cName because
	// it simpler to find the first digit for multidigit number

	// start search from the occurrence of 'chr'
	const char* str = SubStr(cLongName, Abbr, strlen(Abbr));
	if( str )
		for(; *str; str++)
			if( isdigit(*str) || isupper(*str) )
				return str - cLongName;
	return -1;
}

//char* Chrom::LongToShortName(char* name) 
//{
//	short shift = PrefixLength(name);
//	return shift > 0 ? name + shift : NULL;
//}

/************************ end of class Chrom ************************/

#if !defined _WIGREG && !defined _FQSTATN

#if defined _ISCHIP || defined _BEDR_EXT
/************************ class Read ************************/

readlen	Read::Len;				// length of Read
const char	 Read::Strand[] = { '+', '-' };
const char	 Read::NmPosDelimiter = '-';
const char*  Read::NmNumbDelimiter = ":N:";
const string Read::NmSuffMate1 = "/1";
const string Read::NmSuffMate2 = "/2";
BYTE	Read::OutNameLength = 0;
	//Read::Name().length() + 1 +		// + 1 delimiter
	//Chrom::MaxAbbrNameLength +	1 +	// length of chrom's name + 1 delimiter
	//2*CHRLEN_CAPAC + 1 +			// length of PE Read name (the longest) + 1 delimiter
	//NmSuffMate1.length();			// length of Mate suffix

#endif

#ifdef _ISCHIP

char	Read::SeqQuality;		// the quality values for the sequence (ASCII)
BYTE	Read::MapQuality = 1;	// the mapping quality
Read::rNameType	Read::NameType;		// type of name of Read in output files
short	Read::LimitN = vUNDEF;	// maximal permitted number of 'N' in Read or vUNDEF if all
ULONG	Read::MaxCount;			// up limit of writed Reads
ULONG	Read::Count = 0;		// counter of total writed Reads
const char*	Read::NmDelimiter = NULL;
const char Read::ToUp	= 'a' - 'A';
const char Read::Complements[] = {'T',0,'G',0,0,0,'C',0,0,0,0,0,0,'N',0,0,0,0,0,'A'};

void Read::Init(readlen rLen, rNameType nmType, char seqQual, BYTE mapQual, short limN, ULONG maxCnt)
{
	OutNameLength = 
		Read::Name().length() + 1 +		// + 1 delimiter
		Chrom::MaxAbbrNameLength + 1 +	// length of chrom's name + 1 delimiter
		2*CHRLEN_CAPAC + 1 +			// length of PE Read name (the longest) + 1 delimiter
		NmSuffMate1.length();			// length of Mate suffix
	Len = rLen;
	NameType = nmType;
	SeqQuality = seqQual;
	MapQuality = mapQual;
	if(limN < rLen)		LimitN = limN;
	MaxCount = maxCnt;
	if( nmType == nmPos )		NmDelimiter = NmNumbDelimiter + NmDelimiterShift;
	else if( nmType == nmNumb )	NmDelimiter = NmNumbDelimiter;
}

// Copies complemented Read.
void Read::CopyComplement(char* dst, const char* src)
{
	//for(dst += Read::Len-1; *dst != EOL; dst--, src++)
	//	*dst = Complements[*src - 'A' - (*src >= 'a' ? ToUp : 0)];
	for(char i=Read::Len-1; i>=0; i--, src++)
		dst[i] = Complements[*src - 'A' - (*src >= 'a' ? ToUp : 0)];
	//dst[Read::Len] = EOL;
}

// Checks Read for number of 'N'
//	@read: checked Read
//	return: -1 if Read is NULL, 0 if N limit is exceeded, 1 if success
int Read::CheckNLimit(const char* read)
{
	if( !read )		return -1;
	if( LimitN != vUNDEF ) {
		readlen cntN = 0, i = Len;
		for(read += i-1; i; i--, read--)
			if( *read == cN	&& ++cntN > LimitN )
				return 0;
	}
	return 1;
}

// Prints Read values - parameters.
void Read::Print()
{
	cout << "Read: length = " << int(Len)
		 << GroupParSep << "name includes a " << (NameType==nmPos ? "position" : "number")
		 << GroupParSep << "FQ quality = " << SeqQuality
		 << GroupParSep << "map quality = " << int(MapQuality)
		 << GroupParSep << "N-limit = ";
	if( LimitN == vUNDEF )	cout << Options::GetBoolean(false);
	else					cout << LimitN;
	cout << GroupParSep << "limit = " << MaxCount;
}

#endif

/************************ end of struct Read ************************/

#endif
