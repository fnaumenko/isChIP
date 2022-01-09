/**********************************************************
common.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 07.01.2022
-------------------------
Provides common functionality
***********************************************************/

#include "common.h"
#include <sstream>
#include <algorithm>	// std::transform, REAL_SLASH
#ifdef OS_Windows
	#define SLASH '\\'		// standard Windows path separator
	#define REAL_SLASH '/'	// is permitted in Windows too
#else
	#define SLASH '/'		// standard Linux path separator
#endif

/************************ common Functions ************************/

// Returns number of ones in an bynary integer
int OnesCount(int n)
{
	int cnt = 0;
	//for (; n; n >>= 1)  cnt += n & 1;		// 11001: 5 cycles
	for (; n; n &= n - 1)   cnt++;			// 11001: 3 cycles
	return cnt;
}

// Returns right position of right one in an bynary integer
int RightOnePos(int n)
{
	int pos = 0;
	for (; n ^ 1; n >>= 1)	pos++;
	return pos;
}

// Gets number of digist in a integral value
//	@val: integral value
//	@isLocale: if true then adds number of '1000' separators
//	return: number of digist without minus symbol or 0 if value is 0
int DigitsCount (LLONG val, bool isLocale)
{
	int res = 0;
	for(; val; val/=10, res++);
	if(isLocale)	res += (res-1)/3;
	return res;
}

// Returns string represents the percent of part relatively total
//	@percent: value of percent
//	@precision: count of fractional digits; 
//	if count of value's mapped digits is more then that, printed "<X%", or exactly by default
//	@fieldWith: displayed width of value and '%' or '<' sign (excluding parentheses), or exactly if 0;
//	@parentheses: if true then parenthesize the value (not considering fieldWith)
string	sPercent(float val, BYTE precision, BYTE fieldWith, bool parentheses)
{
	float threshold = (float)pow(10., -precision);
	stringstream ss;
	ss << SPACE;
	if(parentheses)		ss << '(';
	//if(fieldWith)	ss << setw(--fieldWith);	// decrease to account '%' sign
	if (val && val < threshold) {
		if (fieldWith) {
			int blankCnt = fieldWith - precision - 5;
			if (blankCnt > 0)	ss << setw(blankCnt) << SPACE;
		}
		ss << '<' << threshold;
	}
	else {
		if (precision && val) {
			if (val >= 100)	precision = 3;
			else if (val < 1)	ss << fixed;
			ss << setprecision(precision);
		}
		//if(val - int(val) != 0)	
		//	ss << fixed << setprecision(precision);
		if (fieldWith)	ss << setw(fieldWith-=2);	// decrease one for SPACE, one for '%'
		ss << val;
	}
	ss << PERS;
	if(parentheses)		ss << ')';
	return ss.str();
}

// Prints horizontal line
//	@w: width of line
void PrintHorLine(int w)
{
#ifdef OS_Windows
	wcout << setw(++w) << setfill(L'\304') << L'\n';
#ifdef _DUP_OUTPUT
	dout.ToFile(string(w, HPH) + LF);
#endif
#else
	for(int i=0; i<w; dout << "─", i++);	dout << LF;
	//dout << setw(++w) << setfill('─') << LF;
#endif
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

#ifdef _DUP_OUTPUT

// Open output file with given name
//	return: true if file is open
bool dostream::OpenFile(const string fname)
{ 
	if(!fname.length())	return false;
	file.open(fname.c_str());
	return true;
}

#endif

/************************ class Options ************************/
#define ENUM_REPLACE '?'	// symbol in description that is replaced by enum value

#define PRINT_IN_PRTHS(v)	cout<<" ["<<(v)<<']'	// prints value in parentheses
#define OPT_TO_STREAM(opt)	HPH<<(opt)

const char*	optTitle = "option ";
const char* Spotteruous = "Spotteruous";
const char*	Default = " Default: ";
const char*	Missing = "missing ";
const char* Warning = "WARNING: ";
const string sValue = "value";

const char* Options::sPrSummary = "print program's summary";
const char* Options::sPrTime = "print run time";
const char* Options::sPrUsage = "print usage information";
const char* Options::sPrVersion = "print program's version";
const char* Options::Booleans [] = {"OFF","ON"};
const char* Options::TypeNames [] = {
	NULL, "<name>", "<char>", "<int>", "<float>", "<long>", NULL, NULL,
	"<[int]:[int]>", "<[float]:[float]>", NULL, NULL, NULL
};

const char Options::Option::EnumDelims [] = {'|', ',', ':'};

// Checks digital value representation. Prints 'wrong val' message in case of failure
//	@str: defined (no NULL) string  representing digital value
//	return: true if digital value representation is correct
bool Options::Option::IsValidFloat(const char *str, bool isInt, bool isPair)
{
	char c = *str;
	const char* str0 = str;
	BYTE dotCnt = 0, eCnt = 0;

	if(!isdigit(c))						// check first char
		if(c == DOT)	dotCnt++;
		else if(c != HPH && c != PLUS)	goto end;
	for(str++; *str; str++)				// check next chars
		if((c=*str) == DOT)	
			if(dotCnt)	goto end;		// more than one dot
			else		dotCnt++;
		else if(tolower(c) == 'e')
			if(eCnt)	goto end;		// more than one 'e'
			else		eCnt++;
		else if(!isdigit(c))
			if(isPair && c==EnumDelims[2])	break;	// don't check the substring after ';'
			else goto end;				// wrong symbol
	if(isInt && dotCnt && !eCnt)		// e.g. 1.5e1 is acceptable as int
		cerr << Warning << ToStr(false) << SepSCl << "float value "
			 << (isPair ? "in " : strEmpty)
			 << str0 << " will be treated is integer\n";
	return true;
end:
	return !PrintWrong(str0);
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

// Recursively prints string with LF inside as a set of left-margin strings
// Used to output aligned option descriptions
// First string is printed from current stdout position.
// Last substring doesn't include LF.
//	@buff: external buffer to copy and output temporary string
//	@str: input string with possible LFs
//	@subStr: substring of input string to the first LF, or NULL if input string is not ended by LF
//	@vals: enum/combi values or NULL for other types
//	@cnt: external counter of enum/combi values
void PrintSubLine(char* buff, const char* str, const char* subStr, const char** vals, short* cnt)
{
	if(subStr) {	// is substring ended by LF exist?
		// form substring
		size_t strLen = subStr - str;
		strncpy(buff, str, strLen);
		buff[strLen] = 0;
		PrintTransformDescr(buff, vals, cnt);	// output enum values
		cout << LF;
		for(BYTE t=0; t<OPT_DESCF_TSHIFT; t++)	cout << TAB;
		str = subStr + 1;		// skip LF
		subStr = strchr(str, LF);
		PrintSubLine(buff, str, subStr, vals, cnt);
	}
	else {			// output rest of initial string without LF
		strcpy(buff, str);
		PrintTransformDescr(buff, vals, cnt);
	}
}

// Sets option value.
//	@opt: option
//	@isword: true if option is a word, false if option is a char
//	@val: value of option
//	@nextItem: next token after opt and val, or NULL
//	@argInd: the current index in argc; increased by 1 if value is accepted
//	Return: 0 if success, -1 if not found, 1 if option or value is wrong
int Options::Option::SetVal(const char* opt, bool isword, char* val, char* nextItem, int& argInd)
{
	if(isword) { 
		if(!Str || strcmp(Str, opt))	return -1;	// not this option
		Sign.MarkAs(eFlag::fWord);
	}
	else if(Char != *opt)	return -1;		// not this option

	if(Sign.Is(eFlag::fTrimmed))	return PrintAmbigOpt(opt, isword, "duplicated");

	//== check actual value
	const bool isValOblig = ValRequired() && !IsValEsc();
	const bool noRealVal =		// true if no actual value is submitted
		val == NULL								// no value
		|| (!nextItem && !isValOblig)			// nextItem is apps parameter
		|| (*val=='-' && !isdigit(val[1]));		// val is not negative value but next option
	
	if(noRealVal) {
		if(isValOblig)	return PrintWrong(NULL, sValue + " required");
		val = NULL;				// in case of 'false' value
	}
	else
		if(ValRequired())	argInd++;
		else 
			if(val && nextItem && *nextItem == HPH) {	// next token after val is option
				cerr << Warning, PrintWrong(NULL, sValue + " prohibited: " + string(val) + " ignored");
				argInd++;
			}

	//== set actual value
	Sign.MarkAs(eFlag::fTrimmed);
	switch(ValType) {
		case tNAME: SVal = noRealVal ? NULL : val;	return 0;
		case tENUM:	return SetEnum(val);
		case tCOMB:	return SetComb(val);
		case tCHAR:	if(NVal != NO_VAL && strlen(val) == 1)
						return SetTriedFloat(*val, MinNVal, MaxNVal);	// value is treated as int
					break;
		case tINT:
		case tFLOAT:
		case tLONG:	return noRealVal?
						0 :				// NVal is default val
						IsValidFloat(val, ValType==tINT) ?
							SetTriedFloat(float(atof(val)), MinNVal, MaxNVal) :
							1;
		case tHELP:	return PrintUsage(true);
		case tSUMM:	return PrintSummary(false);
		case tVERS:	return PrintVersion();
		default:	return SetPair(noRealVal ? NULL : val, ValType==tPR_INT);	// tPR_INT, tPR_FL
	}
	return PrintWrong(val);
}

// Check option for obligatory.
//	return: -1 if option is obligatory but not stated, otherwise 1
int Options::Option::CheckOblig() const
{
	if( Sign.Is(eFlag::fOblig) && ValRequired()
	&& ( (ValType == tNAME && !SVal) || (ValType != tNAME && NVal == NO_VAL) ) ) {
		cerr << Missing << "required option " << NameToStr(false) << LF;
		return -1;
	}
	return 1;
}

// Return option's signature as a string
//	@asPointed: true if returns signature as it was stated by user
inline string Options::Option::NameToStr (bool asPointed) const
{
	ostringstream ss;
	ss << HPH;
	if(asPointed) {
		if(!Sign.Is(eFlag::fWord))	{ ss << Char; return ss.str(); }
	}
	else if(Char != HPH) {
		ss << Char; 
		if(Str)	ss << '|' << HPH;
	}
	if(Str)		ss << HPH << Str;
	return ss.str();
}

// Returns string represented pair of value's separated by delimiter.
const string Options::Option::PairValsToStr(const pairVal* vals) const
{
	static const char* sAuto = "auto";
	ostringstream ss;
	ss << dec;
	if(vals->first == vUNDEF)	ss << sAuto;
	else						ss << vals->first;
	ss << EnumDelims[2];
	if(vals->second == vUNDEF)	ss << sAuto;
	else						ss << vals->second;
	return ss.str();
	//return static_cast<ostringstream & >( ostringstream() << dec
	//	<< vals->first << EnumDelims[2] << vals->second ).str();
}

// Checks limits and set numerical value
//	@val: value
//	return: 1 if limits are exceeded, otherwise 0
int Options::Option::SetTriedFloat(float val, float min, float max)
{
	if(!val && Sign.Is(eFlag::fAllow0))	min = 0;
	if(val < min || val > max) {
		cerr << optTitle << NameToStr(true) << SepSCl << sValue << setprecision(4) << SPACE << val
			 << " is out of available range [" << min << HPH << max << "]\n";
		return 1;
	}
	NVal = val;
	return 0;
}

// Checks and sets enum option value. Prints 'wrong val' message in case of failure
//	@val: input value as string or NULL if value is optional and not defined
//	return: 0 if success, 1 if wrong value
int Options::Option::SetEnum(const char* val)
{
	if (!ValRequired())			// boolean template
		NVal = !NVal;			// invers boolean value
	else if (val) {				// non-boolean template, value is defined
		int ind = GetEnumInd(val);
		if (ind < 0)	return PrintWrong(val);
		NVal = ind + MinNVal;	// add possible minimum value as a shift
	}
	// else val==NULL, so use default
	return 0;
}

// Checks and sets enum option value. Prints 'wrong val' message in case of failure
//	@val: input value as C string
//	return: 0 if success, 1 if wrong value
int Options::Option::SetComb(char* vals)
{
	char* delim;	// a pointer to the first occurrence of delimiter COMMA in vals
	int	ind;

	NVal = 0;	// reset default value
	// run through given values
	for(char* val = vals; true; val = delim + 1) {
		delim = strchr(val, EnumDelims[1]);
		if(delim)	*delim = cNULL;		// temporary cut 'val' to C string with single value
		ind = GetEnumInd(val);
		// set bitwise val, in which running number of each bit (from right side)
		// corresponds to number of detected value
		if (ind >= 0)	NVal = float(int(NVal) ^ (1 << ind));
		if(delim)	*delim = EnumDelims[1];		// restore 'vals' string
		else	break;							// no delimiter: last or single value
	}
	return ind < 0 ? PrintWrong(vals, ind==-2 ? "wrong delimiter in value" : strEmpty) : 0;
}

// Checks and sets pair option value
//	@vals: pair of values as C string or NULL if value isn't set
//	return: 0 if success, 1 if wrong values
int Options::Option::SetPair(const char* vals, bool isInt)
{
	if(!vals)	return 0;	// value isn't stated and it's optional (that was checked before)
	const char* delim = strchr(vals, EnumDelims[2]);	// a pointer to the delimiter ';' in vals

	if(!delim)	return PrintWrong(vals, "missed '" + string(1, EnumDelims[2]) + "' delimiter in value");
	PairVals& lim = *((PairVals*)SVal);

	if(delim != vals) {		// first value is set
		if(!IsValidFloat(vals, isInt, true)
		|| SetTriedFloat(float(atof(vals)), lim.Values(PairVals::MIN).first, lim.Values(PairVals::MAX).first))
			return 1;
		((pairVal*)SVal)->first = NVal;		// set first PairVals element
	}
	if(*(delim+1)) {			// second value is set
		if(!IsValidFloat(delim+1, isInt)
		|| SetTriedFloat(float(atof(delim+1)), lim.Values(PairVals::MIN).second, lim.Values(PairVals::MAX).second))
			return 1;
		((pairVal*)SVal)->second = NVal;		// set first PairVals element
	}
	return 0;
}

// Returns option name and value followed by message
string Options::Option::ToStr(bool prVal) const
{
	string res(optTitle);
	res += NameToStr(true);
	if(prVal)	res += sSPACE + string(SVal);
	return res;
}

// Prints option in full or short way.
//	@descr: if true, prints in full way: signature, description (marks as Required if needed), default value,
//	otherwise signature only
void Options::Option::Print(bool descr) const
{
	if(Sign.Is(fHidden))	return;
	USHORT	len = 0;		// first len is used as counter of printed chars
	bool	fixValType = ValType==tENUM || ValType==tCOMB;
	char*	buffer;

	if(descr)	cout << SPACE, len++;	// full way: double blank first
	// *** signature
	{
	string name = NameToStr(false);
	cout << SPACE << name;
	len += short(1 + name.length());
	}
	// *** option value type
	cout << SPACE, len++;
	if(IsValEsc())	cout << '[', len++;
	if(fixValType)	len += PrintEnumVals();		// print enum values
	else if((buffer = const_cast<char*>(TypeNames[ValType])) != NULL)
		cout << buffer,							// print value type
		len += short(strlen(buffer));
	if(IsValEsc())	cout << ']', len++;
	
	// *** description
	if(!descr)	return;
	// align description 
	short cnt = OPT_DESCF_TSHIFT - len / 8;	// OPT_DESCF_TSHIFT=3 * 8: right boundary of descriptions
	if(cnt <= 0)	cnt = 1;
	if(len + cnt*8 >= (OPT_DESCF_TSHIFT+1)*8)	// too long option,
		cnt = OPT_DESCF_TSHIFT, cout << LF;	// description on the next line
	for(int i=0; i<cnt; i++)	cout << TAB;

	// print description
	cnt = 0;	// use as external enum counter
	len = USHORT(strlen(Descr));	// from now len is used as length of description string
	buffer = new char[len+1];
	PrintSubLine(buffer, Descr, strchr(Descr, LF), fixValType ? (const char**)SVal : NULL, &cnt);
	delete [] buffer;
	if(AddDescr) {
		if(Descr[len-1] != LF)	cout << SPACE;
		cout << AddDescr;
	}
	if(Sign.Is(fOblig))	cout << " Required";
	else if(ValType >= tHELP)	cout << " and exit";

	// print default value
	if(ValRequired() && NVal != NO_VAL) {
		switch(ValType) {
			case tENUM:
			case tCOMB:
				if(NVal >= MinNVal)	// do not print default if it is never set by user
					PRINT_IN_PRTHS(((char**)SVal)
						[int(NVal)-int(MinNVal)] );	// offset by min enum value
				break;
			case tPR_INT:
			case tPR_FL:	
				PRINT_IN_PRTHS(PairValsToStr((pairVal*)SVal)); break;
			case tCHAR:	PRINT_IN_PRTHS(char(NVal)); break;
			default:	PRINT_IN_PRTHS(NVal);
		}
	}
	else if(SVal != NULL)	PRINT_IN_PRTHS(ValType == tENUM ? "NONE" : SVal);
	cout << LF;
}

// Prints enum or combi values
//	return: number of printed symbols
BYTE Options::Option::PrintEnumVals() const
{
	if( !ValRequired() )	return 0;

	char** vals = (char**)SVal;			// array of val images
	BYTE len = BYTE(strlen(vals[0]));	// number of printed chars
	BYTE cnt = BYTE(MaxNVal);			// number of values
	if(MinNVal)							// is range of values limited from below?
		cnt -= BYTE(MinNVal) - 1;
	cout << '<' << vals[0];			// first val image from array of val images
	for(BYTE i=1; i<cnt; i++) {
		cout << EnumDelims[ValType-tENUM] << vals[i];
		len += BYTE(strlen(vals[i]));
	}
	cout << '>';
	return len + cnt + 1;	// here cnt denotes the number of printed delimiters
}

// Performs a case-insensitive search of given string value among enum values.
//	@val: input value as string
//	return: index of finded value in enum,
//	or -1 if the value is not present in enum,
//	or -2 if the wrong delimiter is encountered
int Options::Option::GetEnumInd (const char* val)
{
	int i = 0;

	// check for delimiter
	for(char c = *val; c; c = *(val + ++i))
		if(!isalpha(c))	return -2;
	// detect value
	for(i=0; i<MaxNVal; i++)
		if( !_stricmp(val, ((const char**)SVal)[i]) )
			return i;
	return -1;
}

// Ouptuts option with error message to cerr
//	@val: value or NULL
//	@msg: error message about value
//	@return: always 1
int Options::Option::PrintWrong(const char* val, const string& msg) const
{
	cerr << ToStr(false) << SepSCl << (msg == strEmpty ? "wrong " + sValue : msg);
	if(val) cerr << SPACE << val;
	cerr << LF;
	return 1;
}

#ifdef DEBUG
void Options::Option::Print() const
{
	cout << Str << TAB;
	if(strlen(Str) < 4)	cout << TAB;
	Sign.Print(); cout << LF;
}
#endif

// Prints Usage params
void Options::Usage::Print(Option* opts) const
{
	if(Opt != NO_DEF)	// output option value
		opts[Opt].Print(false);
	else if(Par) {		// output parameter
		if(IsParOblig)	cout << SPACE << Par;
		else			PRINT_IN_PRTHS(Par);
		if(ParDescr)	// output parameter description
			cout << "\n  " << Par << ": " << ParDescr;
	}
	cout << endl;
}

// Check obligatory options and output message about first absent obligatory option.
//	return: -1 if some of obligatory options does not exists, otherwise 1
int Options::CheckObligs()
{
	for(int i=0; i<OptCount; i++)
		if( List[i].CheckOblig() < 0 )	return -1;
	return 1;
}

// Set option [with value] or splitted short options
//	@opt: option without HYPHEN
//	@val: value of option
//	@nextItem: next token after opt and val, or NULL
//	@argInd: the current index in argc; increased by 1 if value is accepted
//	return: 0 if success, 1 otherwise
int Options::SetOption (char* opt, char* val, char* nextItem, int& argInd)
{
	int i, res = 0;
	const bool isWord = opt[0]==HPH;
	char* sopt = opt += isWord;	// sliced option

	// parse possibly united short options
	do {
		for(i=0; i<OptCount; i++)
			if ((res = List[i].SetVal(sopt, isWord, val, nextItem, argInd)) > 0)
				return res;
			else if(!res)	break;
		if(res<0)	return PrintAmbigOpt(sopt, isWord, "unknown", opt);
	}
	while(!isWord && *++sopt);	// next slice of united short options
	return 0;
}

// Returns true if long option opt is defined
bool Options::Find(const char* opt)
{
	for(int i=0; i<OptCount; i++)
		if( List[i].Str && !strcmp(List[i].Str, opt) )
			return true;
	return false; 
}

// Ouptuts ambiguous option with error message to cerr
//	@opt: option
//	@isWord: true if option is a word
//	@headMsg: message at the beginning
//	@inOpt: initial option (in case of ambiguous composite)
//	return: always 1
int Options::PrintAmbigOpt(const char* opt, bool isWord, const char* headMsg, const char* inOpt)
{
	cerr << headMsg << SPACE << optTitle << HPH;
	if(isWord)	cerr << HPH << opt;
	else {
		cerr << *opt;
		if(inOpt && *(inOpt+1))
			cerr << " in " << HPH << inOpt;	// print only non-signle-char splitted short option
	}
	if(inOpt && Find(inOpt))
		cerr << ". Do you mean " << HPH << HPH << inOpt << '?';
	cerr << LF;
	return 1;
}

// Prints version
//	return: always 1
int	Options::PrintVersion()
{
	cout<< Product::Version
#ifndef _NO_ZLIB
		<< "\tzlib " << ZLIB_VERSION
#endif
		<< endl;
	return 1;
}

int Options::PrintSummary(bool prTitle)
{
	if(prTitle)	cout << Product::Title << ": ";
	cout << Product::Descr << LF;
	return 1;
}

// Prints 'usage' information
//	@title: if true prints title before information
//	return: 1 if title is settinf to true, 0 otherwise
int Options::PrintUsage (bool title)
{
	BYTE i, k;
	if(title)		PrintSummary(true), cout << LF;
	
	// output 'Usage' section
	cout << "Usage:";
	for(k=0; k<UsageCount; k++) {
		cout << TAB << Product::Title;	PRINT_IN_PRTHS("options");
		// output available options
		for(i=0; i<OptCount; i++)
			List[i].PrintOblig();
		// input parameters
		Usages[k].Print(List);
	}
	cout << endl;

	// output options section
	cout << "Options:\n";
	for(k=0; k<GroupCount; k++) {
		if(OptGroups[k])	cout << OptGroups[k] << ":\n";
		for(i=0; i<OptCount; i++)
			List[i].PrintGroup(k);
	}
	return int(title);
}

// Returns command line.
//	@argc: count of main() parameters
//	@argv: array of main() parameters
string const Options::CommandLine(int argc, char* argv[])
{
	ostringstream ss;
	int i = 0;
	argc--;
	while (i < argc)	ss << argv[i++] << SPACE;
	ss << argv[i];
	return ss.str();
}

// Parses and checks main() parameters and their values. Output message if some of them is wrong.
//	@argc: count of main() parameters
//	@argv: array of main() parameters
//	@obligPar: name of required application parameter or NULL if not required
//	return: index of first parameter (not option) in argv[],
//	or argc if it is absent,
//	or negative if tokenize complets wrong
int Options::Parse(int argc, char* argv[], const char* obligPar)
{
	if (argc < 2)	{ Options::PrintUsage(true); return -1; }		// output tip	
	int i, res = 1;
	char *token, *nextToken;	// option or parameter

	for (i = 1; i < argc; i++) {		// argv[0] is the path to the program
		token = argv[i];
		nextToken = argv[i+1];
		if( *token != HPH )	{			// token is not an option
			if( i < argc-1 				// not a last token
			&& *nextToken == HPH )		// next token is an option
				cerr << token << ": neither option nor parameter" << LF,	res = -1;
			break;
		}
		if( SetOption(
			token + 1,								// option without HYPHEN
			i + 1 <= argc ? nextToken : NULL,		// option's value
			i + 2 > argc ? NULL : argv[i + 2],		// next item or NULL
			i										// current index in argc
			) )
		{	res = -1; break; }
	}
	if( res > 0 )	res = CheckObligs();			// check required options
	if( res > 0 && obligPar && i == argc ) {
		cerr << Missing << obligPar << LF;
		res = -1;
	}
	return i * res;
}

// Return string value by index: if value is not oblig and is not specified, than defName with given extention
const string Options::GetFileName(int indOpt, const char* defName, const string& ext)
{
	Option opt = List[indOpt];
	if(opt.SVal)
		return string(opt.SVal) + ext;
	if(opt.IsValEsc() && opt.Sign.Is(fTrimmed))
		return FS::ShortFileName(FS::FileNameWithoutExt(defName)) + ext;
	return strEmpty;
}

#ifdef DEBUG
void Options::Print()
{
	for(int i=0; i<OptCount; i++)
	{
		cout << setw(2) << i << "  ";
		List[i].Print();
	}
}
#endif

/************************ end of class Options ************************/

/************************ class Err ************************/

//const char* Err::TREAT_BED_EXT = "after extension";	// clarifying message in the bed stretch operation

const char* Err::_msgs [] = {
/* NONE */		"WARNING",
/* MISSED */	"missing",
/* F_NONE */	"no such file",
/* D_NONE */	"no such folder",
/* FD_NONE */	"no such file or folder",
/* F_MEM */		"memory exceeded",
/* F_OPEN */	"could not open",
/* F_CLOSE */	"could not close",
/* F_READ */	"could not read",
/* F_EMPTY */	"empty",
/* F_BIGLINE */	"buffer is less than length of line",
/* FZ_MEM */	"not enough internal gzip buffer",
/* FZ_OPEN */	"wrong reading mode READ_ANY for zipped file",
/* FZ_BUILD */	"this build does not support zipped files",
/* F_WRITE */	"could not write",
///* F_FORMAT */	"wrong format",
#ifndef _FQSTATN
/* TF_FIELD */	"number of fields is less than required",
/* TF_EMPTY */	"no records",
#endif
/* EMPTY */		""
};

const char* Err::FailOpenOFile = "could not open output file";

// Initializes _outText by cstring contained message kind of "<sender>: <txt> <specTxt>".
void Err::set_message(const char* sender, const char *txt, const char *specTxt)
{
	size_t senderLen = sender != NULL ? strlen(sender) : 0;
	size_t outLen = senderLen + strlen(txt) + strlen(SepSCl) + 2;
	if (specTxt)	outLen += strlen(specTxt) + 1;
	_outText = new char[outLen];		// will be free in destructor
	memset(_outText, cNULL, outLen);
	if (sender) {
		//*_outText = SPACE;
		if (senderLen)	strcat(_outText, sender);
		strcat(_outText, SepCl);
	}
	strcat(_outText, txt);
	if (specTxt) {
		if (*specTxt != ':')		strcat(_outText, sSPACE);
		strcat(_outText, specTxt);
	}
}

// Returns string containing file name and issue number.
//	@issName: name of issue
//	@issNumb: number of issue
//	@fName: file name
//const string Err::IssueNumbToStr(const string& issName, ULONG issNumb, const string& fName)
//{
//	string res = fName;
//	if(fName != strEmpty)	res += SepSCl;
//	return res + issName + SPACE + NSTR(issNumb);
//}

// Gets message "no @fileName.@ext[.gz] files in this directory"
const string Err::MsgNoFiles (const string & fileName, const string ext)
{
	return string("no " + fileName + ext + "[" + ZipFileExt + "] files in this directory");
}


Err::Err(const Err & src)
{
	_code = src._code;
	int size = strlen(src.what()) + 1;
	_outText = new char[size];
	strcpy(_outText, src._outText);
	//_specifyText = src._specifyText;
}

// Throws exception or outputs Err message.
//	@throwExc: if true then throws exception, otherwise outputs Err message
//	if true then carriage return after Err message
void Err::Throw(bool throwExc, bool eol) {
	if (throwExc)
		throw* this;
	else {
		dout << what();
		if(eol)		dout << LF;
		fflush(stdout);
	}
}

// Outputs warning
//	@prefix: output ": " before "WARNING"
//	@eol: if true then carriage return after Err message
void Err::Warning(bool eol, bool prefix)
{
	if (prefix)	dout << SepCl;
	dout << _msgs[ErrWARNING];
	if(*what() != ':') dout << SepCl;	// check if sender is not recorded
	dout << what();
	if(eol)	dout << LF;
	fflush(stdout);
}

/************************ end of class Err ************************/

/************************ class FileSystem ************************/

// Returns true if file system's object exists
//	@name: object's name
//	@st_mode: object's system mode
bool FS::IsExist(const char* name, int st_mode)
{
	struct_stat64 st;
	int res, len = strlen(name) - 1;

	if(name[len] == SLASH) {
		string sname = string(name, len);
		res = _stat64(sname.c_str(), &st);
	}
	else res = _stat64(name, &st);
	return (!res && st.st_mode & st_mode);
	//return ( !_stat64(name, &st) && st.st_mode & st_mode );
}

// Checks if file system's object doesn't exist
//	@name: object's name
//	@st_mode: object's system mode
//	@throwExcept: if true throws exception,
//	otherwise outputs Err message as warning without LF
//	@ecode: error's code
//	return: true if file or directory doesn't exist
bool FS::CheckExist	(const char* name,  int st_mode, bool throwExcept, Err::eCode ecode)
{
	if(IsExist(name, st_mode))	return false;
	Err(ecode, name).Throw(throwExcept);
	return true;
}

// Searches through a file name for the any extention ('/' or '\' insensible).
//	return: the index of the DOT matched extention; otherwise npos
size_t FS::GetLastExtPos(const string &fname) {
	size_t pos = fname.find_last_of(DOT);
	return ( pos!=string::npos &&( pos==0 || ( pos==1 && fname[0]==DOT ))) ?	// ./name || ../name
		string::npos : pos;
}

bool SearchExt(const string& fname, const string& ext, bool isZip, bool composite)
{
	size_t pos = fname.find(ext);
	if(pos == string::npos)		return false;
	if(!composite)				return fname.size() - pos == ext.size();
	return fname.size() - pos - isZip * ZipFileExt.length() == ext.size();
}

// Returns true if file name has specified extension ignoring zip extension. Case insensitive
bool FS::HasCaseInsExt(const string &fname, const string &ext, bool knownZip, bool composite)
{
	bool res = SearchExt(fname, ext, knownZip, composite);
	if(!res) {		// try to case insensitive search
		string str(fname);
		string substr(ext);
		transform(str.begin(), str.end(), str.begin(), ::tolower);
		transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
		res = SearchExt(str, substr, knownZip, composite);
	}
	return res;
}


// Gets size of file or -1 if file doesn't exist
LLONG FS::Size (const char* fname)
{
	struct_stat64 st;
	return _stat64(fname, &st) == -1 ? -1 : st.st_size;
}
//// alternative implementation in C++ style;
// needs #include <fstream>
//std::ifstream::pos_type filesize(const char* fname)
//{
//    std::ifstream ifs(fname, std::ifstream::ate | std::ifstream::binary);
//    return ifs.tellg(); 
//}

// Gets real size of zipped file  or -1 if file cannot open; limited by UINT
LLONG FS::UncomressSize	(const char* fname)
{
	FILE *file = fopen(fname, "rb");	// "read+binary"
	if( file == NULL )		return -1;
	BYTE sz[] = {0,0,0,0};
	_fseeki64(file, -4, SEEK_END);
	fread(sz, 1, 4, file);
	fclose(file);
	return (sz[3] << 3*8) + (sz[2] << 2*8) + (sz[1] << 8) + sz[0];
}

// Throws exsception if file or directory doesn't exist
//	@name: name of file or directory
//	@ext: file extention; if set, check for file first
//	@throwExcept: if true throws exception,
//	otherwise outputs Err message as warning without LF
//	return: true if file or directory doesn't exist
bool FS::CheckFileDirExist(const char* name, const string & ext, bool throwExcept) 
{
	return HasExt(name, ext) ?
		CheckFileExist(name, throwExcept) :
		CheckFileDirExist(name, throwExcept);
}

// Returns a pointer to the file name checked if file exist, otherwise throws exception
//	@optsVal: Options char* value
//	return: pointer to the checked file name
const char* FS::CheckedFileDirName	(const char* name)
{
	if(!IsFileDirExist(name))	Err(Err::FD_NONE, name).Throw();
	return name;
}

// Returns a pointer to the file name checked if file exist, otherwise throws exception
//	@name: pointer to the file name
//	return: pointer to the checked file name
const char* FS::CheckedFileName	(const char* name)
{
	if(name && !IsFileExist(name))	Err(Err::F_NONE, name).Throw();
	return name;
}

// Returns a pointer to the path checked if it exist, otherwise throws exception
//	@opt: Options value
//	return: pointer to the checked path
const char* FS::CheckedDirName	(int opt)
{
	const char* name = Options::GetSVal(opt);
	if(name && !IsDirExist(name))	Err(Err::D_NONE, name).Throw();
	return name;
}


// Returns string containing real file extension (without zip extention).
//	@fname: pointer to the file name
//	return: string containing real file extension or empty string if no real extention
string const FS::GetExt(const char* fname) {
	const char * pdot = strrchr(fname, DOT);
	if( !pdot )		return strEmpty;
	if( strcmp(pdot, ZipFileExt.c_str()) )
		return string(pdot+1);				// no zip extention
	const char * pprevdot = pdot - 1;
	for(; pprevdot >= fname; pprevdot--)	// find previous DOT
		if( *pprevdot == DOT )	
			break;
	return pprevdot+1 == fname ? strEmpty : string(pprevdot+1, pdot-pprevdot-1);
}

// Returns true if file name is short (without path)
bool FS::IsShortFileName(const string& fname)
{
#ifdef OS_Windows
	if(fname.find_first_of(REAL_SLASH) != string::npos)	return false;
	else
#endif
	return fname.find_first_of(SLASH) == string::npos;
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

// Returns the name of last subdirectory
//	@name: long dir name
string const FS::LastDirName (const string& name)
{
	size_t pos = name.find_last_of(SLASH);
	if(pos == string::npos)		return name;
	if(++pos == name.length()) {				// ended by slash?
		pos = name.find_last_of(SLASH, pos - 2) + 1;
		return name.substr(pos, name.length() - 1 - pos);
	}
	return name.substr(pos);
}

// Returns the name ended by slash without checking
string const FS::MakePath(const string& name)
{
#ifdef OS_Windows
	if(name.find(REAL_SLASH) != string::npos) {
		string tmp(name + SLASH);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp;
	}
#endif
	return name[name.length()-1] == SLASH ? name : name + SLASH;
}

#if !defined _WIGREG && !defined _FQSTATN
// Fills external vector of strings by file's names found in given directory
// Implementation depends of OS.
//	@files: external vector of strings that should be filled by file's names
//	@dirName: name of directory
//	@ext: file's extention as a choosing filter
//	@all: true if all files with given extention should be placed into external vector,
//	otherwise only one (any)
//	return: true if files with given extention are found
bool FS::GetFiles	(vector<string>& files, const string& dirName, const string& ext, bool all)
{
#ifdef OS_Windows
	int count = 0;
	string fileTempl = FS::MakePath(dirName) + '*' + ext;
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
		do	files.emplace_back(ffd.cFileName);
		while (FindNextFile(hFind, &ffd));
	}
	else
		files.emplace_back(ffd.cFileName);
	FindClose(hFind);
	return true;
#else
	struct dirent *entry;
	DIR *dir = opendir(dirName.c_str());	// doesn't need to check because of appl options inspection
	string name;

	// fill files
	files.reserve(Chrom::Count);
	while( entry = readdir(dir) ) {
		name = string(entry->d_name);
		if(HasExt(name, ext, false)) {
			files.push_back(name);
			if( !all )	break;
		}
	}
	closedir (dir);
	return files.size() > 0;
#endif	// OS_Windows
}
#endif	// _WIGREG, _FQSTATN
/************************ end of class FileSystem ************************/

/************************  class TimerBasic ************************/

// Prints elapsed time
//	@elapsed: elapsed time in seconds
//	@watch: true if time should be printed as a stopwatch (with decimal places and without empty minutes)
//	@isLF: if true then ended output by LF
void PrintTime(long elapsed, bool watch, bool parentheses, bool isLF)
{
	int hours = elapsed/60;
	int mins = hours%60;
	hours /= 60;
	if(parentheses)	dout << '(';
	dout << setfill('0') << right;		// right couse it may be chanched by previuos output
	if(hours)	dout << setw(2) << hours << COLON;			// hours
	if(mins || !watch)	dout << setw(2) << mins << COLON;	// mins
	// secs
	if(watch)	dout << setw(5) << fixed << setprecision(2) << (elapsed - mins*60);
	else		dout << setw(2) << long(elapsed)%60;
	if(parentheses)	dout << ')';
	if(isLF)	dout << LF, fflush(stdout);
}

bool	TimerBasic::Enabled = false;

// Stops enabled timer and return elapsed wall time in seconds
long TimerBasic::GetElapsed() const
{
	time_t stopTime;
	time( &stopTime );
	return (long)difftime(stopTime, _startTime);
}

// Prints elapsed time interval
//	@elapsed: elapsed time in seconds 
//	@title: string printed before time output
//	@parentheses: if true then output time in parentheses
//	@isLF: if true then ended output by LF
void TimerBasic::Print(long elapsed, const char *title, bool parentheses, bool isLF)
{
	if(title)	dout << title;
	PrintTime(elapsed, false, parentheses, isLF);
}

/************************  end ofclass TimerBasic ************************/

/************************  class Timer ************************/
clock_t	Timer::_StartCPUClock;

// Stops enabled timer and prints elapsed time
//	@offset: space before time output
//	@isLF: if true then ended output by LF
void Timer::Stop(int offset, bool parentheses, bool isLF)
{
	if(_enabled) {
		if(offset)	dout << setw(offset) << SPACE;
		PrintTime(GetElapsed(), NULL, parentheses, isLF);
	}
}

/************************  end of class Timer ************************/
#ifdef _TEST
/************************  class Stopwatch ************************/

// Stops Stopwatch
//	@title: if not empty, and if instance was launched, output sum wall time with title
//	'const' to apply to constant objects
void Stopwatch::Stop(const string title) const
{
	if(!_isStarted)		return;
	_sumTime += GetElapsed();
	if(title!=strEmpty)	PrintTime(_sumTime, (title + sSPACE).c_str(), false, true);
}

/************************  end of class Stopwatch ************************/
#endif	// _TEST
/************************  class StopwatchCPU ************************/

// Stops StopwatchCPU
//	@title: string printed before time output
//	@print: if true time should be printed
//	@isLF: if true then ended output by LF
void StopwatchCPU::Stop(const char* title, bool print, bool isLF)
{
	_sumclock += clock() - _clock;
	if(print) {
		if(title)	dout << title;
		PrintTime(_sumclock/CLOCKS_PER_SEC, false, false, isLF);
	}
}

/************************  end of class StopwatchCPU ************************/

#ifdef _MULTITHREAD
/************************  class Mutex ************************/

bool	Mutex::_active;		// true if the mutex really should work
mutex	Mutex::_mutexes[int(Mutex::eType::NONE)];

/************************  end of class Mutex ************************/
#endif	// _MULTITHREAD

/************************ class Chrom ************************/

const char*		Chrom::Abbr = "chr";
const BYTE		Chrom::MaxAbbrNameLength = BYTE(strlen(Abbr)) + MaxMarkLength;
#ifndef _FQSTATN
const char*		Chrom::Marks = "XYM";
const string	Chrom::UndefName = "UNDEF";
const string	Chrom::Short = "chrom";
const string	Chrom::sTitle = "chromosome";
const BYTE		Chrom::MaxShortNameLength = BYTE(Short.length()) + MaxMarkLength;
const BYTE		Chrom::MaxNamedPosLength = BYTE(strlen(Abbr)) + MaxMarkLength + CHRLEN_CAPAC + 1;
	  BYTE		Chrom::CustomOpt;

chrid Chrom::cID = UnID;		// user-defined chrom ID
chrid Chrom::firstHeteroID = 0;	// first heterosome (X,Y) ID

// Gets somatic (letter) chrom's ID by mark without control, or undefined ID
chrid Chrom::HeteroID	(const char cMark)
{
	if(!IsRelativeID())		return cMark;		// absolute ID
	for(size_t i=0; i<strlen(Marks); i++)
		if(cMark == Marks[i])
			return chrid(firstHeteroID + i);	// relative ID
	return UnID;
}

// Gets chrom ID by case insensitive mark
//	firstHeteroID should be initialized!
chrid Chrom::CaseInsID	(const char* cMark)
{
	if(isdigit(*cMark))	{				// autosome
		chrid id = atoi(cMark) - 1;
		return id < firstHeteroID ? id : UnID;
	}
	return CaseInsHeteroID(*cMark);		// heterosome
}

// Returns a pointer to the first occurrence of C substring. Recurcive.
//	@str: C string to find in
//	@templ: C string to be located
//	@templLen: length of templ (extern because of avoiding recursive recalculate)
const char* SubStr(const char* str, const char* templ, int templLen)
{
	str = strchr(str, *templ);
	if(str)
		for(short i=1; i<templLen; i++)
			if( *++str != templ[i] )
				return SubStr(str, templ, templLen);
	return str;
}

// Returns the length of prefix
//	return: length of substring before short chromosome's name, or -1 if short name is not finded
short Chrom::PrefixLength(const char* cName)
{
	// search from the beginning of the cName because
	// it simpler to find the first digit for multidigit mark

	// start search from the occurrence of 'chr'
	const char* str = SubStr(cName, Abbr, strlen(Abbr));
	if( str )
		for(; *str; str++)
			if( isdigit(*str) || isupper(*str) )
				return str - cName;
	return -1;
}

// Gets chrom's ID by name without control of case insensitivity and undefined ID
//	@cName: chrom's name
//  @prefixLen: length of name prefix
chrid Chrom::ID(const char* cName, size_t prefixLen)
{
	return isdigit(*(cName+=prefixLen)) ? atoi(cName)-1 : HeteroID(*cName);
}

// Validates chrom name and returns chrom ID
//	@cName: string of arbitrary length, starting with chrom's name
//  @prefixLen: length of prefix before mark
chrid Chrom::ValidateID(const char* cName, size_t prefixLen)
{
	if(!cName)					return UnID;
	cName += prefixLen;							// skip prefix
	for(int i=1; i<=MaxMarkLength; i++)
		if(cName[i] == USCORE)	return UnID;	// exclude chroms with '_'

	if(isdigit(*cName))	{						// autosome
		chrid id = atoi(cName);
		if(/*IsRelativeID() && */id > firstHeteroID)	firstHeteroID = id;
		return id - 1;
	}
	return CaseInsHeteroID(*cName);				// heterosome
}

// Sets custom chrom ID with control
//	@prColon: if true then print ": " before exception message
void Chrom::SetCustomID(bool prColon)
{ 
	const char* mark = Options::GetSVal(CustomOpt);		// null if no chrom is set by user
	if (mark && (cID = CaseInsID(mark)) == UnID) {
		ostringstream ss;
		ss << "no such " << sTitle << " in this genome";
		if (prColon)	ss << SepCl;
		ss << Options::OptionToStr(CustomOpt, true);
		Err(ss.str()).Throw();
	}
}

// Sets number of 'custom chrom' progr option
//	@absIDNumb: true if absolute ID numbering discipline is applied
void Chrom::SetCustomOption(int opt/*, bool absIDNumb*/)
{
	CustomOpt = opt;
	if( !(firstHeteroID/* = !absIDNumb*/) )
		cID = ValidateID(Options::GetSVal(opt));	// apply absolute numeration discipline
}

inline const string AutosomeToStr(chrid cid) {	return to_string(cid + 1); }

// Returns mark by ID
const string Chrom::Mark(chrid cid)
{
	if(cid == UnID)		return UndefName;
	if(IsRelativeID())
		return cid < firstHeteroID ? AutosomeToStr(cid) :
			(cid > firstHeteroID + 2) ? UndefName : string(1, Marks[cid - firstHeteroID]);
	return cid < '9' ? AutosomeToStr(cid) : to_string(cid);
}

// Locate chrom mark in string.
//	@str: string checked for chrom number
//	return: pointer to the chrom number in str, or a null pointer if Chrom::Abbr is not part of str.
const char* Chrom::FindMark(const char* str)
{
	const char* substr = strstr(str, Abbr);
	return substr ? (substr + strlen(Abbr)) : NULL;
}

// Gets chrom's abbreviation name by ID
//	@numbSep: if true then separate chrom's number
string Chrom::AbbrName(chrid cid, bool numbSep)
{ 
	return Abbr + (numbSep ? sSPACE : strEmpty) + Mark(cid);
}

#endif	// _FQSTATN

/************************ end of class Chrom ************************/

/************************ struct Region ************************/

// Extends Region with chrom length control.
// If extended Region starts from negative, or ends after chrom length, it is fitted.
//	@extLen: extension length in both directions
//	@cLen: chrom length; if 0 then no check
void Region::Extend(chrlen extLen, chrlen cLen)
{
	Start -= extLen > Start ? Start : extLen;
	End += extLen;
	if (cLen && End > cLen)	End = cLen;
}

/************************ end of struct Region ************************/

/************************ class Regions ************************/

// Geta total length of regions.
//chrlen Regions::Length () const
//{
//	chrlen len = 0;
//	for(Iter it=_regions.begin(); it<_regions.end(); it++)
//		len += it->Length();
//	return len;
//}

//Regions& Regions::operator=(const Regions& rgn)
//{
//	_regions = rgn._regions;
//	return *this;
//}

#if defined _READDENS || defined _BIOCC

// Returns an iterator referring to the past-the-end element, where end is external
//	@curr_it: region's const iterator, from which the search is started
//	@end: external pre-defined end coordinate
Regions::Iter Regions::ExtEnd(Iter curr_it, chrlen end) const
{
	Iter it = curr_it;
	for (; it != _regions.end(); it++)
		if (it->End > end)	break;
	return it;
}

// Initializes this instance by intersection of two Regions.
//	Typically for that purpose is used Interval Tree,
//	but this implementation uses Regions ordering and is simpler.
void Regions::FillOverlap(const Regions& regn1, const Regions& regn2)
{
	chrlen start = 0, end = 0, start1, start2, end1, end2;
	Iter it1 = regn1._regions.begin();
	Iter it2 = regn2._regions.begin();
	Reserve(max(regn1.Count(), regn2.Count()));
	for (; it1 != regn1._regions.end(); it1++) {
		start1 = it1->Start;	end1 = it1->End;
		start2 = it2->Start;	end2 = it2->End;
		if (start1 < start2) {
			if (end1 > start2) {
				start = start2;
				if (end1 > end2) { end = end2;	it2++; it1--; }
				else				end = end1;
			}
		}
		else
			if (start1 >= end2) { it2++; it1--; }
			else {
				start = max(start1, start2);
				if (end1 > end2) { end = end2;	it2++; it1--; }
				else				end = end1;
			}
		if (end) {
			Add(start, end);
			if (it2 == regn2._regions.end())		return;
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
void Regions::FillInvert(const Regions& regn, chrlen maxEnd)
{
	Region rgn;
	Iter it = regn._regions.begin();

	Reserve(regn.Count() + 1);
	for (; it != regn._regions.end(); it++) {
		rgn.End = it->Start - 1;
		Add(rgn);
		rgn.Start = it->End + 1;
	}
	if (rgn.Start < (maxEnd + 1)) {
		rgn.End = maxEnd;
		Add(rgn);
	}
}

#endif	// _READDENS, _BIOCC

#ifdef DEBUG
void Regions::Print() const
{
	int i = 0;
	cout << "Regions:\n";
	for (Iter it = _regions.begin(); it < _regions.end(); it++)
		cout <<
		//setw(2) << 
		++i << COLON << TAB <<
		//setw(9) << 
		it->Start << TAB << it->End << endl;
}
#endif
/************************ end of class Regions ************************/

/************************  MemStatus ************************/

//bool	MemStatus::_enable;
//LLONG	MemStatus::_startVolume;
//
//LLONG MemStatus::getAvailMemory()
//{
//#ifdef OS_Windows
//	MEMORYSTATUSEX status;
//	status.dwLength = sizeof(status);
//	GlobalMemoryStatusEx(&status);
//	//return (size_t)status.ullAvailPhys;
//	return (size_t)status.ullAvailVirtual;
//#else
//	long pages = sysconf(_SC_PHYS_PAGES);
//	long page_size = sysconf(_SC_PAGE_SIZE);
//	return pages * page_size;
//#endif
//}
//
//void MemStatus::StartObserve(bool enable)
//{
//	if(_enable = enable)	
//		_startVolume = getAvailMemory(); 
//}
//
//void MemStatus::StopObserve()
//{
//	if(_enable) {
//		size_t rem = getAvailMemory();
//		cout << "Unallocated memory: " << (_startVolume - rem) << LF;
//	}
//}
/************************  end of MemStatus ************************/