#pragma once

#include "def.h"
#include <string>
#include <iostream>	
#include <sstream>		// for NSTR( x )
#include <iomanip>		// setprecision(), setw()
#include <vector>
#include <sys/stat.h>	// struct stat

#ifdef __unix__
	#include <unistd.h>
	#include <stdlib.h>
	#include <memory.h>
	#include <limits.h>
	#include <math.h>
	#include <dirent.h>
	#include <stdio.h>
#ifdef _MULTITHREAD
	#include <pthread.h>
	#define InterlockedExchangeAdd	__sync_fetch_and_add
	#define InterlockedIncrement(p)	__sync_add_and_fetch(p, 1)
#endif
	#include <string.h>		// strerror_r()
	#include <stdexcept>	// throw std exceptions
	//#include <errno.h>	// pthread_join errcodes

	typedef unsigned char	BYTE;
	typedef	unsigned short	USHORT;
	typedef	unsigned int	UINT;
	typedef unsigned long	ULONG;
	typedef unsigned long long ULLONG;
	typedef long long LLONG;
	typedef struct stat struct_stat64;
	//typedef off64_t __int64;
	typedef void*	retThreadValType;
	#define retThreadValTrue	(void*)1

	#define rand()	rand_r(&_seed)
	#define _stricmp strcasecmp	// case-sensitive comparison
	#define _fseeki64 fseeko64
	#define _ftelli64 ftello64
	#define _fileno	fileno
	#define _gcvt	gcvt
	#define _stat64 stat
	//#define _itoa	itoa
	#define fopen	fopen64
	#define TRUE	1
	#define FALSE	0
	//#define SLASH '/'	// standard Linux path separator
	//#define F_READ_MODE O_RDONLY
#elif defined _WIN32
	#define OS_Windows
	#include <windows.h>
#ifdef _MULTITHREAD
	#include <process.h>	    // _beginthread, _endthread
	#define pthread_t HANDLE
	#define pthread_mutex_t CRITICAL_SECTION
#endif
	typedef unsigned __int64 ULLONG;
	typedef __int64 LLONG;
	typedef struct __stat64 struct_stat64;
	typedef UINT	retThreadValType;
	#define retThreadValTrue	1

	#define atol _atoi64
	//#define SLASH '\\'		// standard Windows path separator
	//#define REAL_SLASH '/'	// is permitted in Windows too
#endif

#ifndef _NO_ZLIB
	#ifdef OS_Windows
		#define ZLIB_WINAPI
	#endif
	#include "zlib.h"
	#if ZLIB_VERNUM >= 0x1240
		#define ZLIB_NEW
	#else	
		#define ZLIB_OLD
	#endif
#endif	// _NO_ZLIB

// specific types
typedef char threadnumb;	// type number of thread
typedef BYTE chrid;			// type number of chromosome
typedef BYTE readlen;		// type length of Read
typedef UINT chrlen;		// type length of chromosome
typedef ULONG genlen;		// type length of genome
typedef float readscr;		// type score of Read

#define CHRLEN_CAPAC	10	// capacity of max chrom length;
							// may be count by DigitsCount() every time,
							// but is the same if chrlen defined as int or long,
							// so is defined as static value
#define	CHRLEN_UNDEF	UINT_MAX	// undefined length of chromosome
#define retThreadValFalse	0

#define TAB	'\t'
#define cN	'N'
#define HPH	'-'
#define USCORE	'_'
#define BLANK	' '
#define sBLANK	" "
#define sBACK	"\b"
#define QUOT	'\''
#define DOT		'.'
//#define HASH	'#'
#define COLON	':'
#define PERS	'%'

#define EOL '\n'		// 10
#define CR	'\r'		// 13
#define MSGSEP_TAB ":\t"	// message tab separator between subject and explanation
#define MSGSEP_BLANK ": "	// message blank separator between subject and explanation

using namespace std;
static const char* WARNING = "NOTICE: ";

static const string ZipFileExt = ".gz";
static const string StrEmpty = "";
static const char* Total = "total";
static const char* Version = "version";
#ifdef _ISCHIP
static const char* GroupParSep = ";  ";		// option print: separator of option values in a group
static const char* MsgFileAbsent = " is not exist. Generate...";
#endif	// _ISCHIP
#ifndef _WIGREG
static const char* Template = "template";
#endif
static const char* MsgDone = " done\t";

/*** COMMON MACROS & FUNCTION ***/

#define vUNDEF	-1	// undefined value

// Digital to STRing
// Returns value's string representation. http://rootdirectory.de/wiki/NSTR()
// Instead of std::to_string(x) [C++11] because of compatibility
#define NSTR( x ) static_cast<ostringstream & >( ostringstream() << dec << (x) ).str()

// Returns two int value's separated by delimiter string representation.
#define NNSTR( x, delim, y ) static_cast<ostringstream & >( ostringstream() << dec << int(x) << delim << int(y) ).str()

// Byte to STRing
#define BSTR(x) static_cast<ostringstream & >( ostringstream() << dec << int(x) ).str()

// Float to STRing
// Returns float's string representation.
//	@p: precision (number of digits after separator)
#define FSTR( x, p ) static_cast<ostringstream & >( ostringstream() << dec << setprecision(p)<< x ).str()

// Returns number of value's symbols, including minus & separator
#define SCNT( x ) NSTR(x).length()

// Returns number of float's symbols, including minus and dot
//	@p: precision (number of digits after separator)
#define SFCNT( x, p ) FSTR(x, p).length()

// Gets number of digist in a integral value
//	@val: integral value
//	return: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCount (LLONG val);

// Returns percent of @part relatively @total
inline float	Percent(ULLONG part, ULLONG total) { 
	return total ? 100.f * part / total : 0.f;
}

// Returns string represents the percent of part relatively total
//	@percent: value of percent
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default
//	@parentheses: if true parenthesize the value (default)
string	sPercent(float percent, BYTE precision=0, BYTE fieldWith=0, bool parentheses=true);

// Returns string represents the percent of part relatively total
//	@part: value represents desired %
//	@total: value represents 100%
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default
//	@parentheses: if true parenthesize the value (default)
inline string	sPercent(ULLONG part, ULLONG total,
	BYTE precision=0, BYTE fieldWith=0, bool parentheses=true) {
		return sPercent(Percent(part, total), precision, fieldWith, parentheses);
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
chrlen AlignPos(chrlen pos, BYTE res, BYTE relative);

#endif

// Gets available system memory
//size_t	getAvailSystemMemory();

//#ifdef OS_Windows
//string	Wchar_tToString(const wchar_t* wchar);
//wchar_t* StringToWchar_t(const string &str, wchar_t* wchar);
//#endif

/*** end of COMMON FUNCTION ***/

static struct Product
/*
 * struct 'Product' keeps Product info
 */
{
	static const string	Title;
	static const string	Version;
	static const string	Descr;

	// Gets length of product name (title plus version).
	static inline BYTE NameLength() {
		return BYTE(Title.length());		// + strlen(Version) + 1);
	}

	// Gets title plus version.
	static inline const string& Name() {
		return Title;	// + string(1, HPH) + string(Version);
	}
} product;

static class Options
/*
 * Class 'Options' implements main() options and parameters treatment.
 */
{
public:
	static const char* Booleans [];		// boolean values

	// Gets pointer to the C string contained 'ON' or 'OFF'
	//	@val: true of false
	static inline const char* GetBoolean(bool val) { return Options::Booleans[int(val)]; }

	// Gets pointer to the C string contained 'ON' or 'OFF'
	//	@i: option index
	static inline const char* GetBoolean(int i) { return Options::Booleans[Options::GetBVal(i)]; }

	// Prints 'usage' information
	//	@title: if true prints title before information
	//	return: 1 if title is settinf to true, 0 otherwise
	static int PrintUsage (bool title);
	
	// Returns command line.
	//	@argc: count of main() parameters
	//	@argv: array of main() parameters
	static const string CommandLine(int argc, char* argv[]);

	// Parses and checks main() parameters and their values.
	//	Output message if some of them is wrong.
	//	@argc: count of main() pearmeters
	//	@argv: array of main() pearmeters
	//	@obligPar: name of required application parameter or NULL if not required
	//	return: index of first parameter (not option) in argv[], argc if it is absent,
	//	negative if  tokenize complets wrong
	static int Tokenize(int argc, char* argv[], const char* obligPar=NULL);
	
	//static void GetOpt(int i);

	// Get double value by index
	static inline double GetDVal(int i)	{ return _Options[i].NVal; }
	// Get float value by index
	static inline float GetFVal	(int i)	{ return float(_Options[i].NVal); }
	// Get string value by index
	static inline char* GetSVal (int i)	{ return _Options[i].SVal; }
	//Get booling value by index
	static inline bool GetBVal	(int i)	{ return _Options[i].NVal != 0; }
	// Get UINT value by index
	static inline UINT GetIVal	(int i)	{ return UINT(_Options[i].NVal); }
	// True if maximal enunm value is setting
	static inline bool IsMaxEnum(int i)	{ 
		return _Options[i].NVal == _Options[i].MaxNVal - 1;
	}
	// Get maximal permissible numeric value by index
	static inline UINT GetMaxIVal(int i){ return UINT(_Options[i].MaxNVal); }
	// Get descriptor by index
	//static inline const char* GetDescr	(int i) { return _Options[i].Descr; }

private:
	// types of option values
	// do not forget to support the correlation with Options::_TypeNames []
	enum eValType	{ tUNDEF, tNAME, tCHAR, tINT, tFLOAT, tLONG, tENUM, tCOMB, tHELP, tVERS };
	
	struct Signs {
	private:
		BYTE	signs;

	public:
		static const BYTE Oblig		= 0x1;	// obligatory option sign
		static const BYTE Recogn	= 0x2;	// recognized option sign 
		//static const BYTE Printed	= 0x4;	// printed option mask

		inline Signs(int x)	{ signs = (BYTE)x; }	// to initialize obligatory in main()
		// Returns true if given sign is set
		bool inline Is(BYTE mark)		{ return signs & mark; }
		// Sets given sign
		void inline MarkAs(BYTE mark)	{ signs |= mark; }
	};

	// structure 'Option' keeps options attributes.
	struct Option {
		const char	Char;		// option - character
		const char*	Str;		// option - string
			  Signs	Sign;		// initialize obligatory in main():
								// 1 if option is obligatory, 0 otherwise;
		//const bool	_ValRequired;// true if option's value is required
		const eValType ValType;	// type of value
		const BYTE	OptGroup;	// option's category
			  double NVal;		// default or established numeric or enum value
		const double MinNVal;	// minimal permissible value;
								// for enum should be a first defined value
		const double MaxNVal;	// maximal permissible value or count of enum values
			  char*	SVal;		// string value or pointer to enum values array;
								// in last case should be cast to char**.
								// If enum option hase ValRequired==false,
								// this pointer dousn't use and may be NULL
		const char*	Descr;		// tip string

		static const char EnumDelims[];	// specifies delimiter for enum [0] and combi [1] values

		// returns true if option value is required
		bool inline ValRequired() { return MinNVal != vUNDEF; }

		// Sets option value.
		//	@isWord: true if option is a word, false if option is a char
		//	@opt: option without HYPHENs
		//	@val: value of option
		//	@isNextOpt: true if next parameter is option
		//	return: 0 if success, -1 if not found, 1 1 if option or value is wrong
		int	SetVal(bool isWord, char* opt, char* val, bool isNextOpt);

		// Check option for obligatory.
		//	return: -1 if option is obligatory but not stated, otherwise 1
		int CheckOblig();

		// Prints option with double blank before if descr==true, single blank otherwise.
		void Print(bool descr);

	private:
		// Checks limits and set numerical value
		//	@val: numerical value
		//	@isWord: true if option enters by long name
		//	return: 1 if limits are exceeded, otherwise 0
		int SetTriedDigit(double val, bool isWord);

		// Checks and sets enum option value.
		//	@val: input value as C string
		//	return: true if success
		bool SetEnum(char* val);

		// Checks and sets enum option value.
		//	@val: input value as C string
		//	return: true if success
		bool SetComb(char* val);

		// Prints enum or combi values
		//	return: number of printed symbols
		BYTE PrintEnumVals();
		
		// Performs a case-insensitive search of given string value among enum values.
		//	@val: input value as string
		//	return: index of finded value in enum, or -1 if the value is not present in enum
		int GetEnumInd (char* val);
	};

	// structure 'Usage' is used to output some Usage variants in PrintUsage()
	struct Usage {
		const int		OptVal;	// enum optValue which string option should be printes in Usage
		const string	Text;	// text which string option should be printes in Usage
	};

	static const char*	_TypeNames [];	// names of option value types in help
	static const char*	_OptGroups [];	// names of option groups in help
	static	Option		_Options[];		// options. Option 'help' always shuld be the last one,
										// option 'version' - before last.
	static const Usage	_Usages[];		// content of 'Usage' variants in help
	static const BYTE	_OptCount,		// count of options
						_GroupCount,	// count of option groups in help
						_UsageCount;	// count of 'Usage' variants in help
	
	static inline string OptToStr	(bool isWord, const char* opt) {
		return string(1 + int(isWord), HPH) + opt;
	}
	
	// Check obligatory options and output message about first absent obligatory option.
	//	return: -1 if some of obligatory options does not exists, otherwise 1
	static int	CheckObligs();
	
	// Set value of option 
	//	@isWord: true if option is a word, false if option is a char
	//	@opt: option with HYPHENs
	//	@val: option's value
	//	@isNextOpt: true if next parameter is option
	//	@argIndex: the current index in argc; is increasing by 1 in case of required value
	//	Return: 0 if success, 1 otherwise
	static int	SetOption(bool isWord, char* opt, char* val, bool isNextOpt, int *argIndex);

	// Prints version
	//	return: always 1
	static int PrintVersion();

	// Ouptuts option with error message to cerr
	// Ouptuts option with error message to cerr
	//	@isWord: true if option is long
	//	@opt: option
	//	@val: value or NULL
	//	@msg: error message about value
	static int PrintWrongOpt(bool isWord, const char* opt, const char* val,
		const string msg = StrEmpty);

	// Ouptuts ambiguous option with error message to cerr
	//	@isWord: true if option is long
	//	@opt: option
	//	@headMsg: message at the beginning
	//	@tailMsg: message at the end or NULL
	static int PrintAmbigOpt(bool isWord, const char* opt, const char* headMsg, const char* tailMsg);

} options;
#define ErrWARNING	Err::NONE

class Err
/*
 * Class 'Error' implements error's (warning's) messages and treatment.
 */
{
public:
	enum eCode {
		NONE,
		P_MISSED,
		F_NON,
		FD_NON,
		F_MEM,
		F_OPEN,
		F_CLOSE,
		F_READ,
		F_WRITE,
		F_BIGLINE,
		F_NOREADREC,
		FZ_MEM,
		FZ_OPEN,
		FZ_BUILD,
#ifndef _FQSTATN
		TF_FIELD,
		TF_SPEC,
		TF_EMPTY,
		BP_BADEND,	// bed position: start is is equal or more than end
		BP_NEGPOS,	// bed position: negative position
		BP_EXCEED,	// bed position: exceeded chrom length
#ifdef _BEDR_EXT
		BR_RNAME,	// bed read: wrong read name format
#endif
		//BR_SIZE,
		FA_LONGLEN,
#endif
#if defined(_ISCHIP) || defined(_FQSTATN)
		FQ_HEADER,
		FQ_HEADER2,
#elif defined _DENPRO || defined _BIOCC
		ARR_OUTRANGE,
		SUM_EXCEED,
#endif
		EMPTY
	};
private:
	static const char* _msgs[];
	enum eCode	_code;			// error code
	char * _outText;			// output message
	//const char * _specifyText;		// aditional text to output message
	
	// Initializes _outText by C string contained message kind of
	// "<@sender>: <@text> <@specifyText>".
	void set_message(const char* sender, const char* text, const char* specifyText=NULL);

	//inline void set_message(const string& sender, const char* text, const char* specifyText=NULL) {
	//	set_message(sender==StrEmpty ? NULL : sender.c_str(), text, specifyText);
	//}

public:
	// Gets message "no @fileName.@fileExt[.gz] files in this directory"
	static const string MsgNoFiles (const string & fileName, const string fileExt)
	{
		return string("no " + fileName + fileExt +
			"[" + ZipFileExt + "] files in this directory");
	}

	// Code-attached constructor.
	//	@code: exception/warning message as code
	//	@sender: name of object who has generated exception/warning
	//	@specifyText: aditional text to specify exception/warning message
	inline Err(eCode code, const string& sender=StrEmpty, const char* specifyText=NULL): _code(code) {
		set_message(sender.c_str(), _msgs[code], specifyText);
	}
	// Code-attached constructor.
	//	@code: exception/warning message as code
	//	@sender: name of object who has generated exception/warning
	//	@specifyText: aditional text to specify exception/warning message
	inline Err(eCode code, const string& sender, const string& specifyText) : _code(code) {
		set_message(sender.c_str(), _msgs[code], specifyText.c_str());
	}
	// C-string-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const char* text, const char* sender=NULL) : _code(NONE) {
		set_message(sender, text);
	}
	// String-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const string& text, const char* sender=NULL) : _code(NONE) {
		set_message(sender, text.c_str());
	}
	// String-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const char* text, const string& sender) : _code(NONE) {
		set_message(sender.c_str(), text);
	}
	// String-attached constructor.
	//	@text: exception/warning message
	//	@sender: name of object who has generated exception/warning
	inline Err(const string& text, const string& sender) : _code(NONE) {
		set_message(sender.c_str(), text.c_str());
	}
	// copy constructor
	Err(const Err& src);

	inline ~Err() { if( _outText) delete [] _outText; }

	inline const char* what() const throw() { return _outText;	}

	inline eCode Code() const		{ return _code; }

	// Returns point to the aditional text to specify exception/warning message or NULL
	//inline const char* SpecifyText() const	{ return _specifyText; }

	// Throws exception or outputs warning.
	//	@throwException: if true throws exception, otherwise outputs Err message
	//	@endOfLine: if true print EOL after Err message
	//	return: false
	void	Throw(bool throwException = true, bool endOfLine = true);
	
	// Outputs warning with prefix "WARNING" and additional text, if it is setting.
	void	Warning	(string const & addText = StrEmpty);
};

static class FS
/*
 * Class 'FileSystem' implements common file routines 
 */
{
private:
	// Returns true if file system's object exists
	//	@name: object's name
	//	@st_mode: object's system mode
	static bool IsExist(const char* name, int st_mode);

	// Checks if file system's object doesn't exist
	//	@name: object's name
	//	@st_mode: object's system mode
	//	@throwExcept: if true throws excwption,
	//	@ecode: error's code
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	static bool CheckExist	(const char* name,  int st_mode, bool throwExcept, Err::eCode ecode);

	// Searches through a file name for the any extention (slash|back-slash insensible).
	//	@fname: file name
	//	return: the index of the DOT mathed extention; otherwise npos
	static size_t GetExtPos	(const string &fname);

public:
	// Gets size of file or -1 if file doesn't exist
	static LLONG Size 	(const char*);

	// Gets real size of zipped file  or -1 if file cannot open; limited by UINT
	static LLONG UncomressSize	(const char*);

	// Returns true if file exists
	inline static bool IsFileExist	 (const char* name) { return IsExist(name, S_IFREG); }
	
	// Returns true if directory exists
	inline static bool IsDirExist	 (const char* name) { return IsExist(name, S_IFDIR); }
	
	// Returns true if file or directory exists
	inline static bool IsFileDirExist(const char* name) { return IsExist(name, S_IFDIR|S_IFREG); }
	
	// Checks if file doesn't exist
	//	@name: name of file
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file doesn't exist
	inline static bool CheckFileExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFREG, throwExcept, Err::F_NON);
	}

	// Checks if file or directory doesn't exist
	//	@name: name of file or directory
	//	@throwExcept: if true throws excwption,
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	static inline bool CheckFileDirExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR|S_IFREG, throwExcept, Err::FD_NON);
	}

	// Throws exsception if file or directory doesn't exist
	//	@name: name of file or directory
	//	@ext: file extention; if set, check for file first
	//	@throwExcept: if true throws exception,
	//	otherwise outputs Err message as warning without EOL
	//	return: true if file or directory doesn't exist
	static bool CheckFileDirExist(const char* name, const string & ext, bool throwExcept) 
	{
		return HasExt(name, ext) ?
			CheckFileExist(name, throwExcept) :
			CheckFileDirExist(name, throwExcept);
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optsVal: Options char* value
	//	return: pointer to the checked file name
	static const char* CheckedFileDirName	(const char* name) {
		if( !IsFileDirExist(name) )	Err(Err::FD_NON, name).Throw();
		return name;
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optVal: Options value
	//	return: pointer to the checked file name
	static inline const char* CheckedFileDirName	(int optVal) {
		return CheckedFileDirName(Options::GetSVal(optVal));
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@name: pointer to the file name
	//	return: pointer to the checked file name
	static const char* CheckedFileName	(const char* name) {
		if( !IsFileExist(name) )	Err(Err::F_NON, name).Throw();
		return name;
	}

	// Returns a pointer to the file name checked if file exist, otherwise throws exception
	//	@optVal: Options value
	//	return: pointer to the checked file name
	static inline const char* CheckedFileName	(int optVal) {
		return CheckedFileName(Options::GetSVal(optVal));
	}

	// Returns true if file has any extension.
	//	@fname: file name
	static inline bool HasExt	(const string &fname) { 
		return GetExtPos(fname) != string::npos;
	}

	// Returns true if file has a specified  extension.
	//	@fname: file name
	//	@ext: extension includes dot symbol and can be composite
	static bool HasExt	(const string& fname, const string& ext);

	// Returns true if file has '.gz' extension
	static inline bool HasGzipExt(const string& fname) { return HasExt(fname,ZipFileExt); }

	// Returns string containing real file extension (without zip extention).
	//	@fname: pointer to the file name
	//	return: string containing real file extension or empty string if no real extention
	static string const GetExt(const char* fname);

	// Returns file name without extentiom (slash|back-slash insensible)
	static string const FileNameWithoutExt (const string& fname);

	// Returns short file name by long one
	//	@fname: long file name
	static string const ShortFileName (const string& fname);

	// Returns directory name by long file name
	//	@fname: long file name
	//	@addSlash: true if slash sould be added at the end
	static string const DirName (const string& fname, bool addSlash = false);

	// Returns the name of last subdirectory by long file name
	//	@fname: long file name
	static string const LastSubDirName (const string& fname);

	// Returns the name ended by slash without checking the name
	static string const MakePath(const string& name);

#ifndef _WIGREG
	// Fills external vector of strings by file's names found in given directory
	// Implementation depends of OS.
	//	@files: external vector of strings that should be filled by file's names
	//	@dirName: name of directory
	//	@fileExt: file's extention as a choosing filter
	//	@all: true if all files with given extention should be placed into external vector,
	//	otherwise only one (any)
	//	return: true if files with given extention are found
	static bool GetFiles	(vector<string>& files, const string& dirName, const string& fileExt, bool all = true);
#endif	// _WIGREG

//	static inline void	Delete		(const char* fname) {
//#ifdef OS_Windows
//		DeleteFile(fname);
//#else
//		unlink(fname);
//#endif
//	}

} fs;

class Timer
{
private:
	time_t _startTime;
	static clock_t	_StartCPUClock;

	// Prints elapsed time interval
	//	@title: string printed before time output
	//	@elapsed: elapsed time interval
	//	@parentheses: if true then output time in parentheses
	//	@isCarrgReturn: if true then ended output by EOL
	static void PrintElapsed(const char *title, long elapsed, bool parentheses, bool isCarriageReturn);

public:
	// True if timing is enabled
	static bool		Enabled;

	// Starts enabled CPU timer, if it is enabled
	static inline void StartCPU()	{ if( Enabled ) _StartCPUClock = clock(); }
	
	// Stops enabled CPU timer and print elapsed time
	//	@isCarrgReturn: if true then ended output by EOL
	static void StopCPU(bool isCarrgReturn=true);

	// Creates a new Timer and starts it if timing is enabled
	inline Timer()	{ Start(); }
	
	// Restarts timer, if timing is enabled
	inline void Start()				{ if( Enabled ) time( &_startTime ); }

	// Stops enabled timer and print elapsed time with title
	//	@title: string printed before time output
	//	@parentheses: if true then output time in parentheses
	//	@isCarrgReturn: if true then ended output by EOL
	void Stop(const char *title, bool parentheses, bool isCarrgReturn);
	
	// Stops enabled timer and print elapsed time
	//	@parentheses: if true then output time in parentheses
	//	@isCarrgReturn: if true then ended output by EOL
	inline void Stop(bool parentheses = false, bool isCarrgReturn = true)	{
		Stop(NULL, parentheses, isCarrgReturn); }
};


#ifdef _ISCHIP
	#define dout	cout
#elif defined _WIGREG
	#define dout	cerr
#else
class dostream
/*
 * 'dostream' duplicates outstream to stdout & file
 */
{
    std::ostream &first, &second;
public:
    dostream(std::ostream &f, std::ostream &s) : first(f), second(s) {}
	template <typename T> dostream& operator<< (T val) {
		first << val;
		second << val;
		return *this;
	}
};

extern ofstream outfile;	// file ostream duplicated stdout; inizialised by file in code
extern dostream dout;		// stream's duplicator

/*** end of class dostream ***/

#endif	// _DENPRO || _BIOCC || _FQSTATN

#ifdef _MULTITHREAD
static class Mutex
{
private:
	static const BYTE Count = 2;
	static pthread_mutex_t	_mutexes[];
public:
	enum eType { OUTPUT, WR_FILE };
	static void Init();
	static void Finalize();
	static void Lock(const eType type);
	static void Unlock(const eType type);
} mutex;

class Thread
{
private:
	pthread_t _thread;

public:
	Thread(retThreadValType(
		#ifdef OS_Windows
		__stdcall
		#endif
		*proc)( void*), void *arglist);

	inline ~Thread()
	{
	#ifdef OS_Windows
		CloseHandle(_thread);
	#endif
		//;
	}

	// Waits for a thread to terminate and detaches the thread.
	inline void WaitFor()
	{
#ifdef OS_Windows
		//DWORD res = WaitForSingleObject(_thread, INFINITE);
		//switch(res) {
		//	case WAIT_ABANDONED: cout << " WAIT_ABANDONED\n"; break;
		//	case WAIT_OBJECT_0:	cout << " WAIT_OBJECT_0\n"; break;
		//	case WAIT_TIMEOUT: cout << " WAIT_TIMEOUT\n"; break;
		//	case WAIT_FAILED: cout << " WAIT_FAILED\n"; break;
		//}
		WaitForSingleObject(_thread, INFINITE);
#else
		//int res = pthread_join(_thread, NULL);
		//switch(res) {
		//	case 0: cout << " WELL\n"; break;
		//	case EINVAL: cout << " EINVAL\n"; break;
		//	case ESRCH:	cout << " ESRCH\n"; break;
		//	case EDEADLK: cout << " EDEADLK\n"; break;
		//}
		pthread_join(_thread, NULL);
#endif
	}

// Terminates the calling thread. Does not used because of uninvoke destructors under Windows
//	void Thread::EndThread()
//	{
//#ifdef OS_Windows
//		_endthread();
//#else
//		pthread_exit(NULL);
//#endif
//	}
};
#endif	// _MULTITHREAD

#ifdef _BIOCC

static struct CCkey	// Represents correlation coeficient constants and methods
{
private:
	static const int Undef = 666;	// undefined coefficient

public:
	enum eCC {	// defines types of CC
		ccP = 0x1,	// Pearson correlation coefficient: first bit
		ccS = 0x2,	// signal correlation coefficient: second bit
	};

	// Returns true if parameter is signal coefficient
	inline static bool IsS(eCC ecc)		{ return static_cast<bool>(ecc&ccS); }
	// Returns true if parameter is Pearson coefficient
	inline static bool IsP(eCC ecc)		{ return static_cast<bool>(ecc&ccP); }
	// Returns true if parameter is both coefficients
	inline static bool IsBoth(eCC ecc)	{ return static_cast<bool>(ecc&ccS && ecc&ccP); }
	// Returns correlation coefficient
	inline static double CalcR(double v1, double v2, double cov) {
		return  !v1 && !v2 ? Undef : cov/(sqrt(v1) * sqrt(v2));	// not sqrt(v1*v2) 
	}
	inline static void Print(double val) {
		if( val == Undef )	dout << "UNDEF";
		else				dout << val;
	}
} cc;

typedef pair<double, double> pairDbl;

struct CC : pairDbl
/*
 * 'CC' represents a pair of correlation coeficients: Pearson(first) and signal(second),
 * and provides theirs output.
 */
{
private:
	static const int Empty = -2;	// uninitialised coefficient

public:
	inline CC() { first = second = Empty; }

	inline CC(double p, double s) { first = p; second = s; }

	inline bool operator < (const CC & ccres) const { 
		return first != Empty ? first < ccres.first : second < ccres.second; 
	}

	// Sets Pearson coefficient
	inline void SetP(double val)	{ first = val;	}
	// Sets signal coefficient; undefined value by default
	inline void SetS(double val=0)	{ second = val;	}

	// Returns true if even one value in pair is setting
	inline bool NotEmpty() const { return( first != Empty || second != Empty ); }

	// returns single value
	inline double GetSingleVal() const { return first != Empty ? first : second; }

	// set single negative value to absolute value
	void SetSingleAbsVal() { 
		if( first != Empty )	first = fabs(first);
		else					second = fabs(second);
	}

	void Print() const { 
		if( first != Empty )	{
			CCkey::Print(first);
			dout << TAB;
			if( SCNT(first) < 8 )	dout << TAB;
		}
		if( second != Empty )
			CCkey::Print(second);
		dout << EOL; 
	}
};

struct CCaggr
/*
 * CCaggregate represents accumulative values needed for calculation CC for a set of arrays
 */
{
private:
	static const UINT Undef = -1;	// uninitialised coefficient
	ULLONG check;					// needs for check-up only
	ULLONG VarS[3];								// Signal variances
	double VarP1, VarP2, CovP, Mean1, Mean2;	// Pearson variances

	// Increases S value with check-up for exceeding
	//	@i: index of value on VarS
	//	@val: value to increase
	//	@cID: chrom's ID, needed for exception only
	//	@return: false if all right 
	bool IncreaseSVar(BYTE i, ULLONG val, chrid cID);

public:
	CCaggr() {
		VarP1 = VarP2 = CovP = Mean1 = Mean2 = 0;
		//VarS[0] = VarS[1] = VarS[2] = 0;
		memset(VarS, 0, 3*sizeof(ULLONG));
	}
	// Sets predefined means
	inline void SetMeans(const pairDbl& means) {
		Mean1 = means.first;
		Mean2 = means.second;
	}

	// Increases P values in consideration of means
	//	@x: first value to increase
	//	@y: second value to increase
	void IncreasePVars(double x, double y);

	// Increases S values with check-up for exceeding
	//	@x: first value to increase
	//	@y: second value to increase
	//	@cID: chrom's ID, needed for exception only
	//	@return: false if all right 
	bool IncreaseSVars(ULLONG x, ULLONG y, chrid cID);

	// True if Signal coefficient is undefined (exceeding digital limit)
	inline bool IsUndefS() { return static_cast<bool>(VarS[0] == Undef); }

	// Calculates coefficients
	CC GetR() const {
		CC cc;
		if( VarP1 || VarP2 )	cc.SetP(CCkey::CalcR(VarP1, VarP2, CovP));
		if( VarS[0] == Undef )	cc.SetS(CCkey::CalcR(0, 0, 0));
		else if( VarS[0] || VarS[1] )	
			cc.SetS(CCkey::CalcR(double(VarS[0]), double(VarS[1]), double(VarS[2])));
		return cc;
	}
};

//static const string sOutOfRange = " is out of range ";	// used in Array only

#endif	// _BIOCC

template <typename T> class Array
{
private:
	ULONG	_len;
#ifdef _BIOCC
	mutable double _mean;
#endif // _BIOCC

protected:
	T	*_data;

public:
	// Reserves array capacity, initialized by 0.
	//	@len: capacity
	void Init(long len) {
		if( (_len=len) > 0 ) {	_data = new T[len]; memset(_data, 0, len*sizeof(T)); }
		else					_data = NULL;
	}

	// Creates an instance with capacity, initialized by 0.
	//	@len: capacity
	inline Array(long len=0)
#ifdef _BIOCC
	: _mean(0)
#endif // _BIOCC
	{ Init(len); }

	inline ~Array()		{ 
		//if( _data )	cout << long(_data[0]);
		//cout << "\tdestructor " << _len << EOL;
		if(_data)	delete [] _data; }
	
	inline bool Empty() const	{ return !_len; }

	inline ULONG Length() const	{ return _len; }

	inline T& operator[](const ULONG i)				{ return _data[i]; }

	inline const T& operator[](const ULONG i) const	{ return _data[i]; }

	Array& operator=(const Array& arr) {
		Reserve(arr._len);
		if( arr._len )	copy(arr._data, arr._data+arr._len, _data);
		return *this;
	}

	// Deletes previous data and reserves array capacity, initialized by 0.
	//	@len: capacity
	inline void Reserve(long len)	{ 
		if(_data)	delete [] _data;
		Init(len);
	}

	// CLears an instance: set all members to zero
	inline void	Clear() { if(_data)	memset(_data, 0, _len*sizeof(T)); }

#if defined _DENPRO || defined _BIOCC
	// Adds array to this instance.
	//	@src: added array
	//	@start: position from which array should be added
	//	@startSrc: start position in added array; start of added array by default
	//	@size: length of added region; whole added array by default
	inline void Concat(const Array& src, ULONG start, ULONG startSrc=0, ULONG size=0) {
		if( !size )	size = src._len - startSrc;
		if( start + size > _len )
			Err(Err::ARR_OUTRANGE, "Array.Concat()", NSTR(_len)).Throw();
			//Err("added arrays length " + NSTR(start + size) +
			//sOutOfRange + NSTR(_len), "Array.Concat()").Throw();
		memcpy(_data + start, src._data + startSrc, size*sizeof(T));
	}

	// Fills subarray by value
	//	@begin: first position of subarray
	//	@end: last position of subarray
	//	@val: value to fill
	void Fill(long begin, long end, T val) {
		if( end > _len )
			Err(Err::ARR_OUTRANGE, "Array.Fill()", NSTR(_len)).Throw();
		for(long i=begin; i<=end; i++)
			_data[i] = val;
	}
#endif	// _DENPRO, _BIOCC
#ifdef _BIOCC
	// Gets mean
	//double Mean() const {
	//	//if( !_mean ) {
	//	//	ULONG sum = 0;
	//	//	for (long i=0; i<_len; i++)
	//	//		sum += _data[i];
	//	//	_mean = double(sum) /_len;
	//	//	cout << "mean: " << _mean << EOL;
	//	//}
	//	if( !_mean )	cout << "empty mean!\n";
	//	return _mean; 
	//}
private:
	// Gets a pair of means of subarrays for this instance and given array
	pairDbl GetMeans(bool notGotThis, const Array & arr, long begin, long end) const
	{
		ULONG sum1=0, sum2=0;
		for (long i=begin; i<end; i++) {
			if(notGotThis)	sum1 += _data[i];
			sum2 += arr._data[i];
		}
		return pairDbl( notGotThis ? 
			double(sum1) / (end-begin) : 0,
			double(sum2) / (end-begin) );
	}

public:
	// Gets a pair of means of subarrays for this instance and given array,
	// and set _mean for each array in case of whole arrays
	pairDbl SetMeans(const Array& arr, long begin=0, long end=0) const
	{
		pairDbl means;

		if( end )		// get means for subarrays	
			means = GetMeans(true, arr, begin, end);
		else			// the whole arrays
			if(!arr._mean) {	// is arr._mean not setting?
				means = GetMeans(!_mean, arr, 0, _len);
				// initialize _mean for each array
				if(!_mean)		_mean = means.first;
				arr._mean = means.second;
			}
			else {				// both _means are setting: initialize return pair
				means.first = _mean;
				means.second = arr._mean;
			}

		return means;
	}

	// Returns maximal value in region
	T GetMaxVal(long begin=0, long end=0) const {
		T res = 0, val;
		if( !end )	end = _len;
		for (long i=begin; i<end; i++)
			if( (val=_data[i]) > res )		res = val;
		return res;
	}

	// Multiplys value in region to ratio
	void MultiplyVal(float ratio, long begin=0, long end=0) {
		if( ratio == 1 )		return;
		if( !end )	end = _len;
		for (long i=begin; i<end; i++)
			_data[i] = T(ratio * _data[i]);
	}

	// Multiplys all values in region to ratio synchronously for pair of arrays
	//static void SynchMultiplyVal(Array<T>& arr1, Array<T>& arr2, 
	//	float ratio1, float ratio2, long begin=0, long end=0)
	//{
	//	if( ratio1 == ratio2 == 1 )		return;
	//	if( !end )	end = arr1._len;	// both arrays have the same length
	//	//else if( begin >= end )
	//	//	Err("begin " + NSTR(begin) + " is more or equal end " + NSTR(end),
	//	//		"Array.MultiplyVal").Throw();
	//	for (long i=begin; i<end; i++) {
	//		if( ratio1 != 1 )	arr1._data[i] = T(ratio1 * arr1._data[i]);
	//		if( ratio2 != 1 )	arr2._data[i] = T(ratio2 * arr2._data[i]);
	//	}
	//}

	// Calculates correlation coefficients
	//	@cID: chrom's ID. Needed for exception only
	//	@ecc: identifier what coefficient to calculate
	//	@arr: second array to compare
	//	@begin: low boundary of calculated range
	//	@end: high boundary of calculated range
	//	@return: pair of coefficients
	CC const GetR (chrid cID, CCkey::eCC ecc, const Array & arr, long begin=0, long end=0) 
	{
		CCaggr ccaggr;
		if( CCkey::IsP(ecc) )
			ccaggr.SetMeans( SetMeans(arr, begin, end) );
		AccumVars(cID, ecc, ccaggr, arr, begin, end);
		return ccaggr.GetR();
	}
	
	// Accumulates variances and covariances for calculating CC
	//	@cID: chrom's ID. Needed for exception only
	//	@ecc: identifier what coefficient to calculate
	//	@ccaggr: accumulative aggregate
	//	@arr: second array to compare
	//	@begin: low boundary of calculated range
	//	@end: high boundary of calculated range
	void AccumVars (chrid cID, CCkey::eCC ecc, CCaggr & ccaggr, const Array & arr, long begin=0, long end=0) const
	{
		if( !end )	end = _len;
		else if( begin >= end )
			Err("begin " + NSTR(begin) + " is equal or more than the end " + NSTR(end),
				"Array.GetR").Throw();
		long i;
		if( CCkey::IsS(ecc) )	// signal
			for (i=begin; i<end; i++)
				if( ccaggr.IncreaseSVars(_data[i], arr._data[i], cID) )
					break;				// if exceeded
		if( CCkey::IsP(ecc) )	// Pearson
			for (i=begin; i<end; i++)
				ccaggr.IncreasePVars(_data[i], arr._data[i]);
	}

	// Returns mean value in a range
	//double GetMean(long begin=0, long end=0) const {
	//	ULONG sum = 0;
	//	for (long i=begin; i<end; i++)	sum += _data[i];
	//	return double(sum) / (end - begin);
	//}

	//// Returns sum of values in a range
	//ULONG GetSum(long begin=0, long end=0) const {
	//	ULONG sum = 0;
	//	if( !end )	end = _len;
	//	for (long i=begin; i<end; i++)	sum += _data[i];
	//	return sum;
	//}
#endif	// _BIOCC
#ifdef DEBUG
	void Print(long begin=0, long end=0) const {
		if( !end )	end = _len;
		for (long i=begin; i<end; i++)
			dout << i << TAB << _data[i] << EOL;
	}
#endif	// _DEBUG
};

#if defined _DENPRO || defined _BIOCC
typedef Array<chrlen> arrchrlen;
#endif	// _DENPRO || _BIOCC
#ifdef _BIOCC
typedef Array<BYTE> arrbmapval;
#endif	// _BIOCC

static class Chrom
/*
 * Class 'Chrom' establishes correspondence between inner chromosome's ID and it's name.
 * Terms:
 * Short chromosome's name: 1, ..., Y.
 * Abbreviation chromosome's name: chr1, ..., chrY
 * Title chromosome's name: chromosome X, ..., chromosome Y
 * Prefix: substring before short chromosome's name: chr, chromosome.
 * Chromosome's ID is an integer that corresponds to short name for digital names,
 * and is a letter's ASCII code for literal names (X, Y).
 */
{
	static const char	X = 'X';
	static const char	Y = 'Y';
	static const char*	UndefName;

	static chrid _cID;	// user-defined chrom ID

public:
	static const char	M = 'M';
	static const chrid	UnID = 0;				// Undefined ID
	static const chrid	Count = 24;				// Count of chromosomes by default
	static const BYTE	MaxShortNameLength = 2;	// Maximal length of short chrom's name
	static const BYTE	MaxAbbrNameLength;		// Maximal length of abbreviation chrom's name
	static const BYTE	MaxNamedPosLength;	// Maximal length of named chrom's position 'chrX:12345'
	static const char*	Abbr;				// Chromosome abbreviation
	static const string	Title;				// Chromosome title

	// Gets chromosome's ID stated by user
	static inline chrid StatedID()	{ return _cID; }

	// Returns true if user stated all chromosomes
	static inline bool StatedAll()	{ return _cID == UnID; }

	// Sets chromosome's ID with validation.
	//	@cID: chromosome's ID
	//	return: true if the check has passed, otherwise do not set and print message to cerr
	static bool SetStatedID(chrid cID);

	// Sets chromosome's ID stated by user with validation.
	//	@cName: chromosome's name
	//	return: true if the check has passed, otherwise do not set and print message to cerr
	static inline bool SetStatedID(const char* cName) { return SetStatedID(ID(cName)); }

	// Gets chromosome's ID by name
	//	@cName: any name ('file.fa' etc) or short name by default
	//  @prefixLen: length of prefix or 0 by default
	static chrid ID(const char* cName, size_t prefixLen=0);

	//// Checks chrom abbreviation
	////	@cName: chrom abbreviation name
	////	return: true if name is wrong
	//inline static bool CheckAbbr(const char* cName) { 
	//	return strncmp(cName, Abbr, strlen(Abbr)) != 0;
	//}

	// Locate chrom number in string.
	//	@str: string checked for chrom number
	//	return: pointer to the chrom number in str,
	//	or a null pointer if Chrom::Abbr is not part of str.
	static const char* FindNumb(const char* str);

	// Gets chromosome's short name by ID.
	inline static string Name(chrid cID) {
		return cID != UnID ? (cID < M ? (BSTR(cID)) : string(1, cID)) : UndefName;
	}

	// Gets chromosome's abbreviation name (with prefix 'chr') by its ID.
	inline static string AbbrName(chrid cID)	{ return Abbr + Name(cID); }

	// Gets chromosome's title name (with prefix 'chromosome') or plural
	//	@cID: chromosome's ID or UnID if plural
	inline static string TitleName(chrid cID = UnID)	{ 
		return Title + (cID==UnID ? "s" : (sBLANK + Name(cID)));
	}
	// Gets chromosome's ID by abbreviation name (with prefix 'chr')
	inline static chrid IDbyAbbrName(const char* cAbbrName) { return ID(cAbbrName, strlen(Abbr)); }

	// Gets chromosome's ID by string short (without prefix) name
	inline static chrid ID(const string & cName){ return ID(cName.c_str(), 0); }

	// Gets chromosome's ID by long name
	inline static chrid IDbyLongName(const char* cLongName) { 
		return ID(cLongName, PrefixLength(cLongName)); }

	// Gets the length of short chromosome's name by ID
	static inline BYTE NameLength(chrid cID) {
		return cID != UnID ? (cID < M ? (cID < 10 ? 1 : 2) : 1) : strlen(UndefName);
	}

	//// Gets the length of abbreviation chromosome's name by ID
	//static inline BYTE AbbrNameLength(chrid cID) {
	//	return NameLength(cID) + strlen(Abbr);
	//}

	// Gets the length of prefix (substring before short chromosome's name)
	//	@cLongName: long chromosome's name
	//	return: length of substring before short chromosome's name
	//	or -1 if short name is not finded
	static short PrefixLength(const char* cLongName);

	inline static const string Absent(chrid cID, const string & what) {
		return AbbrName(cID) + " is absent in " + what + " file: skipped";
	}

	//static char* LongToShortName(char* longName) ;
} chrom;

#if !defined _WIGREG && !defined _FQSTATN

struct Read
/*
 * 'Read' represents Read (with name and score in case of _BEDR_EXT) as item
 */
{
#if defined _ISCHIP || defined _BEDR_EXT
	enum rNameType {	// defines types of Read's name
		nmUndef = 0,	// undefined name: never stated by user
		nmNumb	= 1,	// a unique number stated as name
		nmPos	= 2		// a start position stated as name
	};

	static char	SeqQuality;						// the quality values for the sequence (ASCII)
	static BYTE	MapQuality;						// the mapping quality
	static const char	NmPosDelimiter;			// delimiter between two positions in pair
	static const char*	NmNumbDelimiter;		// delimiter between chrom name and number
	static const BYTE	NmDelimiterShift = 2;	// shift to pass ':N' for nmPos type
	static const string	NmSuffMate1,			// suffix of Read's name on pair-end first mate
						NmSuffMate2;			// suffix of Read's name on pair-end second mate
	static BYTE	OutNameLength;			// Maximum length of Read name in output file
#endif	// _ISCHIP or _BEDR_EXT
	static	readlen	Len;						// length of Read
	static const char	Strand[2];				// strand markers: [0] - positive, [1] - negative

#ifdef _ISCHIP

private:
	static short	LimitN;				// maximal permitted number of 'N' in Read or vUNDEF if all
	static ULONG	Count;				// counter of total writed Reads
	static rNameType	NameType;			// type of name of Read in output files
	static const char ToUp;				// shift beween lowercase and uppercase characters
	static const char Complements[];	// template for complementing Read

public:
	static ULONG	MaxCount;	// up limit of writes Reads
	static const char*	NmDelimiter;	// delimiter between chrom name & value: ":N" or ":"

	static void Init(readlen rLen, rNameType name, char seqQual, BYTE mapQual, short limN, ULONG maxCnt);

	// Gets true if start position is stated as name
	//static inline bool IsPositionName () { return NameType == nmPos; }
	// Gets true if Read name keeps its number
	static inline bool IsNameAsNumber () { return NameType == nmNumb; }

	// Gets the common part of Read name in output files':'
	static inline const string& Name () { return Product::Title; }

	// Increments counter of total writed Reads thread-safely.
	//	return: true if limit is exceeded.
	static inline bool IncrementCounter() {	return InterlockedIncrement(&Count) >= MaxCount; }

	// Copies complemented Read.
	static void CopyComplement(char* dst, const char* src);

	// Checks Read for number of 'N'
	//	@read: checked Read
	//	return: -1 if Read is NULL, 0 if N limit is exceeded, 1 if success
	static int CheckNLimit(const char* frag);

	// Prints Read values - parameters.
	static void Print();

#else
	chrlen	Pos;		// Read's actual start position

#ifdef _BEDR_EXT
	chrid	InitCID;	// initial chrom - owner
	size_t	Num;		// Read's number or initial start position
	readscr	Score;		// Read's score

	inline Read(chrlen pos, chrid cid, size_t num, readscr score)
		: Pos(pos), InitCID(cid), Num(num), Score(score) {}
#else
	inline Read(chrlen pos) : Pos(pos) {}

#endif
	// Compares two Reads by position. For sorting a container.
	static inline bool CompareByStartPos(const Read& r1, const Read& r2) {
		return r1.Pos < r2.Pos;
	}

	//static inline bool CompareByNum(const Read& r1, const Read& r2) {
	//	return r1.Num < r2.Num;
	//}

	void inline Print() const {
		cout << Pos 
#ifdef _BEDR_EXT
			 << TAB << Num << TAB << int(Score)
#endif
			 << EOL;
	}
#endif
};

#endif
