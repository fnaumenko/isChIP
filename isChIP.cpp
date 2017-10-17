/*
	isChIP (In-Silico ChIP) is a fast realistic ChIP-seq simulator.
	The model is based on the real protocol of ChIP-seq 
	and was developed by Dr. T.Subkhankuliva (subkhankul@hotmail.com)

	The real protocol is simulated by repeating the basic cycle.
	Each basic cycle corresponds to single cell simulation, and consists of the next phases:
	• “shearing” the chromatin, or random cutting the reference genome in fragments;
	• “extraction” of the fragments overlapping with the binding events;
	• amplification the selected fragments if required in number of cycles;
	• “loss” of selected fragments according to desired percentage;
	• “contamination” with background fragments;
	• size selection: selection of fragments fitted to desirable size;
	• sequencing of the fragments from positive and negative strands.
	Simulated binding events are specified by single optional parameter, called template.

	Path to UCSC reference genome: ftp://hgdownload.soe.ucsc.edu/goldenPath/

	Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)

	This program is free software. It is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	GNU General Public License for more details.
 */

#include "isChIP.h"
#include "Imitator.h"

using namespace std;

const string Product::Title = "isChIP";
const string Product::Version = "1.0";
const string Product::Descr = "ChIP-seq simulator";

const string DefFileName [] = { "mTest", "mInput", "mRegular" };
const string OutFileTip = "location of output files or existing directory\n[test mode: " + 
	DefFileName[TEST] + ".*, " + 
	"control mode: " + DefFileName[CONTROL] + ".*, " +
	"regular mode: " + DefFileName[REGULAR] + ".*]";  

enum eOptGroup	{ oINPUT, oTREAT, oFRAG, oDISTR, oREAD, oOUTPUT, oOTHER };	// oOTHER should be the last 
const BYTE	Options::_GroupCount = oOTHER + 1;

const char* Options::_OptGroups [] = {
	"Input", "Processing", "Fragment", "Fragment's size distribution", "Reads", "Output", "Other"
};

// --smode option
const char* smodes	[] = { "SE", "PE" };		// corresponds to OutFile::eMode
// -read-name option: types of Read name notations
const char* rnames	[] = { "NMB", "POS" };		// corresponds to Read::rNameType; nmUndef is hidden
// --format option: format notations
const char* formats	[] = { "FQ", "BED", "SAM" };// corresponds to OutFile::eFormat	
// --verbose option: verbose notations
const char* verbs	[] = { "CRIT", "RES", "RT", "PAR", "DBG" };

const char* ForTest = "For the test mode only";

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::_Options [] = {
	{ 'a', "amplif",	0,	tINT,	oTREAT, 1, 1, 2000, NULL, "coefficient of amplification", NULL },
	{ 'b', "bg-level",	0,	tFLOAT,	oTREAT, 1, 0, 100, NULL,
	"number of selected fragments outside the features,\nin percent of foreground.", ForTest },
	{ HPH, "fg-level",	0,	tFLOAT,	oTREAT, 100, 0, 100, NULL,
	"in test mode the number of selected fragments within the features,\nin percent;\nin control mode the number of selected fragments, in percent", NULL },
	{ 'n', "cells",		0,	tLONG,	oTREAT, 1, 1, 1e7, NULL, "number of cells", NULL },
	{ 'g', "gen",		1,	tNAME,	oINPUT, vUNDEF, 0, 0, NULL,
	"reference genome library or single nucleotide sequence.", NULL },
	{ 'c', Chrom::Abbr,	0,	tCHAR,	oTREAT, vUNDEF, 0, 0, NULL,
	"generate output for the specified chromosome only", NULL },
	{ HPH, "frag-len",	0,	tINT,	oFRAG, 200, 50, 400, NULL, "average size of selected fragments", NULL },
	{ HPH, "frag-dev",	0,	tINT,	oFRAG, 20, 0, 200, NULL, "deviation of selected fragments", NULL },
	{ HPH, "bg-all",	0,	tENUM,	oTREAT, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off generation background for all chromosomes.\n", ForTest },
	{ HPH, "bind-len",	0,	tINT,	oTREAT, 1, 1, 100, NULL, "minimum binding length.", ForTest },
	{ HPH, "flat-len",	0,	tINT,	oTREAT, 0, 0, 200, NULL, "boundary flattening length.", ForTest },
	{ HPH, "mean",		0,	tINT,	oDISTR, 200, 0, 1500, NULL,
	"expectation of the based normal distribution", NULL },
	{ HPH, "sigma",		0,	tINT,	oDISTR, 200, 1, 700, NULL,
	"standard deviation of the based normal distribution", NULL },
	{ HPH, "ln-factor",	0,	tINT,	oDISTR, 500, 10, 900, NULL,
	"power multiplication factor in lognormal distribution", NULL },
	{ HPH, "ln-term",	0,	tFLOAT,	oDISTR, 5.1, 1, 7, NULL,
	"power summand in lognormal distribution", NULL },
	{ HPH, "let-N",		0,	tENUM,	oTREAT, FALSE, vUNDEF, 2, NULL,
	"include the ambiguous reference characters (N) on the beginning\nand on the end of chromosome", NULL },
	{ 'r', "read-len",	0,	tINT,	oREAD, 50, 20, 200, NULL, "length of output read", NULL },
	{ HPH, "read-name",	0,	tENUM,	oREAD, Read::nmPos, Read::nmNumb, Read::nmPos, (char*)rnames,
	"name of read in output files includes:\n? - read`s unique number within chromosome\n? - read`s true start position", NULL },
	{ HPH,"read-Nlimit",0,	tINT,	oREAD, vUNDEF, 0, 100, NULL,
	"maximum permitted number of ambiguous characters (N) in read [--read-len]", NULL },
	{ HPH,"reads-limit",0,	tLONG,	oREAD, 2e8, 1e5, (float)ULONG_MAX, NULL,
	"maximum permitted number of total written reads", NULL },
	{ HPH, "fq-qual",	0,	tCHAR,	oREAD, '~', '!', '~', NULL,
	"the quality values for the read in FQ output", NULL },
	{ HPH, "map-qual",	0,	tINT,	oREAD, 42, 0, 42, NULL,
	"the mapping quality for the read in SAM output", NULL },
	{ HPH, "sz-sel",	0,	tENUM,	oFRAG, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off fragment's size selection", NULL },
	{ HPH,"sz-sel-sigma",0,	tINT,	oFRAG, 20, 1, 100, NULL,
	"standard deviation of the size selection normal distribution", NULL },
//{ HPH, "stdout",	0, tBOOL,	FALSE, 0, 0, NULL, "Write output to the 'standard output' filehandle instead of 'standart error'" },
	{ 'm', "smode",		0,	tENUM,	oTREAT, OutFile::mSE, OutFile::mSE, OutFile::mEmpty, (char*)smodes,
	"sequencing mode: ? - single end, ? - paired end", NULL },
	{ HPH,"strand-admix",0, tENUM, oTREAT, FALSE, 0, 2, (char*)Options::Booleans,
	"turn on/off opposite strand admixture at the bound\nof binding site.", ForTest },
	{ HPH,"ts-uni",	0,	tENUM,	oTREAT, FALSE, vUNDEF, 2, NULL, "uniform template score.", ForTest },
	{ 'p',"threads",0,	tINT,	oTREAT, 1, 1, 50, NULL, "number of threads", NULL },
	{ HPH,"fix",	0,	tENUM,	oTREAT, FALSE, vUNDEF, 2, NULL,
	"fix random emission to get repetitive results", NULL },
	{ 'R',"regular",0,	tINT,	oTREAT, vUNDEF, 1, 400, NULL,
	"regular mode: write each read on starting position\nincreased by stated shift", NULL },
	{ 'f',"format",	0,	tCOMB,	oOUTPUT, OutFile::ofFQ, OutFile::ofFQ, 3, (char*)formats,
	"format of output sequences/alignment, in any combination", NULL },
	{ 'o',"out",	0,	tNAME,	oOUTPUT, vUNDEF, 0, 0, NULL, OutFileTip.c_str()	},
#ifndef _NO_ZLIB
	{ 'z',"gzip",	0,	tENUM,	oOUTPUT, FALSE, vUNDEF, 2, NULL, "compress output files with gzip", NULL},
#endif
	{ 't', "time",	0,	tENUM,	oOTHER,	FALSE,	vUNDEF, 2, NULL, "print run time", NULL },
	{ 'V',"verbose",0,	tENUM,	oOTHER, vRT, vCRIT, vDEBUG+1, (char*)verbs,
	"\tset verbose level:\n? -\tshow critical messages only (silent mode)\n? -\tshow result summary\n?  -\tshow run-time information\n? -\tshow parameters\n? -\tshow debug messages", NULL },
	{ 'v', Version,	0,	tVERS,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print program's version", NULL },
	{ 'h', "help",	0,	tHELP,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print usage information", NULL }
};

const BYTE	Options::_OptCount = oHELP + 1;
const BYTE	Options::_UsageCount = 1;		// count of 'Usage' variants in help
const Options::Usage Options::_Usages[] = {	// content of 'Usage' variants in help
	{	vUNDEF,	" [template]\n\n  template - bed file whose features specify binding sites"	},
};

static const char* SignPar = "# ";	// Marker of output parameter in Usage

// Returns common name of output files
string GetOutFileName();
void PrintImitParams(const ChromFiles& cFiles, const char* templName, OutFile& oFile);
void PrintReadInfo();

/*****************************************/
int main(int argc, char* argv[])
{
	if (argc < 2)	return Options::PrintUsage(false);			// output tip
	int fileInd = Options::Tokenize(argc, argv);
	if( fileInd < 0 )	return 1;								// wrong otpion
	if(!Chrom::SetStatedID(Options::GetSVal(oCHROM))) return 1;	// wrong chrom name

	int ret = 0;
	BedF* templ = NULL;
	ChromSizes* cSizes = NULL;
	const char* fBedName = fileInd==argc ? NULL : argv[fileInd];	// template name

	Read::Init(
		Options::GetIVal(oREAD_LEN),
		Read::rNameType(Options::GetIVal(oREAD_NAME)),
		char(Options::GetIVal(oFQ_QUAL)),
		Options::GetIVal(oMAP_QUAL),
		Options::GetIVal(oREAD_LIMIT_N),
		// reduce limit because of thread's independent limit control
		(ULONG(Options::GetDVal(oREAD_LIMIT)) - (THREADS_CNT() >> 1))
	);
	Imitator::Init(
		RGL_SHIFT() > vUNDEF ? REGULAR : fBedName != NULL ? TEST : CONTROL,	// fBedName should be set
		ULONG(Options::GetDVal(oNUMB_CELLS)),
		Options::GetIVal(oVERB),
		Options::GetBVal(oBG_ALL) * bool(SAMPLE_BG()),	// bg all: always false if bg level is not 1.0
		Options::GetBVal(oLET_N),
		Options::GetBVal(oTS_UNIFORM),
		Options::GetBVal(oSTRAND_MIX),
		readlen(Options::GetIVal(oFLAT_LEN))
	);
	Imitator::InitFragLen(
		Options::GetIVal(oFRAG_LEN),
		Options::GetIVal(oFRAG_DEV),
		SZ_SEL()
	);
	LognormDistribution::Init(
		Options::GetFVal(oMEAN),
		Options::GetFVal(oSIGMA),
		Options::GetFVal(oLN_FACTOR),
		Options::GetFVal(oLN_TERM),
		Options::GetFVal(oSZ_SEL_SIGMA)
	);
	string outFileName = GetOutFileName();	// Imitator::Init() should be called before to define mode,
											// on which depends default name
	Amplification::Coefficient = short(Options::GetDVal(oAMPL));
	Random::SetSeed(!Options::GetBVal(oFIX));
	//setlocale(LC_ALL, strEmpty);

	// execution
	Mutex::Init();
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer::StartCPU();
	Timer timer;
	try {
		if(fBedName)	fBedName = FS::CheckedFileName(fBedName);
		ChromFiles cFiles(FS::CheckedFileDirName(oGFILE), Imitator::All);
		OutFile oFile(outFileName,
			OutFile::eFormat(Options::GetIVal(oFORMAT)),
			OutFile::eMode(Options::GetIVal(oSMODE)),
#ifdef _NO_ZLIB
			false
#else
			Options::GetBVal(oGZIP)
#endif
		);
		PrintImitParams(cFiles, fBedName, oFile);

		if(fBedName || oFile.IsSamSet())
			cSizes = new ChromSizes(cFiles);
		if(fBedName) {
			Obj::eInfo info = Imitator::Verbose(vDEBUG) ? Obj::iSTAT : Obj::iLAC;
			Obj::Ambig::SetSupplAlarm(Obj::Ambig::SHORT, "for given binding length");
			templ = new BedF(Template, fBedName, cSizes, info, false,
				readlen(Options::GetIVal(oBS_LEN)),					// binding length
				Imitator::Verbose(vDEBUG)							// print ambigs alarm
			);
			templ->Extend(1-Options::GetIVal(oBS_LEN), cSizes, info);
		}
		if( !cFiles.SetTreated(templ) )
			Err(Err::TF_EMPTY, fBedName, "features per selected chromosomes").Throw();
		oFile.Init(cSizes, Options::CommandLine(argc, argv));
		if(cSizes)	{ delete cSizes; cSizes = NULL; }
		
		Imitator(cFiles, oFile, templ).Execute();
	}
	catch(Err &e)				{ ret = 1; cerr << e.what() << endl; }
	catch(const exception &e)	{ ret = 1; cerr << e.what() << EOL; }
	catch(...)					{ ret = 1; cerr << "Unregistered error" << endl; }
	if(templ)	delete templ;
	if(cSizes)	delete cSizes;
	Timer::StopCPU(true);
	timer.Stop("wall-clock: ", false, true);
	Mutex::Finalize();
	//if( OPT_DEBUG() && Imitator::Verbose(vDEBUG))
	//	cout << "Locked memory: " << (mem0 - getAvailSystemMemory())/1024 << " Kb\n";
	return ret;
}

void PrintReadInfo()
{
	cout << SignPar;
	Read::Print();
	cout<< EOL << SignPar << "Include N along the edges" 
		<< SepCl << Options::GetBoolean(oLET_N) << endl;
}

const char* Equel = " = ";

void PrintImitParams(const ChromFiles& cFiles, const char* templName, OutFile& oFile)
{
	if( !Imitator::Verbose(vPAR) )	return;
	if( RegularMode )
		cout << SignPar << "REGULAR MODE\n";
	cout << SignPar << "Reference: genome" << SepCl << cFiles.Path() 
		 << SepCm << Chrom::Title << 's' << SepCl;
	if(cFiles.ChromsCount() == 1)	cout << Chrom::Name(CID(cFiles.cBegin()));
	else if(Chrom::StatedAll())		cout << "all";
	else							cout << Chrom::Name(Chrom::StatedID());
	cout << EOL;
	if(templName)
		cout << SignPar << "Template" << SepCl << templName << EOL;
	oFile.Print(SignPar);
	if( RegularMode ) {
		cout << SignPar << "Shift" << SepCl << RGL_SHIFT() << EOL << EOL;
		PrintReadInfo();
	}
	else {
		cout << SignPar << "Count of cells" << SepCl << ULONG(Options::GetDVal(oNUMB_CELLS)) << EOL;
		cout << SignPar << "Amplification" << SepCl;
		if(NoAmplification)	cout << Options::GetBoolean(false) << EOL;
		else				cout << Amplification::Coefficient << EOL;
		PrintReadInfo();
		if( TestMode ) {
			cout << SignPar << "Background for all chromosomes" << SepCl;
			if(bool(SAMPLE_BG()))	cout << Options::GetBoolean(oBG_ALL) << EOL;
			else					cout << "negligible due to zero background sample\n";
			cout << SignPar << "Input sample: foreground" << Equel << SAMPLE_FG()
				 << PERS << SepGroup << "background" << Equel << SAMPLE_BG() << PERS << EOL;
			cout << SignPar << "Binding length" << SepCl << Options::GetIVal(oBS_LEN) << EOL;
			cout << SignPar << "Boundary flattening length" << SepCl << Options::GetIVal(oFLAT_LEN) << EOL;
			cout << SignPar << "Strand admixture" << SepCl << Options::GetBoolean(oSTRAND_MIX) << EOL;
			cout << SignPar << "Uniform template score"  << SepCl
				 << (Options::GetBVal(oTS_UNIFORM) ? "YES" : "NO") << EOL;
		}
		cout << SignPar << "Lognormal distribution" << SepCl
			 << "sigma" << Equel << Options::GetIVal(oSIGMA)
			 << SepGroup << "mean" << Equel << Options::GetIVal(oMEAN) 
			 << SepGroup << "lnFactor" << Equel << Options::GetIVal(oLN_FACTOR)
			 << SepGroup << "lnTerm" << Equel << Options::GetDVal(oLN_TERM) << EOL;
		cout << SignPar << "Frag's size selection" << SepCl;
		if(SZ_SEL())
			cout << "length" << Equel << Options::GetIVal(oFRAG_LEN) 
				 << SepGroup << "deviation" << Equel << Options::GetIVal(oFRAG_DEV)
				 << SepGroup << "SD of the normal distribution" << Equel
				 << Options::GetFVal(oSZ_SEL_SIGMA) << EOL;
		else
			cout << Options::GetBoolean(false) << EOL;
		if(THREADS_CNT() > 1)
			cout << SignPar << int(THREADS_CNT()) << " threads\n";
	}
	cout << endl;
}

// Returns common name of output files
string GetOutFileName()
{
	const char* outName_c = Options::GetSVal(oOUT_FILE);
	
	if( !outName_c )				return DefFileName[Imitator::Mode];
	if( FS::IsDirExist(outName_c) )	return FS::MakePath(string(outName_c)) + DefFileName[Imitator::Mode];

	string outName = string(outName_c);
	if( FS::HasExt(outName) ) {
		cout << "discarded extention in " << outName << EOL;
		outName = FS::FileNameWithoutExt(outName);
	}
	return outName;
}

//#ifdef DEBUG
//void	CheckFaReadWrite()
//{
//	Mutex::Init();
//	CPU_Timer::Enabled = true;
//	CPU_Timer::Start();
//	try {
//		//TxtFile file("\\documents\\prof\\icl\\data\\aaa.txt", TxtFile::READ, 1);
//		//char buff[100];
//		////memset(buff, 0, 100);
//		//char *read;
//		//int i=0;
//		//while( read=file.GetRecord(NULL) ) {
//		//	memset(buff, 0, 100);
//		//	memcpy(buff,read,file.RecordLength()-file._EOLSize);
//		//	cout << (int)file._EOLSize << '\t' << (int)(file.RecordLength()-file._EOLSize) << '\t' << buff << endl;
//		//	i++;
//		//}
//		//cout << "===========\n";
//		//cout << "count " << i << endl;
//		//cout << file.RecordCount() << endl;
//		string fname = "\\documents\\prof\\icl\\data\\mousegenome9\\chr";
//		Chrom chr(fname + "1.fa.gz", false);
//		//Chrom chr(fname + "1.fa", false);
//		cout << "full =\t" << chr.Length() << endl;
//		cout << "short =\t" << (chr.End() - chr.Start()) << endl;
//		cout << "n =\t" << chr.CountN() << endl;
//		//chr.Write(fname + "1_.fa", "1");
//
//		//Sleep(1);
//		//Chrom chr1(fname + "1_.fa", false);
//		//if( chr.Length() != chr1.Length() )		cout << "Error length!\n";
//		//else {
//		//	for(ULONG i=0; i<chr.Length(); i++)
//		//		if( chr._nts[i] != chr1._nts[i] ) {
//		//			cout << "Error!\n";
//		//			break;
//		//		}
//		//	cout << "Done.\n";
//		//}
//	}
//	catch(const Msg&e) { cerr << e.what() << endl; }
//	catch(...) { cerr << "error\n"; }
//	CPU_Timer::Stop();
//	Mutex::Finalize();
//}
//#endif

	//Timer tm(true);
	//int i, lim = 20000000;
	//int x = 12345;
	//string str;
	//char buf[10];
	//tm.Start();
	//for(i=0; i<lim; i++)
	//	_itoa(x, buf, 10);
	//cout << buf << endl;
	//tm.Stop();	// out: 12345 Windows: 0:01
	//
	//tm.Start();
	//x++;
	//for(i=0; i<lim; i++)
	//	_gcvt(x, DigitsCount(x), buf);
	//cout << buf << endl;
	//tm.Stop();	// out: 12346 Windows: 0:18  Linux: 00:11

	//tm.Start();
	//x++;
	//for(i=0; i<lim; i++)
	//	str = NSTR(x);
	//cout << str << endl;
	//tm.Stop();	// out: 12347 Windows: 0:44  Linux: 00:11
	//return 1;
