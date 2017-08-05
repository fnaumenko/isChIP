 /*
  *  Copyright (C) 2017  Fedor Naumenko
  *
  *   This program is free software: you can redistribute it and/or modify
  *   it under the terms of the GNU General Public License as published by
  *   the Free Software Foundation, either version 3 of the License, or
  *   (at your option) any later version.
  *
  *   This program is distributed in the hope that it will be useful,
  *   but WITHOUT ANY WARRANTY; without even the implied warranty of
  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *   GNU General Public License for more details.
  *   Path to UCSC reference genome: ftp://hgdownload.soe.ucsc.edu/goldenPath/
  *
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

const char* smodes	[] = { "SE", "PE" };			// values of --smode option
const char* rnames	[] = { "NMB", "POS" };			// values of --read-name option
const char* formats	[] = { "FQ", "BED", "SAM" };	// values of --format option
const char* verbs	[] = { "CRIT", "RES", "RT", "PAR", "DBG" };	// values of --verbose option

//{ char, str, optOblig, valRequired, type, group, defVal, minVal, maxVal, strVal, descr }
Options::Option Options::_Options [] = {
	{ 'a',"amplif",	0, true, tINT,	oTREAT, 1, 1, 2000, NULL, "coefficient of amplification" },
	{ 'b',"bg-level",0,true, tFLOAT,oTREAT, 1, 0, 100, NULL,
	"background of out-of-features fragments, in percent of foreground.\nFor test mode only" },
	{ HPH,"fg-level",0,true, tFLOAT,oTREAT, 100, 0, 100, NULL,
	"foreground of in-features fragments, in percent.\nIn control mode fragments sample" },
	{ 'n',"cells",	0, true, tLONG,	oTREAT, 1, 1, 1e7, NULL, "number of cells" },
	{ 'g',"gen",	1, true, tNAME, oINPUT, vUNDEF, 0, 0, NULL,
	"reference genome (file name or directory). Required" },
	{ 'c',Chrom::Abbr,0,true,tCHAR, oTREAT, vUNDEF, 0, 0, NULL,
	"generate output for stated chromosome only" },
	{ HPH,"flat-len",0,true, tINT,	oTREAT, 0, 0, 200, NULL,
	"boundary flattening length. For test mode only" },
	{ HPH,"frag-len",0,true, tINT,	oFRAG, 200, 50, 400, NULL, "average size of selected fragments"	},
	{ HPH,"frag-dev",0,true, tINT,	oFRAG, 20, 0, 200, NULL, "deviation of selected fragments" },
	{ HPH,"bg-all",	0, true, tENUM,	oTREAT, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off generation background for all chromosomes.\nFor test mode only" },
	{ HPH,"bind-len",0, true, tINT,	oTREAT, 1, 1, 100, NULL,
	"minimum binding length. For test mode only" },
	{ HPH,"mean",	0, true, tINT,	oDISTR, 200, 0, 1500, NULL,
	"expectation of the based normal distribution" },
	{ HPH,"sigma",	0, true, tINT,	oDISTR, 200, 1, 700, NULL,
	"standard deviation of the based normal distribution" },
	{ HPH,"ln-factor",	0, true, tINT,	oDISTR, 500, 10, 900, NULL,
	"power multiplication factor in lognormal distribution" },
	{ HPH,"ln-term",	0, true, tFLOAT,oDISTR, 5.1, 1, 7, NULL,
	"power summand in lognormal distribution" },
	{ HPH,"let-N",		0, false,tENUM,	oTREAT, FALSE, 0, 2, NULL,
	"include the ambiguous reference characters (N) on the beginning\nand on the end of chromosome" },
	{ 'r',"read-len",	0, true, tINT,	oREAD, 50, 20, 200, NULL, "length of output read" },
	{ HPH,"read-name",	0, true, tENUM, oREAD, Read::nmPos, Read::nmNumb, Read::nmUndef, (char*)rnames,
	"name of read in output files includes:\n? - read`s unique number within chromosome\n? - read`s true start position" },
	{ HPH,"read-Nlimit",0, true, tINT,	oREAD, vUNDEF, 0, 100, NULL,
	"maximum permitted number of ambiguous characters (N) in read [--read-len]" },
	{ HPH,"reads-limit",0, true, tLONG,	oREAD, 2e8, 1e5, (float)ULONG_MAX, NULL,
	"maximum permitted number of total written reads" },
	{ HPH,"fq-qual",	0, true, tCHAR,	oREAD, '~', '!', '~', NULL,
	"the quality values for the read in FQ output" },
	{ HPH,"map-qual",	0, true, tINT,	oREAD, 42, 0, 42, NULL,
	"the mapping quality for the read in SAM output" },
	{ HPH,"sz-sel",		0, true, tENUM,	oFRAG, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off fragment's size selection" },
	{ HPH,"sz-sel-sigma",0,true,tINT,oFRAG, 20, 1, 100, NULL,
	"standard deviation of the size selection normal distribution" },
//{ HPH, "stdout",	0, false,tBOOL,	FALSE, 0, 0, NULL, "Write output to the 'standard output' filehandle instead of 'standart error'" },
	{ HPH,"smode",	0, true,tENUM,	oTREAT, 0, 0, 2, (char*)smodes,
	"sequencing mode: ? - single-end, ? - paired-end" },
	{ HPH,"strand-admix",0, true,tENUM,oTREAT, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off opposite strand admixture at the bound\nof binding site. For test mode only" },
	{ HPH,"ts-uni",	0, false,tENUM, oTREAT, FALSE, 0, 2, NULL, "uniform template score" },
	{ 'p',"threads",0, true, tINT,	oTREAT, 1, 1, 50, NULL, "number of threads"},
	{ HPH,"debug",	0, false,tENUM, oTREAT, FALSE, 0, 2, NULL,
	"fix random emission to get repetitive results" },
	{ 'R',"regular",	0, true, tINT,	oTREAT, vUNDEF, 1, 400, NULL,
	"regular mode: write each read on starting position\nincreased by stated shift" },
	{ 'f',"format",	0, true, tCOMB,	oOUTPUT, 1, 1, 3, (char*)formats,
	"format of output sequences/alignment, in any combination"},
	{ 'o',"out",	0, true, tNAME,	oOUTPUT, vUNDEF, 0, 0, NULL, OutFileTip.c_str()	},
	{ 'z',"gzip",	0, false,tENUM,	oOUTPUT, FALSE, 0, 2, NULL, "compress output files with gzip"},
	{ 't',"time",	0, false,tENUM, oOTHER, FALSE, 0, 2, NULL, "print run time" },
	{ 'V',"verbose",0, true, tENUM, oOTHER, V_RT, V_CRIT, V_DEBUG+1, (char*)verbs,
	"\tset verbose level:\n? -\tshow critical messages only (silent mode)\n? -\tshow result summary\n?  -\tshow run time information\n? -\tshow process parameters\n? -\tshow debug messages" },
	{ 'v',"version",0, false,tVERS, oOTHER, vUNDEF, 0, 0, NULL,
	"print program and ZLib version, and quit" },
	{ 'h',"help",	0, false,tHELP, oOTHER, vUNDEF, 0, 0, NULL,
	"print usage information and quit" }
};

const BYTE	Options::_OptCount = oHELP + 1;
const BYTE	Options::_UsageCount = 1;		// count of 'Usage' variants in help
const Options::Usage Options::_Usages[] = {	// content of 'Usage' variants in help
	{	vUNDEF,	" [template]"	},
};

// Returns default output file name without extention
inline string GetDefOutFileName()	{ return DefFileName[Imitator::Mode]; }

void PrintImitParams(const char* genFileName, OutFile& oFile);

/*****************************************/
int main(int argc, char* argv[])
{
	if (argc < 2)	return Options::PrintUsage(false);			// output tip
	short fileInd = Options::Tokenize(argc, argv);
	if( fileInd < 0 )	return 1;								// wrong otpion
	
	chrid cID = Chrom::ID(Options::GetSVal(oCHROM));
	if(!Chrom::CheckID(cID)) {
		cerr << Options::GetSVal(oCHROM) << ": wrong " << Chrom::Title << " name\n";
		return 1;
	}

	int ret = 0;
	const char* fBedName = fileInd==argc ? NULL : argv[fileInd];// template name
	BedF* templ = NULL;

	Read::Init(
		Options::GetIVal(oREAD_LEN),
		static_cast<Read::eName>(Options::GetIVal(oREAD_NAME)),
		char(Options::GetIVal(oFQ_QUAL)),
		Options::GetIVal(oMAP_QUAL),
		Options::GetIVal(oREAD_LIMIT_N),
		// reduce limit because of thread's independent limit control
		(ULONG(Options::GetDVal(oREAD_LIMIT)) - (THREADS_CNT() >> 1))
	);

	// set imitator parameters; used fBedName
	Imitator::Init(
		RGL_SHIFT() > vUNDEF ? 
			REGULAR :
			fBedName != NULL ? TEST : CONTROL,
		ULONG(Options::GetDVal(oNUMB_CELLS)),
		Options::GetIVal(oVERB),
		// bg all: always false if bg level is not 1.0
		Options::GetBVal(oBG_ALL) * bool(SAMPLE_BG()),
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
	Amplification::Coefficient = short(Options::GetDVal(oAMPL));

	// set out file name
	// Imitator::Init() should be called before to define default name
	string foutName;
	{	// set fq file name
		const char* foutName_c = Options::GetSVal(oOUT_FILE);
		if( foutName_c ) {
			foutName = string(foutName_c);
			if( FS::IsDirExist(foutName_c) )
				foutName += SLASH + GetDefOutFileName();
			else if( FS::HasExt(foutName) ) {
				cout << "discarded extention in " << foutName << EOL;
				foutName = FS::FileNameWithoutExt(foutName);
			}
		}
		else
			foutName = GetDefOutFileName();
	}
	Random::SetSeed(!Options::GetBVal(oDBG));
	//setlocale(LC_ALL, StrEmpty);

	// execution
	Timer timer(Options::GetBVal(oTIME));
	Timer::StartCPU();
	timer.Start();
	Mutex::Init();
	try {
		//if(fBedName)	SamFragDistr(string(fBedName));
	
		const char* genName	= FS::CheckedFileDirName(oGFILE);
		if(fBedName) {
			bool verbParam = Imitator::Verbose(V_PAR);
			Bed::Ambig::SetTotalAlarm(Bed::Ambig::SHORT, "for given binding length");
			if(verbParam)					cout << SignPar;
			if(Imitator::Verbose(V_RT))		cout << Template << MSGSEP_BLANK;
			templ = new BedF(
				FS::CheckedFileName(fBedName), cID,
				verbParam,							// print file name
				readlen(Options::GetIVal(oBS_LEN)),	// binding length
				Imitator::Verbose(V_DEBUG),			// print ambigs alarm
				verbParam,							// print statistics
				Imitator::Verbose(V_RT)				// print number of readed features
			);
			templ->Expand(1-Options::GetIVal(oBS_LEN));
			//templ->Print(3);
		}
		ChromFiles cFiles(genName, cID, Imitator::All);
		if( !cFiles.SetTreated(templ) )
			Err(Err::TF_EMPTY, fBedName, "features in selected chromosomes").Throw();
		OutFile oFile(foutName,
			Options::GetIVal(oFORMAT),
			Options::GetIVal(oSMODE),
			Options::GetBVal(oGZIP)
		);
		PrintImitParams(genName, oFile);

		oFile.Init(cFiles, Options::CommandLine(argc, argv));
		Imitator(cFiles, oFile, templ).Execute();
	}
	catch(Err &e)				{ ret = 1; cerr << e.what() << endl; }
	catch(const exception &e)	{ ret = 1; dout << e.what() << EOL; }
	catch(...)					{ ret = 1; cerr << "Unregistered error" << endl; }
	if( templ )	delete templ;
	Timer::StopCPU(true);
	timer.Stop("wall-clock: ", false, true);
	Mutex::Finalize();
	//if( OPT_DEBUG() && Imitator::Verbose(V_DEBUG))
	//	cout << "Locked memory: " << (mem0 - getAvailSystemMemory())/1024 << " Kb\n";
	return ret;
}

void PrintReadInfo()
{
	cout << SignPar;
	Read::Print();
	cout << EOL << SignPar << "Include N along the edges: " 
			<< Options::GetBoolean(oLET_N) << endl;
}

void PrintImitParams(const char* genFileName, OutFile& oFile)
{
	if( !Imitator::Verbose(V_PAR) )	return;
	if( RegularMode )
		cout << SignPar << "REGULAR MODE\n";
	if( !FS::IsDirExist(genFileName) )
		cout << SignPar << "Input Fa file: " << genFileName << EOL;
	oFile.Print(SignPar);
	if( RegularMode ) {
		cout << SignPar << "Shift: " << RGL_SHIFT() << EOL << EOL;
		PrintReadInfo();
	}
	else {
		cout << SignPar << "Count of cells: " << ULONG(Options::GetDVal(oNUMB_CELLS)) << EOL;
		cout << SignPar << "Amplification: ";
		if(NoAmplification)	cout << Options::GetBoolean(false) << EOL;
		else				cout << Amplification::Coefficient << EOL;
		PrintReadInfo();
		if( TestMode ) {
			cout << SignPar;
			if( Options::GetIVal(oFLAT_LEN) )
				cout << "Boundary flattening length: " << Options::GetIVal(oFLAT_LEN);
			else
				cout << "Binding length: " << Options::GetIVal(oBS_LEN);
			cout << EOL << SignPar << "Input sample: foreground = " << SAMPLE_FG()
				 << PERS << GroupParSep << "background = " << SAMPLE_BG() << PERS << EOL;
			cout << SignPar
				 << "Background for all chromosomes: ";
			if(bool(SAMPLE_BG()))
				cout << Options::GetBoolean(oBG_ALL) << EOL;
			else
				cout << "negligible due to zero background sample\n";
		}
		cout << SignPar << "Lognormal distribution: "
			 << "sigma = " << Options::GetIVal(oSIGMA)
			 << GroupParSep << "mean = " << Options::GetIVal(oMEAN) 
			 << GroupParSep << "lnFactor = " << Options::GetIVal(oLN_FACTOR)
			 << GroupParSep << "lnTerm = " << Options::GetDVal(oLN_TERM) << EOL;
		cout << SignPar << "Frag's size selection: ";
		if(SZ_SEL())
			cout << "length = " << Options::GetIVal(oFRAG_LEN) 
				 << GroupParSep << "deviation = " << Options::GetIVal(oFRAG_DEV) << EOL;
		else
			cout << Options::GetBoolean(false) << EOL;
		if(THREADS_CNT() > 1)
			cout << SignPar << int(THREADS_CNT()) << " threads\n";
	}
	cout << endl;
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
