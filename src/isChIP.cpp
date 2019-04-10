/************************************************************************************
isChIP (In-Silico ChIP) is a fast realistic ChIP-seq simulator.
The model is based on the real protocol of ChIP-seq 

Copyright (C) 2014-2019 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 10.04.2019
-------------------------

This program is free software.
It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
GNU General Public License for more details.

Path to UCSC reference genome: ftp://hgdownload.soe.ucsc.edu/goldenPath/
************************************************************************************/

#include "isChIP.h"
#include "Imitator.h"

using namespace std;

const string Product::Title = "isChIP";
const string Product::Version = "1.0";
const string Product::Descr = "ChIP-seq simulator";

const string DefFileName [] = { "mTest", "mInput" };
const string OutFileTip = "location of output files or existing directory\n[TEST mode: " + 
	DefFileName[TEST] + ".*, " + 
	"CONTROL mode: " + DefFileName[CONTROL] + ".*]";

// options groups
enum optGroup { gTREAT, gFRAG, gREAD, gOUTPUT, gOTHER };	// gOTHER should be the last 

const char* Options::OptGroups [] = {
	"Processing", "Fragment distribution", "Reads", "Output", "Other"
};
const BYTE Options::GroupCount = sizeof(Options::OptGroups)/sizeof(char*);

// --ampl option
const char* ampls	[] = { "MDA","PCR" };				// corresponds to OutFile::eMode
// --smode option
const char* smodes	[] = { "SE","PE" };					// corresponds to OutFile::eMode
// -read-name option: types of Read name notations
const char* rnames	[] = { "NONE","NUMB","POS" };		// corresponds to Read::rNameType
// --format option: format notations
const char* formats	[] = { "FQ","BED","SAM","WIG","FREQ" };	// corresponds to OutFile::oFormat	
// --verbose option: verbose notations
const char* verbs	[] = { "SL","RES","RT","PAR","DBG" };
// --flat-len option
pairVal flattens = make_pair(0, 0);

const char* ForTest = "For the test mode only";

// { char, str, Signs, type, group, defVal, minVal, maxVal, strVal, descr, addDescr }
// defVal: vUNDEF if no default value should be printed
// minVal: vUNDEF if value is prohibited
Options::Option Options::List [] = {
	{ 'g',"gen",	1,	tNAME,	gTREAT, vUNDEF, 0, 0, NULL, "reference genome library.", NULL },
	{ 'n',"cells",	0,	tINT,	gTREAT, 1, 1, 2e6, NULL, "number of cells", NULL },
	{ HPH,"ampl",	2,	tENUM,	gTREAT, 1, 0, 2, (char*)ampls, "amplification technique", NULL },
	{ 'a',"ampl-c",	0,	tINT,	gTREAT, 0, 0, 100, NULL, "number of PCR cycles", NULL },
	//"number of displacement reactions (MDA) | cycles (PCR)", NULL },
	{ 'b', "bg",	0,	tFLOAT,	gTREAT, 1, 0, 100, NULL,
	"number of selected fragments outside BSs,\nin percent of foreground.", ForTest },
	{ HPH, "fg",	0,	tFLOAT,	gTREAT, 100, 0, 100, NULL,
	"in test mode the number of selected fragments within BSs, in percent;\nin control mode the number of selected fragments, in percent", NULL },
	{ HPH,"bg-all",	0,	tENUM,	gTREAT, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off generation background for all chromosomes.\n", ForTest },
	{ 'c',Chrom::Abbr,0,tNAME,	gTREAT, vUNDEF, 0, 0, NULL,
	"generate output for the specified chromosome only", NULL },
	{ HPH,"fr-mean",0,	tFLOAT,	gFRAG, 5.46, 2, 9, NULL,
	"mean of fragment lognormal distribution", NULL },
	{ HPH,"fr-sd",	0,	tFLOAT,	gFRAG, 0.4, 0.1, 1., NULL,
	"standard deviation of fragment lognormal distribution", NULL },
	{ HPH,"ss-mean",0,	tINT,	gFRAG, vUNDEF, 50, 1000, NULL,
	"mean of size selection normal distribution [auto]", NULL },
	{ HPH, "ss-sd",	0,	tINT,	gFRAG, 30, 1, 200, NULL,
	"standard deviation of size selection normal distribution", NULL },
	{ HPH, "ss",	0,	tENUM,	gFRAG, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off fragment's size selection", NULL },
	{ 'N', "full-size",	0,	tENUM,	gTREAT, FALSE, vUNDEF, 2, NULL,
	"process the entire reference chromosome (including gaps at each end)", NULL },
	//{ HPH, "bind-len",	0,	tINT,	gTREAT, 1, 1, 100, NULL, "minimum binding length.", ForTest },
	{ HPH, "flat-len",	2,	tPAIR,	gTREAT, 0, 0, 200, (char*)&flattens,
	"inside, outside BS boundary flattening length.", ForTest },
	{ 'r', "rd-len",	0,	tINT,	gREAD, 50, 20, 500, NULL, "length of output read", NULL },
	{ HPH, "rd-name",	0,	tENUM,	gREAD, Read::nmNone, Read::nmNone, Read::nmPos, (char*)rnames,
	"info added to read's name in output files:\n? - nothing\n? - read`s unique number across genome\n?  - read`s actual start position", NULL },
	{ HPH, "rd-Nlim",	0,	tINT,	gREAD, vUNDEF, 0, 500, NULL,
	"maximum permitted number of ambiguous code N in read [OFF]", NULL },
	{ HPH, "rd-lim",	0,	tLONG,	gREAD, 2e8, 1e5, (float)ULONG_MAX, NULL,
	"maximum permitted number of total recorded reads", NULL },
	{ HPH, "rd-ql",		0,	tCHAR,	gREAD, '~', '!', '~', NULL, "uniform read quality score", NULL },
	{ HPH,"rd-ql-patt",	0,	tNAME,	gREAD, vUNDEF, 0, 0, NULL, "read quality scores pattern", NULL },
	{ HPH,"rd-ql-map",	0,	tINT,	gREAD, 255, 0, 255, NULL,
	"read mapping quality in SAM and BED output", NULL },
//{ HPH, "stdout",	0, tBOOL,	FALSE, 0, 0, NULL, "Write output to the 'standard output' filehandle instead of 'standart error'" },
	{ 'm', "smode",		0,	tENUM,	gTREAT, Seq::SE, Seq::SE, Seq::Undef, (char*)smodes,
	"sequencing mode: ? - single end, ? - paired end", NULL },
	{ HPH,"strand-err",	2,	tFLOAT, gTREAT, 0, 0, 100, NULL,
	"percentage of reads with wrong strand", NULL },
	{ 'u',"ts-uni",	0,	tENUM,	gTREAT, FALSE, vUNDEF, 2, NULL, "uniform template score.", ForTest },
	{ 'p',"threads",0,	tINT,	gTREAT, 1, 1, 20, NULL, "number of threads", NULL },
	{ HPH, "seed",	0,	tINT,	gTREAT, 0, 0, 100, NULL,
	"fix random emission with given seed, or 0 if don't fix", NULL },
	{ 'f',"format",	0,	tCOMB,	gOUTPUT, OutFiles::ofFQ, OutFiles::ofFQ, 5, (char*)formats,
	"format of output data, in any order", NULL },
	{ 'C',"control",0,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL, 
	"generate control simultaneously with test", NULL },
	{ 'S', "strand",0,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL,
	"generate two additional wig files, each one per strand", NULL },
	{ 's', "sep",	0,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL,
	"display reads number thousands separator", NULL },
	{ 'o', "out",	0,	tNAME,	gOUTPUT, vUNDEF, 0, 0, NULL, OutFileTip.c_str()	},
#ifndef _NO_ZLIB
	{ 'z',"gzip",	0,	tENUM,	gOUTPUT, FALSE, vUNDEF, 2, NULL, "compress the output", NULL},
#endif
	{ 't',"time",	0,	tENUM,	gOTHER,	FALSE,	vUNDEF, 2, NULL, "print run time", NULL },
	{ 'V',"verbose",0,	tENUM,	gOTHER, vRT, vCRIT, vDBG+1, (char*)verbs,
	"set verbose level:\n? -\tsilent mode (show critical messages only)\n? -\tshow result summary\n?  -\tshow run-time information\n? -\tshow actual parameters\n? -\tshow debug messages", NULL },
	{ 'v',Version,	0,	tVERS,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, "print program's version", NULL },
	{ 'h',"help",	0,	tHELP,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, "print usage information", NULL }
};
const BYTE	Options::OptCount = sizeof(Options::List)/sizeof(Options::Option);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, "<template>", false, "bed file whose features specify binding sites (BS)"	},
};
const BYTE Options::UsageCount = sizeof(Options::Usages)/sizeof(Options::Usage);

// Returns common name of output files
string GetOutFileName();
void PrintImitParams(const ChromSizesExt& cSizes, const char* templName, const BedF* templ, const OutFiles& oFile);

/*****************************************/
int main(int argc, char* argv[])
{
	if (argc < 2)	return Options::PrintUsage(true);			// output tip
	int fileInd = Options::Parse(argc, argv);
	if( fileInd < 0 )	return 1;								// wrong otpion

	//MemStatus::StartObserve(Options::GetUIVal(oVERB)==vDBG);

	int ret = 0;
	BedF* templ = NULL;
	const char* fBedName = fileInd==argc ? NULL : argv[fileInd];	// template name

	// initialize Seq before Read::Init()
	Seq::Init(Seq::sMode(Options::GetIVal(oSMODE)), ULONG(Options::GetFVal(oRD_LIMIT)));
	// initialize Read before OutFiles::Init()
	Read::Init(
		Options::GetUIVal(oRD_LEN),
		Read::rNameType(Options::GetIVal(oRD_NAME)),
		Seq::IsPE() && !Options::Assigned(oRD_NAME),
		char(Options::GetFVal(oRD_QUAL)),
		Options::GetUIVal(oRD_LIMIT_N)
	);
	OutFiles::Init(
		OutFiles::oFormat(Options::GetIVal(oFORMAT)),
		Options::GetUIVal(oMAP_QUAL),
		Options::GetBVal(oSTRAND),
		Options::GetFVal(oSTRAND_ERR) / 100,
#ifdef _NO_ZLIB
		false
#else
		Options::GetBVal(oGZIP)
#endif
	);
	// initialize DistrParams before Imitator::Init()
	DistrParams::Init(
		Options::GetFVal(oFRAG_MEAN),
		Options::GetFVal(oFRAG_SIGMA),
		Options::GetFVal(oSS_MEAN),
		Options::GetFVal(oSS_SIGMA),
		Options::GetBVal(oSS)
	);
	Imitator::Init(
		fBedName != NULL ? TEST : CONTROL,	// fBedName should be set
		Options::GetBVal(oMAKE_INPUT),
		ULONG(Options::GetFVal(oNUMB_CELLS)),
		Options::GetBVal(oAMPL),
		a_coeff(Options::GetFVal(oAMPL_COEF)),
		Options::GetUIVal(oVERB),
		Options::GetBVal(oBG_ALL) * bool(SAMPLE_BG()),	// bg all: always false if bg level is 0
		Options::GetBVal(oTS_UNIFORM),
		//readlen(Options::GetUIVal(oBS_LEN)),
		flattens
	);
	RefSeq::LetGaps = Options::GetBVal(oLET_GAPS);
	RefSeq::StatGaps = Imitator::Verbose(vPAR);	// Imitator::Init() should be called before
	Random::SetSeed(Options::GetUIVal(oSEED));
	if(Options::GetBVal(oLOCALE))	cout.imbue(locale(LOCALE_ENG));

	// execution
	Mutex::Init(Options::GetIVal(oNUMB_THREAD)>1);
	Timer::Enabled = Options::GetBVal(oTIME) && Imitator::Verbose(vRES);
	Timer::StartCPU();
	Timer timer;
	try {
		// check file names first of all
		FS::CheckedFileName(fBedName);
		const char* qlPatt = FS::CheckedFileName(oRD_QUAL_PATT);	// read quality pattern

		ChromSizesExt cSizes(FS::CheckedFileDirName(oGFILE), Imitator::Verbose(vRT));
		if(!Chrom::SetCustomID(Options::GetSVal(oCHROM)))		// should be invoke after ChromSizes()
			Err(Err::OPT_CHROM, Options::OptionToStr(oCHROM).c_str()).Throw();

		if(fBedName) {
			Obj::eInfo info = Obj::iNONE;
			if(Imitator::Verbose(vPAR))		info = Obj::iLAC;
			if(Imitator::Verbose(vDBG))		info = Obj::iSTAT;
			//Obj::Ambig::SetSupplAlarm(Obj::Ambig::SHORT, "for given binding length");
			templ = new BedF(Template, fBedName, &cSizes, info, false,
				0, //readlen(Options::GetUIVal(oBS_LEN)),		// binding length
				Imitator::Verbose(vDBG)							// print ambigs alarm
			);
			if(flattens.second)
				templ->Extend(flattens.second, &cSizes, Obj::iLAC);
		}
		// set treated chroms
		if( !cSizes.SetTreated(Imitator::All, templ) )
			Err(Err::TF_EMPTY, fBedName, "features per stated " + Chrom::ShortName(Chrom::CustomID())).Throw();
		
		Imitator::SetThreadNumb( min(chrid(Options::GetFVal(oNUMB_THREAD)), cSizes.TreatedCount()) );

		OutFiles oFile(GetOutFileName(), Imitator::IsControl(),
			cSizes, qlPatt, Options::CommandLine(argc, argv));

		PrintImitParams(cSizes, fBedName, templ, oFile);
		Imitator(cSizes, oFile).Execute(templ);
	}
	catch(Err &e)				{ ret = 1; cerr << e.what() << endl; }
	catch(const exception &e)	{ ret = 1; cerr << e.what() << EOL; }
	catch(...)					{ ret = 1; cerr << "Unregistered error" << endl; }
	if(templ)	delete templ;
	Timer::StopCPU(false);
	timer.Stop("\twall-clock: ", false, true);
	Mutex::Finalize();
	//MemStatus::StopObserve();
	return ret;
}

void PrintImitParams(const ChromSizesExt& cSizes, const char* templName, 
	const BedF* templ, const OutFiles& oFile)
{
	if( !Imitator::Verbose(vPAR) )	return;
	cout << SignPar << "Reference" << SepDCl << "genome" << SepCl << cSizes.Path();
	cout << SepCm << Chrom::TitleName() << SepCl;
	if(!Chrom::NoCustom())		
		cout << Chrom::Mark(Chrom::CustomID());
	else if(cSizes.TreatedCount() < cSizes.ChromCount())	
		cSizes.PrintTreatedChroms();
	else	
		cout << "all";
 	cout << EOL;

	if(templName) {
		cout << SignPar << Template << SepCl << templName << SepCl;
		templ->PrintItemCount();
	}
	oFile.PrintFormat(SignPar);		// print output formats, sequencing mode
	cout << SignPar; Seq::Print();
	cout << SignPar << "Count of cells" << SepCl << ULONG(Options::GetFVal(oNUMB_CELLS)) << EOL;
	Imitator::PrintAmpl(ampls);
	cout << SignPar;
	Read::Print();
	oFile.PrintReadQual(SignPar);
	cout << SignPar << "Optimization: process the entire ref. chromosome"
		 << SepCl << Options::BoolToStr(oLET_GAPS) << EOL;
	cout << SignPar << "Stated sample: ";
	if(ControlMode)		cout << SAMPLE_FG() << PERS << EOL;
	else {
		cout << "foreground" << Equel << SAMPLE_FG()
				<< PERS << SepGroup << "background" << Equel << SAMPLE_BG() << PERS << EOL;
		if(cSizes.TreatedCount() > 1) {
			cout << SignPar << "Background for all chromosomes" << SepCl;
			if(bool(SAMPLE_BG()))	cout << Options::BoolToStr(oBG_ALL) << EOL;
			else					cout << "negligible due to zero background sample\n";
		}
		//cout << SignPar << "Binding length" << SepCl << Options::GetUIVal(oBS_LEN) << EOL;
		if(flattens.first || flattens.second)
			cout << SignPar << "BS boundary flattening length: inside" << Equel << flattens.first
				 << SepGroup << "outside" << Equel << flattens.second << EOL;
		cout << SignPar << "Uniform template score"  << SepCl
				<< (Options::GetBVal(oTS_UNIFORM) ? "YES" : "NO") << EOL;
	}
	float ws = Options::GetFVal(oSTRAND_ERR);
	if(ws)	cout << SignPar << "Reads with wrong strand"
					<< SepCl << ws << PERS << EOL;
	if(!Imitator::IsSingleThread())
		cout << SignPar << "Actual threads" << SepCl << int(Imitator::ThrCnt) << EOL;
	cout << SignPar << "Fragment lognormal distribution"
			<< SepCl << "mean" << Equel << Options::GetFVal(oFRAG_MEAN) 
			<< SepGroup << "stand.dev" << Equel << Options::GetFVal(oFRAG_SIGMA) << EOL;
	cout << SignPar << "Fragment size selection" << SepCl;
	if(DistrParams::IsSS())
		cout << "mean" << Equel << setprecision(5) << DistrParams::ssMean 
				<< SepGroup << "stand.dev" << Equel << DistrParams::ssSigma << EOL;
	else
		cout << Options::BoolToStr(false) << EOL;
}

// Returns common name of output files
string GetOutFileName()
{
	const char* outName_c = Options::GetSVal(oOUT_FILE);
	
	if( !outName_c )				return DefFileName[Imitator::TMode];
	if( FS::IsDirExist(outName_c) )	return FS::MakePath(string(outName_c)) + DefFileName[Imitator::TMode];
	return string(outName_c);
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
//		//		if( chr._seq[i] != chr1._seq[i] ) {
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


	//cout.imbue(locale(""));
	//cout.precision(0);
	//cout << fixed;

	//setlocale(LC_ALL, "");
	//cout << "LC_ALL: " << setlocale(LC_ALL, NULL) << endl;
	//cout << "LC_NUMERIC: " << setlocale(LC_NUMERIC, NULL) << endl;

	//locale loc;
	////std::cout << "The global locale is: " << loc.name() << EOL;
	////loc = getloc();
	//cout << "The global locale is: " << loc.name() << EOL;
