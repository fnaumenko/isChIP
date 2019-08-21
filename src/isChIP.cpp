/************************************************************************************
isChIP (In-Silico ChIP) is a fast realistic ChIP-seq simulator.
The model is based on the real protocol of ChIP-seq 

Copyright (C) 2014-2019 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 21.08.2019
-------------------------

This program is free software.
It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
GNU General Public License for more details.

UCSC reference genome: ftp://hgdownload.soe.ucsc.edu/goldenPath/
************************************************************************************/

#include "isChIP.h"
#include "Imitator.h"

using namespace std;

const string Product::Title = "isChIP";
const string Product::Version = "1.0";
const string Product::Descr = "ChIP-seq simulator";

const string DefFileName [] = { "mTest", "mInput" };
const string OutFileTip = "location of output files or existing folder\n[TEST mode: " + 
	DefFileName[TEST] + ".*, " + 
	"CONTROL mode: " + DefFileName[CONTROL] + ".*]";

// options groups
enum { gTREAT, gFRAG, gREAD, gOUTPUT, gOTHER };	// gOTHER should be the last one

const char* Options::OptGroups [] = {
	"Processing", "Fragment size distribution", "Reads", "Output", "Other"
};
const BYTE Options::GroupCount = sizeof(Options::OptGroups)/sizeof(char*);

// --smode option
const char* smodes	[] = { "SE","PE" };					// corresponds to OutFile::eMode
// --format option: format notations
const char* formats	[] = { "FQ","BED","SAM","WIG","FREQ" };	// corresponds to OutFile::oFormat	
// --verbose option: verbose notations
const char* verbs	[] = { "SL","RES","RT","PAR","DBG" };
// --ground option
const Options::PairVals grounds(100, 1, 0, 0, 100, 100);	// defFG, defBG, minFG, minBG, maxFG, maxFG
// --lndist option
const Options::PairVals lnd(5.46, 0.4, 2, 0.01, 9, 1);	// defMean, defSD, minMean, minSD, maxMean, maxSD
// --ssddist option
const Options::PairVals ssd(vUNDEF, 30, 50, 2, 5000, 500);	// "auto", defSD, minMean, minSD, maxMean, maxSD
// --flat-len option
//const Options::PairVals flattens(0, 0, 0, 0, 50, 50);

const char* ForTest = "For the test mode only";
const string UnstableBSLen = "unstable binding length";

// { char, str, Signs, type, group, defVal, minVal, maxVal, strVal, descr, addDescr }
// defVal: vUNDEF if no default value should be printed
// minVal: vUNDEF if value is prohibited
Options::Option Options::List [] = {
	{ 'g',"gen",	1,	tNAME,	gTREAT, vUNDEF, 0, 0, NULL, "reference genome library.", NULL },
	{ 'n',"cells",	0,	tINT,	gTREAT, 1, 1, 2e6, NULL, "number of nominal cells", NULL },
	{ 'G', "ground",0,	tPAIR_FL,	gTREAT, 0, 0, 0, (char*)&grounds,
	"fore- and background levels:\nnumber of selected fragments inside | outside binding sites,\n\
in percent of all | foreground.\nIn control mode background is ignored ", NULL },
	{ 'D',"mda",	0,	tENUM,	gTREAT,	FALSE,	vUNDEF, 2, NULL, "apply MDA technique", NULL },
	{ 'a',"pcr",	0,	tINT,	gTREAT, 0, 0, 500, NULL, "number of PCR cycles", NULL },
	{ 'c',Chrom::Abbr,0,tNAME,	gTREAT, vUNDEF, 0, 0, NULL,
	"generate output for the specified chromosome only", NULL },
	{ HPH,"bg-all",	0,	tENUM,	gTREAT, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off generation background for all chromosomes.\n", ForTest },
	//{ HPH, "bind-len",	0,	tINT,	gTREAT, 1, 1, 100, NULL, "minimum binding length.", ForTest },
	{ 'm', "smode",		0,	tENUM,	gTREAT, Seq::SE, Seq::SE, Seq::Undef, (char*)smodes,
	"sequencing mode: ? - single end, ? - paired end", NULL },
	{ 's', "bscore",	4,	tINT,	gTREAT, 5, 4, 12, NULL,
	"index of the template field used to score each binding event.\n\
Specify '0' to ignore template scores.", ForTest },
	//{ HPH, "flat-len",	0,	tPAIR_INT,	gTREAT, 0, 0, 0, (char*)&flattens,
	//"inside, outside BS boundary flattening length.", ForTest },
	{ HPH, "edge-len",	0,	tINT,	gTREAT, 0, 0, 10, NULL,
	"unstable binding length (BS edge effect).", ForTest },
	{ HPH,"strand-err",	8,	tFLOAT, gTREAT, 0, 0, 100, NULL,
	"percentage of reads with wrong strand", NULL },
	{ 'N', "full-gen",	0,	tENUM,	gTREAT, FALSE, vUNDEF, 2, NULL,
	"process the entire reference chromosomes (including marginal gaps)", NULL },
	{ 'P', "threads",	0,	tINT,	gTREAT, 1, 1, 20, NULL, "number of threads", NULL },
	{ HPH, "serv",	0,	tNAME,	gTREAT, vUNDEF, 0, 0, NULL,
	"folder to store service files [-g|--gen]", NULL },
	{ HPH, "seed",	0,	tINT,	gTREAT, 0, 0, 1000, NULL,
	"fix random emission with given seed, or 0 if don't fix", NULL },
	{ 'L',"ln",	0,	tPAIR_FL,	gFRAG, 0, 0, 0, (char*)&lnd,
	"mean and stand dev of fragment lognormal distribution", NULL },
	{ 'S',"ss",	2,	tPAIR_INT,	gFRAG, 0, 0, 0, (char*)&ssd,
	"mean and stand dev of size selection normal distribution.\n\
If not specified, then disabled", NULL },
	{ 'r', "rd-len",	0,	tINT,	gREAD, 50, 20, 500, NULL, "length of output read", NULL },
	{ 'p', "rd-pos",	0,	tENUM,	gREAD, FALSE, vUNDEF, 2, NULL,"add read position to its name", NULL },
	{ HPH, "rd-Nlim",	0,	tINT,	gREAD, vUNDEF, 0, 500, NULL,
	"maximum permitted number of ambiguous code N in read [OFF]", NULL },
	{ HPH, "rd-lim",	0,	tLONG,	gREAD, 2e8, 1e5, (float)ULONG_MAX, NULL,
	"maximum permitted number of total recorded reads", NULL },
	{ HPH, "rd-ql",		0,	tCHAR,	gREAD, '~', '!', '~', NULL, "uniform read quality score", NULL },
	{ HPH,"rd-ql-patt",	0,	tNAME,	gREAD, vUNDEF, 0, 0, NULL, "read quality scores pattern", NULL },
	{ HPH,"rd-ql-map",	0,	tINT,	gREAD, 255, 0, 255, NULL,
	"read mapping quality in SAM and BED output", NULL },
	{ 'f',"format",	0,	tCOMB,	gOUTPUT, Output::ofFQ, Output::ofFQ, 5, (char*)formats,
	"format of output data, in any order", NULL },
	{ 'C',"control",0,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL, 
	"generate control simultaneously with test", NULL },
	{ 'x', "strand",0,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL,
	"generate two additional wig files, each one per strand", NULL },
	{ 'T', "sep",	0,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL,
	"display number thousands separator", NULL },
	{ 'o', "output",0,	tNAME,	gOUTPUT, vUNDEF, 0, 0, NULL, OutFileTip.c_str()	},
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
	{ vUNDEF, "<template>", false, "bed file whose features specify binding site"	},
};
const BYTE Options::UsageCount = sizeof(Options::Usages)/sizeof(Options::Usage);

// Returns common name of output files
string GetOutFileName();
void PrintParams(const ChromSizesExt& cSizes, const char* templName, const Features* templ, const Output& oFile);


/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv);
	if( fileInd < 0 )	return 1;		// wrong option or tip output

	int ret = 0;						// main() return code
	Features* templ = NULL;
	const char* fBedName = fileInd==argc ? NULL : argv[fileInd];	// template name

	Chrom::SetCustomOption(oCHROM);
	// initialize Seq before Read::Init()
	Seq::Init(Options::GetIVal(oSMODE), Options::GetFVal(oRD_LIMIT));
	// initialize Read before Output::Init()
	Read::Init(
		Options::GetUIVal(oRD_LEN),
		Options::GetBVal(oRD_NAME),
		Options::GetFVal(oRD_QUAL),
		Options::GetUIVal(oRD_LIMIT_N)
	);
	Output::Init(
		Options::GetIVal(oFORMAT),
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
	DistrParams::Init(lnd.Values(), ssd.Values(), Options::Assigned(oSS_DIST));
	Imitator::Init(
		fBedName != NULL ? TEST : CONTROL,	// fBedName should be set
		Options::GetBVal(oMAKE_INPUT),
		Options::GetFVal(oNUMB_CELLS),
		Options::GetBVal(oMDA),
		Options::GetFVal(oPCR_CYCLES),
		Options::GetUIVal(oVERB),
		Options::GetBVal(oBG_ALL) * bool(SAMPLE_BG()),	// bg all: always false if bg level is 0
		!Options::GetBVal(oBS_SCORE),
		//readlen(Options::GetUIVal(oBS_LEN)),
		//flattens.Values()
		Options::GetUIVal(oFLAT_LEN)
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

		ChromSizesExt cSizes(
			Options::GetSVal(oGEN), Options::GetSVal(oSERV), Imitator::Verbose(vRT) );

		if(fBedName) {
			//Obj::Ambig::SetSupplAlarm(Obj::Ambig::SHORT, "for given binding length");
			templ = new Features(sTemplate, fBedName, cSizes, 
				Imitator::Verbose(vDBG) ? Obj::iSTAT : Obj::iNONE, 
				false, Options::GetIVal(oBS_SCORE),
				0, //readlen(Options::GetUIVal(oBS_LEN)),	// binding length
				Imitator::Verbose(vDBG)								// print ambigs alarm
			);
			//if(flattens.Values().second)
			//	templ->Extend(flattens.Values().second, cSizes, Obj::iLAC);
			chrlen halfMinFLen = templ->GetMinFeatureLength()/2;
			if(halfMinFLen < Imitator::FlatLen) {
				Imitator::FlatLen = halfMinFLen;
				if( Imitator::Verbose(vRT) )
					Err(UnstableBSLen + " reduced to " + NSTR(halfMinFLen)
					+ " due to minimum template feature length").Warning();

			}
		}
		// set treated chroms
		if( !cSizes.SetTreated(Imitator::All, templ) )
			Err(Err::TF_EMPTY, fBedName, "features per stated " + Chrom::ShortName(Chrom::CustomID()))
				.Throw();
		
		Imitator::SetThreadNumb( min(chrid(Options::GetFVal(oNUMB_THREAD)), cSizes.TreatedCount()) );

		Output oFile(GetOutFileName(), Imitator::IsControl(),
			cSizes, qlPatt, Options::CommandLine(argc, argv));

		PrintParams(cSizes, fBedName, templ, oFile);
		Imitator(cSizes, oFile).Execute(templ);
	}
	catch(Err &e)				{ ret = 1; cerr << e.what() << endl; }
	catch(const exception &e)	{ ret = 1; cerr << e.what() << EOL; }
	//catch(...)					{ ret = 1; cerr << "Unregistered error" << endl; }
	if(templ)	delete templ;
	Mutex::Finalize();
	if(!ret)
		Timer::StopCPU(false),
		timer.Stop("\twall-clock: ", false, true);
	//MemStatus::StopObserve();
	return ret;
}


void PrintParams(const ChromSizesExt& cSizes, const char* templName, 
	const Features* templ, const Output& oFile)
{
	if( !Imitator::Verbose(vPAR) )	return;
	cout << SignPar << "Reference" << SepDCl << "genome" << SepCl << cSizes.RefPath();
	cout << SepCm << Chrom::TitleName() << SepCl;
	if(!Chrom::NoCustom())		
		cout << Chrom::Mark(Chrom::CustomID());
	else if(cSizes.TreatedCount() < cSizes.ChromCount())	
		cSizes.PrintTreatedChroms();
	else	
		cout << "all";
 	cout << EOL;
	if(!cSizes.IsServAsRef())
		cout << SignPar << "Service folder" << SepCl << cSizes.ServPath() << EOL;

	if(templName) {
		cout << SignPar << sTemplate << SepCl << templName << SepCl;
		templ->PrintItemCount();
	}
	oFile.PrintFormat(SignPar);		// print output formats, sequencing mode
	cout << SignPar; Seq::Print();
	cout << SignPar << "Count of cells" << SepCl << ULONG(Options::GetFVal(oNUMB_CELLS)) << EOL;
	Imitator::PrintAmpl();
	cout << SignPar;
	Read::Print();
	oFile.PrintReadQual(SignPar);
	cout << SignPar << "Optimization: process the entire ref. chromosome"
		 << SepCl << Options::BoolToStr(oLET_GAPS) << EOL;
	cout << SignPar << "Stated sample: ";
	if(TestMode) {
		cout << "foreground" << Equel << SAMPLE_FG()
				<< PERS << SepSCl << "background" << Equel << SAMPLE_BG() << PERS << EOL;
		if(cSizes.TreatedCount() > 1) {
			cout << SignPar << "Background for all chromosomes" << SepCl;
			if(bool(SAMPLE_BG()))	cout << Options::BoolToStr(oBG_ALL) << EOL;
			else					cout << "negligible due to zero background sample\n";
		}
		//cout << SignPar << "Binding length" << SepCl << Options::GetUIVal(oBS_LEN) << EOL;
		//if(flattens.Values().first || flattens.Values().second)
		//	cout << SignPar << "BS boundary flattening length: inside" << Equel << flattens.Values().first
		//		 << SepSCl << "outside" << Equel << flattens.Values().second << EOL;
		if(Imitator::FlatLen)
			cout << SignPar << UnstableBSLen << SepCl << Imitator::FlatLen << EOL;
		cout << SignPar << "Score index in template"  << SepCl;
		if(Imitator::UniformScore)	cout << "uniform score\n";
		else	cout << Options::GetIVal(oBS_SCORE) << EOL;
	}
	else	cout << SAMPLE_FG() << PERS << EOL;		// control mode
	float ws = Options::GetFVal(oSTRAND_ERR);
	if(ws)	cout << SignPar << "Reads with wrong strand"
					<< SepCl << ws << PERS << EOL;
	if(!Imitator::IsSingleThread())
		cout << SignPar << "Actual threads" << SepCl << int(Imitator::ThrCnt) << EOL;
	cout << SignPar << "Fragment lognorm distribution"
			<< SepCl << "mean" << Equel << DistrParams::lnMean 
			<< SepSCl << "stand.dev" << Equel << DistrParams::lnSigma
			<< SepSCl << "Mean" << Equel << setprecision(5) << DistrParams::LnMean() 
			<< SepSCl << "Mode" << Equel << DistrParams::LnMode() << EOL;
	cout << SignPar << "Fragment size selection" << SepCl;
	if(DistrParams::IsSS())
		cout << "mean" << Equel << setprecision(5) << DistrParams::ssMean 
				<< SepSCl << "stand.dev" << Equel << DistrParams::ssSigma << EOL;
	else
		cout << Options::BoolToStr(false) << EOL;
}

// Returns common name of output files
string GetOutFileName()
{
	const char* outName = Options::GetSVal(oOUT_FILE);
	
	if( !outName )
		return DefFileName[Imitator::TMode];
	if( FS::IsDirExist(outName) )
		return FS::MakePath(string(outName)) + DefFileName[Imitator::TMode];
	return string(outName);
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
