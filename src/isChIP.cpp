/************************************************************************************
isChIP (In-Silico ChIP) is a fast realistic ChIP-seq simulator.
The model is based on the real protocol of ChIP-seq

Copyright (C) 2014-2022 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 02.01.2022
-------------------------

This program is free software.
It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
GNU General Public License for more details.

UCSC reference genome:		ftp://hgdownload.soe.ucsc.edu/goldenPath/
Ensembl reference genome:	ftp://ftp.ensembl.org/pub/release-100/fasta/
Ensembl hg38 genome:	ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/
************************************************************************************/

#include "isChIP.h"
#include "Imitator.h"

using namespace std;

const string Product::Title = "isChIP";
const string Product::Version = "2.0";
const string Product::Descr = "ChIP-seq simulator";

const string DefFileName[] = { "mTest", "mInput" };
const string OutFileTip = "location of output files or existing folder\n[TEST mode: " +
	DefFileName[TEST] + ".*, " +
	"CONTROL mode: " + DefFileName[CONTROL] + ".*]";

// --smode option
const char* smodes[] = { "SE","PE" };						// corresponds to OutFile::eMode
// --format option: format notations
const char* formats[] = { "FQ","BED","SAM","BG","FDENS","RDENS","FDIST","RDIST" };	// corresponds to Output::oFormat	
// --verbose option: verbose notations
const char* verbs[] = { "SL","RES","RT","PAR","DBG" };
// --ground option
const Options::PairVals grounds(50, 1, 0, 0, 100, 100);		// defFG, defBG, minFG, minBG, maxFG, maxFG
// --lndist option
const Options::PairVals lnd(5.26f, 0.3f, 3, 0.01f, 9, 1);	// defMean, defSD, minMean, minSD, maxMean, maxSD
// --ssddist option
const Options::PairVals ssd(vUNDEF, 30, 50, 2, 2000, 500);	// "auto", defSD, minMean, minSD, maxMean, maxSD
// --rd-dist option
const Options::PairVals rdd(200, 20, 50, 2, 1000, 300);		// defMean, defSD, minMean, minSD, maxMean, maxSD
// --flat-len option
//const Options::PairVals flattens(0, 0, 0, 0, 50, 50);

const char* ForTest = "For the test mode only";
const string UnstableBSLen = "unstable binding length";

//using flags = Options::Signs::oSign;
//typedef Options::Signs::Flags flag;

// *** Options definition

enum { gTREAT, gTEMPL, gFRAG, gREAD, gOUTPUT, gOTHER };	// gOTHER should be the last one
const char* Options::OptGroups[] = {
	"Processing", "Template", "Fragment size distribution", "Reads", "Output", "Other"
};
const BYTE Options::GroupCount = ArrCnt(Options::OptGroups);

// { char, str, Signs, type, group, defVal, minVal, maxVal, strVal, descr, addDescr }
// defVal: vUNDEF if no default value should be printed
// minVal: vUNDEF if value is prohibited
Options::Option Options::List[] = {
	{ 'g',	sGen,	fOblig,	tNAME,	gTREAT, vUNDEF, 0, 0, NULL, "reference genome library.", NULL },
	//{ 'p',"prot",	fNone,	tENUM,	gTREAT, 0, 0, 1, (char*)prots,
	//"protocol: ? - Illumina, ? - Ion Torent", NULL },
	{ 'n',"cells",	fNone,	tINT,	gTREAT, 1, 1, 2e6, NULL, "number of nominal cells", NULL },
	{ 'G',"ground",	fNone,	tPR_FL,	gTREAT, 0, 0, 0, (char*)&grounds,
	"fore- and background levels:\nnumber of selected fragments inside/outside binding sites,\n\
in percent of all/foreground.\nIn control mode background is ignored ", NULL },
	{ 'E',"exo",	fOptnal,tINT,	gTREAT,	6,	0, 20, NULL,
	"apply ChIP-exo protocol with specified exonuclease 'headroom' length.\n\
If the option is not specified, ChIP-exo is not applied", NULL },
	{ 'D',"mda",	fNone,	tENUM,	gTREAT,	FALSE,	vUNDEF, 2, NULL, "apply MDA technique", NULL },
	{ 'a',"pcr",	fNone,	tINT,	gTREAT, 0, 0, 500, NULL, "number of PCR cycles", NULL },
	{ 'c',Chrom::Abbr,fNone,tNAME,	gTREAT, vUNDEF, 0, 0, NULL,
	"generate output for the specified chromosome only", NULL },
	{ HPH,"bg-all",	fNone,	tENUM,	gTREAT, TRUE, 0, 2, (char*)Options::Booleans,
	"turn on/off generation background for all chromosomes.\n", ForTest },
	//{ HPH, "bind-len", fNone,	tINT,	gTREAT, 1, 1, 100, NULL, "minimum binding length.", ForTest },
	{ 'm', "smode",	fNone,	tENUM,	gTREAT, Seq::SE, Seq::SE, ArrCnt(smodes), (char*)smodes,
	"sequencing mode: ? - single end, ? - paired end", NULL },
	//{ HPH, "flat-len",	fNone,	tPR_INT,	gTREAT, 0, 0, 0, (char*)&flattens,
	//"inside, outside BS boundary flattening length.", ForTest },
	{ HPH, "edge-len",	fNone,	tINT,	gTREAT, 0, 0, 10, NULL,
	"unstable binding length (BS edge effect).", ForTest },
	{ HPH,"strand-err",	fHidden,tFLOAT, gTREAT, 0, 0, 100, NULL,
	"percentage of reads with wrong strand", NULL },
	{ 'N', "full-gen",	fNone,	tENUM,	gTREAT, FALSE, vUNDEF, 2, NULL,
	"process the entire reference chromosomes (including marginal gaps)", NULL },
	{ 'P',"threads",fNone,	tINT,	gTREAT, 1, 1, 20, NULL, "number of threads", NULL },
	{ HPH, "serv",	fNone,	tNAME,	gTREAT, vUNDEF, 0, 0, NULL,
	"folder to store service files [-g|--gen]", NULL },
	{ HPH, "seed",	fNone,	tINT,	gTREAT, 0, 0, 1000, NULL,
	"fix random emission with given seed, or 0 if don't fix", NULL },
	{ 'o', "overl",	fNone,	tENUM,	gTEMPL, FALSE,	0, 2, (char*)Options::Booleans,
	"allow (and merge) overlapping template features", NULL },
	{ 's',"bscore",	fAllow0,tINT,	gTEMPL, 5, 4, 12, NULL,
	"index of the template field used to score each binding event.\n\
Value '0' means ignore template scores.", ForTest },
	{ 'L',"ln",	fNone,	tPR_FL,	gFRAG, 0, 0, 0, (char*)&lnd,
	"mean and stand dev of fragment lognormal distribution", NULL },
	{ 'S',"ss",	fOptnal,tPR_INT,gFRAG, 0, 0, 0, (char*)&ssd,
	"apply size selection normal distribution with specified mean and stand dev.\n\
If the option is not specified, size selection is disabled", NULL },
	{ 'r',"rd-len",	fNone,	tINT,	gREAD, 50, 20, 1000, NULL, 
	"fixed length of output read, or minimum length of variable reads", NULL },
	{ 'R',"rd-dist",fOptnal,tPR_INT,gREAD, 0, 0, 0, (char*)&rdd,
	"mean and stand dev of variable read normal distribution,\naccording to Ion Torrent/Roche454 protocol.\n\
If the option is not specified, the read length is fixed", NULL },
	{ 'l',"rd-pos",	fNone,	tENUM,	gREAD, FALSE, vUNDEF, 2, NULL,"add read position (location) to its name", NULL },
	{ HPH,"rd-Nlim",fNone,	tINT,	gREAD, vUNDEF, 0, 500, NULL,
	"maximum permitted number of ambiguous code N in read [OFF]", NULL },
	{ HPH,"rd-lim",	fNone,	tLONG,	gREAD, 2e8, 1e5, (float)ULONG_MAX, NULL,
	"maximum permitted number of total recorded reads", NULL },
	{ HPH,"rd-ql",	fNone,	tCHAR,	gREAD, '~', '!', '~', NULL, "uniform read quality score", NULL },
	{ HPH,"rd-ql-patt",	fNone,	tNAME,	gREAD, vUNDEF, 0, 0, NULL, "read quality scores pattern", NULL },
	{ HPH,"rd-ql-map",	fNone,	tINT,	gREAD, 255, 0, 255, NULL,
	"read mapping quality in SAM and BED output", NULL },
	{ 'f',"format",	fNone,	tCOMB,	gOUTPUT, float(Output::eFormat::FG), float(Output::eFormat::FG), ArrCnt(formats),
	(char*)formats,	"format of output data, in any order", NULL },
	{ 'C',"control",fNone,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL,
	"generate control simultaneously with test", NULL },
	{ 'x',"strand",	fNone,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL,
	"generate two additional wig files, each one per strand", NULL },
	{ 'O', sOutput,	fNone,	tNAME,	gOUTPUT, vUNDEF, 0, 0, NULL, OutFileTip.c_str()	},
	{ 'T', "sep",	fNone,	tENUM,	gOUTPUT, FALSE,	vUNDEF, 2, NULL, "use 1000 separator in output", NULL },
#ifndef _NO_ZLIB
	{ 'z',"gzip",	fNone,	tENUM,	gOUTPUT, FALSE, vUNDEF, 2, NULL, "compress the output", NULL},
#endif
	{ 't',	sTime,	fNone,	tENUM,	gOTHER,	FALSE,	vUNDEF, 2, NULL, sPrTime, NULL },
	{ 'V',"verbose",fNone,	tENUM,	gOTHER, float(eVerb::PAR), float(eVerb::CRIT), ArrCnt(verbs), (char*)verbs,
	"set verbose level:\n? -\tsilent mode (show critical messages only)\n? -\tshow result summary\n?  -\tshow run-time information\n? -\tshow actual parameters\n? -\tshow debug messages", NULL },
	{ 'v',	sVers,	fNone,	tVERS,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrVersion, NULL },
	{ 'h',	sHelp,	fNone,	tHELP,	gOTHER,	vUNDEF, vUNDEF, 0, NULL, sPrUsage, NULL }
};
const BYTE	Options::OptCount = ArrCnt(Options::List);

const Options::Usage Options::Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, "<template>", false, "bed file whose features specify binding site"	},
};
const BYTE Options::UsageCount = ArrCnt(Options::Usages);

// Returns common name of output files
string GetOutFileName();
void PrintParams(const ChromSizesExt& cSizes, const char* templName, const Features* templ, const Output& oFile);

/*****************************************/
int main(int argc, char* argv[])
{
	int fileInd = Options::Parse(argc, argv);
	if (fileInd < 0)	return 1;		// wrong option or tip output

	int ret = 0;						// main() return code
	Features* templ = NULL;
	const char* fBedName = fileInd == argc ? NULL : argv[fileInd];	// template name

	// initialize Seq before Read::Init()
	Seq::Init(Options::GetIVal(oSMODE), ULONG(Options::GetFVal(oRD_LIMIT)));
	// initialize Read before Output::Init()
	Read::Init(
		Options::GetUIVal(oRD_LEN),
		Options::GetBVal(oRD_NAME),
		char(Options::GetFVal(oRD_QUAL)),
		Options::GetUIVal(oRD_LIMIT_N)
	);
	// initialize DistrParams before Imitator and Output!!
	DistrParams::Init(
		lnd.Values(), Options::Assigned(oFR_DIST),
		ssd.Values(), Options::Assigned(oSS_DIST),
		rdd.Values(), Options::Assigned(oRD_DIST)
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
	Imitator::Init(
		fBedName != NULL ? TEST : CONTROL,	// fBedName should be set
		Options::GetBVal(oMAKE_INPUT),
		ULONG(Options::GetFVal(oNUMB_CELLS)),
		Options::Assigned(oEXO),
		Options::Assigned(oRD_LEN),
		Options::GetBVal(oMDA),
		a_coeff(Options::GetFVal(oPCR_CYCLES)),
		Options::GetUIVal(oVERB),
		Options::GetBVal(oBG_ALL) * bool(SAMPLE_BG()),	// bg all: always false if bg level is 0
		//!Options::GetBVal(oBS_SCORE),
		//readlen(Options::GetUIVal(oBS_LEN)),
		//flattens.Values()
		Options::GetUIVal(oFLAT_LEN)
	);
	RefSeq::LetGaps = Options::GetBVal(oLET_GAPS);
	RefSeq::StatGaps = Imitator::Verbose(eVerb::PAR);		// Imitator::Init() should be called before
	Random::SetSeed(Options::GetUIVal(oSEED), Options::GetUIVal(oEXO));
	if (Options::GetBVal(oLOCALE))	cout.imbue(locale(LOCALE_ENG));

	// execution
	Mutex::Init(Options::GetIVal(oNUMB_THREAD) > 1);
	Timer::Enabled = Options::GetBVal(oTIME) && Imitator::Verbose(eVerb::RES);
	Timer::StartCPU();
	Timer timer;
	try {
		// check file names first of all
		FS::CheckedFileName(fBedName);
		Output::SetReadQualPatt(FS::CheckedFileName(oRD_QUAL_PATT));	// read quality pattern file name

		ChromSizesExt cSizes(
			Options::GetSVal(oGEN), oCHROM, Imitator::Verbose(eVerb::RT), Options::GetSVal(oSERV));

		if (fBedName) {
			//Obj::Ambig::SetSupplAlarm(Obj::Ambig::SHORT, "for given binding Imitator::FlatLenlength");
			templ = new Features(fBedName, cSizes,
				Options::GetBVal(oOVERL),
				Options::GetIVal(oBS_SCORE),
				0, //readlen(Options::GetUIVal(oBS_LEN)),		// binding length
				true											// print name
				//Imitator::Verbose(eVerb::DBG)						// print ambigs alarm
			);
			//templ->Print();
			//if(flattens.Values().second)
			//	templ->Extend(flattens.Values().second, cSizes, Obj::iLAC);
			//templ->Extend(10000, cSizes, Obj::iSTAT);

			chrlen halfMinFLen = templ->GetMinFeatureLength() / 2;
			if (halfMinFLen < chrlen(Imitator::FlatLen)) {
				Imitator::FlatLen = halfMinFLen;
				if (Imitator::Verbose(eVerb::RT))
					Err(UnstableBSLen + " reduced to " + to_string(halfMinFLen)
						+ " due to minimum template feature length").Warning();
			}
		}

		// set treated chroms
		if (!cSizes.SetTreated(Imitator::All, templ))
			Err(Err::TF_EMPTY, fBedName, "features per stated " + Chrom::ShortName(Chrom::CustomID()))
			.Throw();

		Imitator::SetThreadNumb(min(chrid(Options::GetFVal(oNUMB_THREAD)), cSizes.TreatedCount()));
		Output oFile(GetOutFileName(), Imitator::IsControl(), Options::CommandLine(argc, argv), cSizes);

		PrintParams(cSizes, fBedName, templ, oFile);
		Imitator(cSizes, oFile).Execute(templ);
	}
	catch (Err & e) { ret = 1; cerr << e.what() << LF; }
	catch (const exception & e) { ret = 1; cerr << e.what() << LF; }
	//catch(...)					{ ret = 1; cerr << "Unregistered error" << endl; }
	if (templ)	delete templ;
	//Mutex::Finalize();
	if (!ret)
		Timer::StopCPU(false),
		timer.Stop("\twall-clock: ", false, true);
	//MemStatus::StopObserve();
	return ret;
}

// Prints simulation parameters in -V par mode
void PrintParams(const ChromSizesExt& cSizes, const char* templName,
	const Features* templ, const Output& oFile)
{
	if (!Imitator::Verbose(eVerb::PAR))		return;
	cout << SignPar << "Reference" << SepDCl << "genome" << SepCl << cSizes.RefPath();
	cout << SepCm << Chrom::TitleName() << COLON;
	if (!Chrom::IsCustom())
		cout << SPACE << Chrom::Mark(Chrom::CustomID());
	else 
		cSizes.PrintTreatedChroms();
	cout << LF;

	if (!cSizes.IsServAsRef())
		cout << SignPar << "Service folder" << SepCl << cSizes.ServPath() << LF;

	if (templName) {
		cout << SignPar << sTemplate << SepCl << templName << SepCl;
		templ->PrintItemCount(FT::eType::BED, false);
		cout << SepSCl;
		if (templ->IsUniScore())		cout << "uniform score\n";
		else	cout << "score index" << Equel << Options::GetIVal(oBS_SCORE) << LF;
	}
	oFile.PrintFormat(SignPar);		// print output formats, sequencing mode
	Seq::Print(SignPar);			// print sequencing modes
	cout << SignPar << "Sequencing modification" << SepCl << "ChIP-";
	if (Imitator::IsExo)
		cout << "exo" << SepSCl << "exonuclease 'headroom' length" << Equel << Options::GetUIVal(oEXO) << LF;
	else
		cout << "seq\n";
	cout << SignPar << "Count of cells" << SepCl << ULONG(Options::GetFVal(oNUMB_CELLS)) << LF;
	Imitator::PrintAmpl	(SignPar);
	Read::PrintParams	(SignPar, DistrParams::IsRVL());
	oFile.PrintReadQual	(SignPar);
	DistrParams::PrintReadDistr(cout, SignPar, Read::Title);
	cout << SignPar << "Optimization: process the entire ref. " << Chrom::Title()
		<< SepCl << Options::BoolToStr(oLET_GAPS) << LF;
	cout << SignPar << "Stated sample: ";
	if (TestMode) {
		cout << "foreground" << Equel << SAMPLE_FG()
			<< PERS << SepSCl << "background" << Equel << SAMPLE_BG() << PERS
			<< " (relative to the foreground)\n";
		if (cSizes.TreatedCount() > 1) {
			cout << SignPar << "Background for all " << Chrom::Title(true) << SepCl;
			if (bool(SAMPLE_BG()))	cout << Options::BoolToStr(oBG_ALL) << LF;
			else					cout << "negligible due to zero background sample\n";
		}
		//cout << SignPar << "Binding length" << SepCl << Options::GetUIVal(oBS_LEN) << LF;
		//if(flattens.Values().first || flattens.Values().second)
		//	cout << SignPar << "BS boundary flattening length: inside" << Equel << flattens.Values().first
		//		 << SepSCl << "outside" << Equel << flattens.Values().second << LF;
		if (Imitator::FlatLen)
			cout << SignPar << UnstableBSLen << SepCl << Imitator::FlatLen << LF;
	}
	else	cout << SAMPLE_FG() << PERS << LF;		// control mode
	float ws = Options::GetFVal(oSTRAND_ERR);
	if (ws)	cout << SignPar << Read::Title << "s with wrong strand"
		<< SepCl << ws << PERS << LF;
	if (!Imitator::IsSingleThread())
		cout << SignPar << "Actual threads" << SepCl << int(Imitator::ThrCnt) << LF;
	DistrParams::PrintFragDistr(cout, SignPar, true);
}

// Returns common name of output files
string GetOutFileName()
{
	const char* outName = Options::GetSVal(oOUT_FILE);

	if (!outName)
		return DefFileName[Imitator::TMode];
	if (FS::IsDirExist(outName))
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
