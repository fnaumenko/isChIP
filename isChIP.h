#pragma once
#include "def.h"

#define	RGL_SHIFT()		short(Options::GetDVal(oREG_MODE))
#define	THREADS_CNT()	BYTE(Options::GetDVal(oNUMB_THREAD))
#define	SAMPLE_FG()		Options::GetDVal(oFG_VEVEL)
#define	SAMPLE_BG()		Options::GetDVal(oBG_LEVEL)
#define	SZ_SEL()		Options::GetBVal(oSZ_SEL)

enum eVerb {	// verbose level
	vCRIT,		// print critical messages only
	vRES,		// print vCRIT + results
	vRT,		// print vRES + runtime info
	vPAR,		// print vRT + parameters
	vDEBUG		// print vPAR + additional info
};

enum optValue {
	oAMPL,
	oBG_LEVEL,
	oFG_VEVEL,
	oNUMB_CELLS,
	oGFILE,
	oCHROM,
	oFRAG_LEN,
	oFRAG_DEV,
	oBG_ALL,
	oBS_LEN,
	oFLAT_LEN,
	oMEAN,
	oSIGMA,
	oLN_FACTOR,
	oLN_TERM,
	oLET_N,
	oREAD_LEN,
	oREAD_NAME,
	oREAD_LIMIT_N,
	oREAD_LIMIT,
	oFQ_QUAL,
	oFQ_QUAL_PATT,
	oMAP_QUAL,
	oSZ_SEL,
	oSZ_SEL_SIGMA,
	//oSTDOUT,
	oSMODE,
	oSTRAND_MIX,
	oTS_UNIFORM,
	oNUMB_THREAD,
	oFIX,
	oREG_MODE,
	oFORMAT,
	oOUT_FILE,
#ifndef _NO_ZLIB
	oGZIP,
#endif
	oTIME,
	oVERB,
	oVERSION,
	oHELP
};

