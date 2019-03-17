#pragma once
#include "def.h"

#define	SAMPLE_FG()		Options::GetFVal(oFG_VEVEL)
#define	SAMPLE_BG()		Options::GetFVal(oBG_LEVEL)

enum eVerb {	// verbose level
	vCRIT,		// print critical messages only
	vRES,		// print vCRIT + results
	vRT,		// print vRES + runtime info
	vPAR,		// print vRT + parameters
	vDBG		// print vPAR + additional info
};

enum optValue {
	oGFILE,
	oNUMB_CELLS,
	oAMPL,
	oAMPL_COEF,
	oBG_LEVEL,
	oFG_VEVEL,
	oBG_ALL,
	oCHROM,
	oFRAG_MEAN,
	oFRAG_SIGMA,
	oSS_MEAN,
	oSS_SIGMA,
	oSS,
	oLET_N,
	//oBS_LEN,
	oFLAT_LEN,
	oREAD_LEN,
	oREAD_NAME,
	oREAD_LIMIT_N,
	oREAD_LIMIT,
	oFQ_QUAL,
	oFQ_QUAL_PATT,
	oMAP_QUAL,
	//oSTDOUT,
	oSMODE,
	oSTRAND_ERR,
	oTS_UNIFORM,
	oNUMB_THREAD,
	oSEED,
	oFORMAT,
	oMAKE_INPUT,
	oSTRAND,
	oLOCALE,
	oOUT_FILE,
#ifndef _NO_ZLIB
	oGZIP,
#endif
	oTIME,
	oVERB,
	oVERSION,
	oHELP
};

