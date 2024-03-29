/**********************************************************
isChIP.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 02.01.2022
-------------------------
Provides option emum
***********************************************************/
#pragma once
#include "def.h"

#define	SAMPLE_FG()		grounds.Values().first
#define	SAMPLE_BG()		grounds.Values().second

enum optValue {
	oGEN,
	//oPROT,
	oNUMB_CELLS,
	oGR_LEVEL,
	oEXO,
	oMDA,
	oPCR_CYCLES,
	oCHROM,
	oBG_ALL,
	//oBS_LEN,
	oSMODE,
	oFLAT_LEN,
	oSTRAND_ERR,
	oLET_GAPS,
	oNUMB_THREAD,
	oSERV,
	oSEED,
	oOVERL,
	oBS_SCORE,
	oFR_DIST,
	oSS_DIST,
	oRD_LEN,
	oRD_DIST,
	oRD_NAME,
	oRD_LIMIT_N,
	oRD_LIMIT,
	oRD_QUAL,
	oRD_QUAL_PATT,
	oMAP_QUAL,
	oFORMAT,
	oMAKE_INPUT,
	oSTRAND,
	oOUT_FILE,
	oLOCALE,
#ifndef _NO_ZLIB
	oGZIP,
#endif
	oTIME,
	oVERB,
	oVERSION,
	oHELP
};