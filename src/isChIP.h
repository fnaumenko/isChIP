/**********************************************************
isChIP.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 05/01/2024
-------------------------
Provides option emum
***********************************************************/
#pragma once

static const char* Equel = " = ";

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
	oOUTFILE,
	oLOCALE,
#ifdef _ZLIB
	oGZIP,
#endif
	oTIME,
	oVERB,
	oVERSION,
	oHELP
};