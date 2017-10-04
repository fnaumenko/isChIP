#pragma once

#define V_CRIT	0
#define V_RES	1
#define V_RT	2
#define V_PAR	3
#define V_DEBUG	4
#define	RGL_SHIFT()		short(Options::GetDVal(oREG_MODE))
#define	THREADS_CNT()	BYTE(Options::GetDVal(oNUMB_THREAD))
#define	SAMPLE_FG()		Options::GetDVal(oFG_VEVEL)
#define	SAMPLE_BG()		Options::GetDVal(oBG_LEVEL)
#define	SZ_SEL()		Options::GetBVal(oSZ_SEL)

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
	oMAP_QUAL,
	oSZ_SEL,
	oSZ_SEL_SIGMA,
	//oSTDOUT,
	oSMODE,
	oSTRAND_MIX,
	oTS_UNIFORM,
	oNUMB_THREAD,
	oDBG,
	oREG_MODE,
	oFORMAT,
	oOUT_FILE,
	oGZIP,
	oTIME,
	oVERB,
	oVERSION,
	oHELP
};

