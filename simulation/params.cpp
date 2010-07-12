/*
 * params.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "params.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

Params* Params::_pInstance = NULL;

const char* Params::_description[][4] =
{
	{ "diffusivityTNF", "TNF diffusivity", "cm^2/s", "GR" },
	{ "diffusivityShedTNFR2", "Shed sTNF/TNFR2 comples diffusivity", "cm^2/s", "GR" },
	{ "diffusivityChemokines", "Chemokine diffusivity", "cm^2/s", "GR" },
	{ "degRateTNF", "TNF degradation rate (per 6s)", "", "GR" },
	{ "degRateChemokines", "Chemokine degradation rate (per 6s)", "", "GR" },
	{ "minChemotaxis", "Chemotaxis sensivity range (lower bound)", "#molecules", "GR" },
	{ "maxChemotaxis", "Chemotaxis sensivity range (upper bound)", "#molecules", "GR" },
	// molecular TNF-associated parameters
	{ "kSynthMac", "Rate of mTNF synthesis by a macrophage", "#/cell.sec", "GR" },
	{ "kSynthTcell", "Rate of mTNF synthesis by a T cell", "#/cell.sec", "GR" },
	{ "kTaceMac", "Rate of mTNF release from a macrophage", "1/sec", "GR" },
	{ "kTaceTcell", "Rate of mTNF release from a T cell", "1/sec", "GR" },
	{ "KD1", "sTNF/TNFR1 dissociation equilibrium constant", "M", "GR" },
	{ "KD2", "sTNF/TNFR2 dissociation equilibrium constant", "M", "GR" },
	{ "kOn1", "sTNF/TNFR1 association rate constant", "1/M.s", "GR" },
	{ "kOn2", "sTNF/TNFR2 association rate constant", "1/M.s", "GR" },
	{ "kInt1", "sTNF/TNFR1 complex internaliation rate constant", "1/s", "GR" },
	{ "kInt2", "sTNF/TNFR2 complex internaliation rate constant", "1/s", "GR" },
	{ "kShed", "sTNF/TNFR2 complex shedding rate constant", "1/s", "GR" },
	{ "kRec1", "TNFR1 recycling rate constant", "1/s", "GR" },
	{ "kRec2", "TNFR2 recycling rate constant", "1/s", "GR" },
	{ "kT1", "TNFR1 turn-over rate constant", "1/s", "GR" },
	{ "kT2", "TNFR2 turn-over rate constant", "1/s", "GR" },
	{ "kDeg1", "TNFR1 degradtion rate constant", "1/s", "GR" },
	{ "kDeg2", "TNFR2 degradtion rate constant", "1/s", "GR" },
	{ "maxTNFR1Mac", "Maximum density of TNFR1 on macrophages", "#/cell", "GR" },
	{ "minTNFR1Mac", "Minimum density of TNFR1 on macrophages", "#/cell", "GR" },
	{ "maxTNFR2Mac", "Maximum density of TNFR2 on macrophages", "#/cell", "GR" },
	{ "minTNFR2Mac", "Minimum density of TNFR2 on macrophages", "#/cell", "GR" },
	{ "maxTNFR1Tcell", "Maximum density of TNFR1 on T cells", "#/cell", "GR" },
	{ "minTNFR1Tcell", "Minimum density of TNFR1 on T cells", "#/cell", "GR" },
	{ "maxTNFR2Tcell", "Maximum density of TNFR2 on T cells", "#/cell", "GR" },
	{ "minTNFR2Tcell", "Minimum density of TNFR2 on T cells", "#/cell", "GR" },
	// end of molecular TNF-associated parameters
	{ "thresholdApoptosisTNF", "TNF threshold for TNF-induced apoptosis", "#molecules", "GR" },
	{ "kApoptosis", "Rate of apoptosis happening", "1/s", "GR" },
	{ "probApoptosisTNF", "Probability of TNF-induced apoptosis", "", "GR" },
	{ "effectRecTNF", "Effect of TNF on recruitment", "", "GR" },
	{ "effectRecCCL2", "Effect of CCL2 on recruitment", "", "GR" },
	{ "effectRecCCL5", "Effect of CCL5 on recruitment", "", "GR" },
	{ "effectRecCXCL9", "Effect of CXCL9 on recruitment", "", "GR" },
	{ "dAttractant", "Secretion rate of macrophage attractant secreted by debris", "", "GR" },
	{ "dCCL2", "Secretion rate of CCL2", "#molecules/6s", "Mac" },
	{ "dCCL5", "Secretion rate of CCL5", "#molecules/6s", "Mac" },
	{ "dCXCL9", "Secretion rate of CXCL9", "#molecules/6s", "Mac" },
	{ "dTNF", "Secretion rate of TNF", "#molecules/6s", "Mac" },
	{ "thresholdNFkBTNF", "TNF threshold for NFkB activation", "#molecules", "Mac" },
	{ "kNFkB", "Rate of NFkB activation", "1/s", "Mac" },
	{ "nrExtMtbUptakeRest", "Number of extracellular bacteria a resting macrophage can uptake/kill", "#bacteria", "Mac" },
	{ "probKillExtMtbRest", "Probability of a resting macrophage to kill some extracellular bacteria", "", "Mac" },
	{ "nrExtMtbNFkB", "Number of extracellular bacteria for NFkB activation in an infected macrophage", "#bacteria", "Mac" },
	{ "nrIntMtbCInf", "Number of intracellular bacteria for an infected macrophage to become chronically infected", "#bacteria", "Mac" },
	{ "nrIntMtbBurstCInf", "Number of intracellular bacteria for an chronically infected macrophage to burst", "#bacteria", "Mac" },
	{ "probStat1Tgam", "Probability of STAT1 activation due to a Tgam cell", "", "Mac" },
	{ "nrExtMtbUptakeAct", "Number of extracellular bacteria an active macrophage can uptake and subsequently kill", "#bacteria", "Mac" },
	{ "thresholdRec", "TNF/chemokine threshold for macrophage recruitment", "", "Mac" },
	{ "probRec", "Probability of recruiting a macrophage", "", "Mac" },
	{ "probRec", "Probability of recruiting a T cell", "", "Tcell" },
	{ "probMoveToMac", "Probability of a T cell moving onto a compartment already containing a macrophage", "", "Tcell" },
	{ "probMoveToTcell", "Probability of a T cell moving onto a compartment already containing another T cell", "", "Tcell" },
	{ "probApoptosisFasFasL", "Probability of Fas/FasL induced apoptosis by a Tgam cell", "", "Tgam" },
	{ "thresholdRec", "TNF/chemokine threshold for Tgam recruitment", "", "Tgam" },
	{ "probRec", "Probability of recruiting a Tgam cell", "", "Tgam" },
	{ "thresholdRec", "TNF/chemokine threshold for Tcyt recruitment", "", "Tcyt" },
	{ "probRec", "Probability of recruiting a Tcyt cell", "", "Tcyt" },
	{ "probKillMac", "Probability of a Tcyt cell killing a (chronically) infected macrophage", "", "Tcyt" },
	{ "probKillMacCleanly", "Probability of a Tcyt cell killing a chronically infected macrophage cleanly", "", "Tcyt" },
	{ "thresholdRec", "TNF/chemokine threshold for Treg recruitment", "", "Treg" },
	{ "probRec", "Probability of recruiting a Treg cell", "", "Treg" },
	{ "growthRateIntMtb", "Growth rate of intracellular bacteria", "", "Mtb" },
	{ "growthRateExtMtb", "Growth rate of extracellular bacteria", "", "Mtb" },
	{ "growthExtMtbBound", "Upper bound on the number of extracellular bacteria used in growth function", "#bacteria", "Mtb" },
	{ "muMDC_LN", "ODE stuff", "", "GR" },
	{ "sn4", "ODE stuff", "", "GR" },
	{ "muN4", "ODE stuff", "", "GR" },
	{ "k13", "ODE stuff", "", "GR" },
	{ "hs13", "ODE stuff", "", "GR" },
	{ "k14", "ODE stuff", "", "GR" },
	{ "k15", "ODE stuff", "", "GR" },
	{ "rho2", "ODE stuff", "", "GR" },
	{ "k20a", "ODE stuff", "", "GR" },
	{ "hs20a", "ODE stuff", "", "GR" },
	{ "csi1", "ODE stuff", "", "GR" },
	{ "csi1a", "ODE stuff", "", "GR" },
	{ "sn8", "ODE stuff", "", "GR" },
	{ "muN8", "ODE stuff", "", "GR" },
	{ "wT80", "ODE stuff", "", "GR" },
	{ "k16", "ODE stuff", "", "GR" },
	{ "hs16", "ODE stuff", "", "GR" },
	{ "k17", "ODE stuff", "", "GR" },
	{ "hs17", "ODE stuff", "", "GR" },
	{ "k18", "ODE stuff", "", "GR" },
	{ "rho3", "ODE stuff", "", "GR" },
	{ "k24a", "ODE stuff", "", "GR" },
	{ "hs24a", "ODE stuff", "", "GR" },
	{ "csi2", "ODE stuff", "", "GR" },
	{ "csi2a", "ODE stuff", "", "GR" },
	{ "csi2b", "ODE stuff", "", "GR" },
	{ "scaling", "ODE stuff", "", "GR" },
	{ "m", "ODE stuff", "", "GR" },
	/* INT */
	{ "nrSources", "Number of vascular sources on the grid", "", "GR" },
	{ "nrKillingsCaseation", "Number of killings for a compartment to become caseated", "", "GR" },
	{ "maxAge", "Maximal macrophage age", "#timesteps", "Mac" },
	{ "maxAgeAct", "Maximal activated macrophage age", "#timesteps", "Mac" },
	{ "initNumber", "Initial number of resting macrophages on the grid", "", "Mac" },
	{ "maxAge", "Maximal T cell age", "#timesteps", "Tcell" },
	{ "timeRecEnabled", "Time after which T cell recruitment is enabled", "#timesteps", "Tcell" },
	{ "maxTimeReg", "Time span during which a Tgam cell remains down-regulated", "#timesteps", "Tgam" },
	{ "maxTimeReg", "Time span during which a Tcyt cell remains down-regulated", "#timesteps", "Tcyt" },
	{ "movementRest", "Time required for a resting macrophage to move one micro-compartment", "#timesteps", "Mac" },
	{ "movementInf", "Time required for an infected macrophage to move one micro-compartment", "#timesteps", "Mac" },
	{ "movementAct", "Time required for an active macrophage to move one micro-compartment", "#timesteps", "Mac" },
	{ "maxTimeReg", "Time span during which a macrophage remains down-regulated", "#timesteps", "Mac" },
};

Params::Params(bool ode)
	: _useRecruitmentWeights(false)
	, _ode(ode)
	, _doubleParam()
	, _intParam()
	, _initialMacs()
	, _initialExtMtb()
{
	memset(_doubleParam, 0, sizeof(double) * PARAM_DOUBLE_COUNT);
	memset(_intParam, 0, sizeof(int) * PARAM_INT_COUNT);
}

Params::~Params()
{
}

bool Params::reinit(const char* filename)
{
	Params* pNewInstance = new Params();
	if (!pNewInstance->fromXml(filename))
	{
		delete pNewInstance;
		return false;
	}
	else
	{
		delete _pInstance;
		_pInstance = pNewInstance;
		return true;
	}
}

bool Params::readParam(const TiXmlElement* pElement, const char* paramName, double* pVar, bool prob)
{
	switch (pElement->QueryDoubleAttribute(paramName, pVar))
	{
	case TIXML_NO_ATTRIBUTE:
		std::cerr << "Expected attribute '" << pElement->Value() << "/@"
			<< paramName << "'" << std::endl;
		return false;
	case TIXML_WRONG_TYPE:
		std::cerr << "Value of attribute '" << pElement->Value() << "/@" 
			<< paramName << "' must be a double" << std::endl;
		return false;
	}

	if (prob && !(0 <= *pVar && *pVar <= 1))
	{
		std::cerr << "Value of attribute '" << pElement->Value() << "/@" 
			<< paramName << "' must be in the range [0,1]" << std::endl;
		return false;
	}
	
	return true;
}

bool Params::readParam(const TiXmlElement* pElement, const char* paramName, double* pVar, double defaultVal, bool prob)
{
	switch (pElement->QueryDoubleAttribute(paramName, pVar))
	{
	case TIXML_NO_ATTRIBUTE:
		*pVar = defaultVal;
		return true;
	case TIXML_WRONG_TYPE:
		std::cerr << "Value of attribute '" << pElement->Value() << "/@" 
			<< paramName << "' must be a double" << std::endl;
		return false;
	}
	
	if (prob && !(0 <= *pVar && *pVar <= 1))
	{
		std::cerr << "Value of attribute '" << pElement->Value() << "/@" 
			<< paramName << "' must be in the range [0,1]" << std::endl;
		return false;
	}

	return true;
}

bool Params::readParam(const TiXmlElement* pElement, const char* paramName, int* pVar, bool pos)
{
	switch (pElement->QueryIntAttribute(paramName, pVar))
	{
	case TIXML_NO_ATTRIBUTE:
		std::cerr << "Expected attribute '" << pElement->Value() << "/@"
			<< paramName << "'" << std::endl;
		return false;
	case TIXML_WRONG_TYPE:
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a" << (pos ? "n " : " positive ")
			<< "integer" << std::endl;
		return false;
	}

	if (pos && *pVar < 0)
	{
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a positive integer" << std::endl;
		return false;
	}

	return true;
}

bool Params::readParam(const TiXmlElement* pElement, const char* paramName, int* pVar, int defaultVal, bool pos)
{
	switch (pElement->QueryIntAttribute(paramName, pVar))
	{
	case TIXML_NO_ATTRIBUTE:
		*pVar = defaultVal;
		return true;
	case TIXML_WRONG_TYPE:
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a" << (pos ? "n " : " positive ")
			<< "integer" << std::endl;
		return false;
	}

	if (pos && *pVar < 0)
	{
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a positive integer" << std::endl;
		return false;
	}

	return true;
}

bool Params::readParam(const TiXmlElement* pElement, ParamDoubleType param, bool prob)
{
	const char* paramName = getName(param);
	double* pVar = _doubleParam + param;

	return readParam(pElement, paramName, pVar, prob);
}

bool Params::readParam(const TiXmlElement* pElement, ParamDoubleType param,
	double defaultVal, bool prob)
{
	const char* paramName = getName(param);
	double* pVar = _doubleParam + param;

	return readParam(pElement, paramName, pVar, defaultVal, prob);
}

bool Params::readParam(const TiXmlElement* pElement, ParamIntType param, bool pos)
{
	const char* paramName = getName(param);
	int* pVar = _intParam + param;

	return readParam(pElement, paramName, pVar, pos);
}

bool Params::readParam(const TiXmlElement* pElement, ParamIntType param, int defaultVal, bool pos)
{
	const char* paramName = getName(param);
	int* pVar = _intParam + param;

	return readParam(pElement, paramName, pVar, defaultVal, pos);
}

void Params::writeParam(std::ostream& out, ParamDoubleType param,
	int indentation, ClosingTagType close) const
{
	for (int i = 0; i < indentation; i++)
		out << "\t";

	out << getName(param) << " = \"" << getParam(param);

	switch (close)
	{
	case CLOSE_START_TAG:
		out << "\">\n";
		break;
	case CLOSE_END_TAG:
		out << "\"/>\n";
		break;
	case CLOSE_NONE:
		out << "\"\n";
	}
}

void Params::writeParam(std::ostream& out, ParamIntType param,
	int indentation, ClosingTagType close) const
{
	for (int i = 0; i < indentation; i++)
		out << "\t";

	out << getName(param) << " = \"" << getParam(param);

	switch (close)
	{
	case CLOSE_START_TAG:
		out << "\">\n";
		break;
	case CLOSE_END_TAG:
		out << "\"/>\n";
		break;
	case CLOSE_NONE:
		out << "\"\n";
	}
}

bool Params::toXml(const char* filename) const
{
	std::ofstream outFile(filename);
	if (!outFile.good())
	{
		return false;
	}

	outFile << "<GR\n";
	writeParam(outFile, PARAM_GR_NR_KILLINGS_FOR_CASEATION, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_NR_SOURCES, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_D_TNF, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_D_SHED_TNFR2, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_D_CHEMOKINES, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_DEG_TNF, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_DEG_CHEMOKINES, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_THRESHOLD_APOPTOSIS_TNF, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_APOPTOSIS, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_PROB_APOPTOSIS_TNF, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MIN_CHEMOTAXIS, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_SEC_RATE_ATTRACTANT, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MAX_CHEMOTAXIS, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_SYNTH_MAC, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_SYNTH_TCELL, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_TACE_MAC, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_TACE_TCELL, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_KD1, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_KD2, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_ON1, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_ON2, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_INT1, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_INT2, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_SHED, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_REC1, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_REC2, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_T1, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_T2, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_DEG1, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_K_DEG2, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MAX_TNFR1_MAC, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MIN_TNFR1_MAC, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MAX_TNFR2_MAC, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MIN_TNFR2_MAC, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MAX_TNFR1_TCELL, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MIN_TNFR1_TCELL, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MAX_TNFR2_TCELL, 1, CLOSE_NONE);
	writeParam(outFile, PARAM_GR_MIN_TNFR2_TCELL, 1,
		_useRecruitmentWeights ? CLOSE_NONE : CLOSE_START_TAG);

	if (_useRecruitmentWeights)
	{
		writeParam(outFile, PARAM_GR_WEIGHT_TNF_RECRUITMENT, 1, CLOSE_NONE);
		writeParam(outFile, PARAM_GR_WEIGHT_CCL2_RECRUITMENT, 1, CLOSE_NONE);
		writeParam(outFile, PARAM_GR_WEIGHT_CCL5_RECRUITMENT, 1, CLOSE_NONE);
		writeParam(outFile, PARAM_GR_WEIGHT_CXCL9_RECRUITMENT, 1, CLOSE_START_TAG);
	}

	outFile << "\t<Mac\n";
	writeParam(outFile, PARAM_MAC_INIT_NUMBER, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_AGE, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_A_AGE, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_TIMESPAN_REGULATED, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_MOVEMENT_RESTING, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_MOVEMENT_ACTIVE, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_MOVEMENT_INFECTED, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_SEC_RATE_TNF, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_SEC_RATE_CCL2, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_SEC_RATE_CCL5, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_SEC_RATE_CXCL9, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_THRESHOLD_NFKB_TNF, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_K_NFKB, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_NR_UPTAKE_RI_EXTMTB, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_PROB_KILL_R_EXTMTB, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_THRESHOLD_NFKB_EXTMTB, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_THRESHOLD_BURST_CI_INTMTB, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_PROB_STAT1_TGAM, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_NR_UPTAKE_A_EXTMTB, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_THRESHOLD_RECRUITMENT, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_MAC_PROB_RECRUITMENT, 2, CLOSE_END_TAG);

	outFile << "\t<Tcell\n";
	writeParam(outFile, PARAM_TCELL_AGE, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_TCELL_TIME_RECRUITMENT_ENABLED, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_TCELL_PROB_MOVE_TO_MAC, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_TCELL_PROB_MOVE_TO_TCELL, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_TCELL_PROB_RECRUITMENT, 2, CLOSE_START_TAG);

	outFile << "\t\t<Tgam\n";
	writeParam(outFile, PARAM_TGAM_TIMESPAN_REGULATED, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TGAM_PROB_APOPTOSIS_FAS_FASL, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TGAM_THRESHOLD_RECRUITMENT, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TGAM_PROB_RECRUITMENT, 3, CLOSE_END_TAG);

	outFile << "\t\t<Tcyt\n";
	writeParam(outFile, PARAM_TCYT_TIMESPAN_REGULATED, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TCYT_PROB_KILL_MAC, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TCYT_PROB_KILL_MAC_CLEANLY, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TCYT_THRESHOLD_RECRUITMENT, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TCYT_PROB_RECRUITMENT, 3, CLOSE_END_TAG);

	outFile << "\t\t<Treg\n";
	writeParam(outFile, PARAM_TREG_THRESHOLD_RECRUITMENT, 3, CLOSE_NONE);
	writeParam(outFile, PARAM_TREG_PROB_RECRUITMENT, 3, CLOSE_END_TAG);

	outFile << "\t</Tcell>\n";

	outFile << "\t<Mtb\n";
	writeParam(outFile, PARAM_INTMTB_GROWTH_RATE, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_EXTMTB_GROWTH_RATE, 2, CLOSE_NONE);
	writeParam(outFile, PARAM_EXTMTB_UPPER_BOUND, 2, CLOSE_END_TAG);

	outFile << "\t<Init>\n";
	for (PosVector::const_iterator it = _initialMacs.begin(); it != _initialMacs.end(); it++)
	{
		outFile << "\t\t<Mac row = \"" << it->first << "\" col = \"" << it->second << "\"/>\n";
	}
	for (PosVector::const_iterator it = _initialExtMtb.begin(); it != _initialExtMtb.end(); it++)
	{
		outFile << "\t\t<ExtMtb row = \"" << it->first << "\" col = \"" << it->second << "\"/>\n";
	}
	outFile << "\t</Init>\n" << "</GR>\n";

	outFile.flush();
	outFile.close();

	return true;
}

bool Params::readGRElement(const TiXmlElement* pGrElement)
{
	assert(pGrElement);

	const TiXmlElement* pMtbElement = pGrElement->FirstChildElement("Mtb");
	if (!pMtbElement)
	{
		std::cerr << "Expected element '/GR/Mtb'" << std::endl;
		return false;
	}

	const TiXmlElement* pMacElement = pGrElement->FirstChildElement("Mac");
	if (!pMacElement)
	{
		std::cerr << "Expected element '/GR/Mac'" << std::endl;
		return false;
	}

	const TiXmlElement* pTcellElement = pGrElement->FirstChildElement("Tcell");
	if (!pTcellElement)
	{
		std::cerr << "Expected element '/GR/Tcell'" << std::endl;
		return false;
	}

	bool res = true;

	res &= readParam(pGrElement, PARAM_GR_NR_KILLINGS_FOR_CASEATION, true);
	res &= readParam(pGrElement, PARAM_GR_NR_SOURCES, true);
	res &= readParam(pGrElement, PARAM_GR_D_TNF, false);
	res &= readParam(pGrElement, PARAM_GR_D_SHED_TNFR2, false);
	res &= readParam(pGrElement, PARAM_GR_D_CHEMOKINES, false);
	res &= readParam(pGrElement, PARAM_GR_DEG_TNF, false);
	res &= readParam(pGrElement, PARAM_GR_DEG_CHEMOKINES, false);
	res &= readParam(pGrElement, PARAM_GR_THRESHOLD_APOPTOSIS_TNF, false);
	res &= readParam(pGrElement, PARAM_GR_K_APOPTOSIS, false);
	res &= readParam(pGrElement, PARAM_GR_PROB_APOPTOSIS_TNF, true);
	res &= readParam(pGrElement, PARAM_GR_MIN_CHEMOTAXIS, false);
	res &= readParam(pGrElement, PARAM_GR_MAX_CHEMOTAXIS, false);
	res &= readParam(pGrElement, PARAM_GR_K_SYNTH_MAC, false);
	res &= readParam(pGrElement, PARAM_GR_K_SYNTH_TCELL, false);
	res &= readParam(pGrElement, PARAM_GR_K_TACE_MAC, false);
	res &= readParam(pGrElement, PARAM_GR_K_TACE_TCELL, false);
	res &= readParam(pGrElement, PARAM_GR_KD1, false);
	res &= readParam(pGrElement, PARAM_GR_KD2, false);
	res &= readParam(pGrElement, PARAM_GR_K_ON1, false);
	res &= readParam(pGrElement, PARAM_GR_K_ON2, false);
	res &= readParam(pGrElement, PARAM_GR_K_INT1, false);
	res &= readParam(pGrElement, PARAM_GR_K_INT2, false);
	res &= readParam(pGrElement, PARAM_GR_K_SHED, false);
	res &= readParam(pGrElement, PARAM_GR_K_REC1, false);
	res &= readParam(pGrElement, PARAM_GR_K_REC2, false);
	res &= readParam(pGrElement, PARAM_GR_K_T1, false);
	res &= readParam(pGrElement, PARAM_GR_K_T2, false);
	res &= readParam(pGrElement, PARAM_GR_K_DEG1, false);
	res &= readParam(pGrElement, PARAM_GR_K_DEG2, false);
	res &= readParam(pGrElement, PARAM_GR_MAX_TNFR1_MAC, false);
	res &= readParam(pGrElement, PARAM_GR_MIN_TNFR1_MAC, false);
	res &= readParam(pGrElement, PARAM_GR_MAX_TNFR2_MAC, false);
	res &= readParam(pGrElement, PARAM_GR_MIN_TNFR2_MAC, false);
	res &= readParam(pGrElement, PARAM_GR_MAX_TNFR1_TCELL, false);
	res &= readParam(pGrElement, PARAM_GR_MIN_TNFR1_TCELL, false);
	res &= readParam(pGrElement, PARAM_GR_MAX_TNFR2_TCELL, false);
	res &= readParam(pGrElement, PARAM_GR_MIN_TNFR2_TCELL, false);
	res &= readParam(pGrElement, PARAM_GR_SEC_RATE_ATTRACTANT, false);

	res &= readMacElement(pMacElement);
	res &= readTcellElement(pTcellElement);
	res &= readMtbElement(pMtbElement);
	res &= readInitElement(pGrElement);

	// the following four parameters are optional
	_useRecruitmentWeights = pGrElement->Attribute(getName(PARAM_GR_WEIGHT_TNF_RECRUITMENT)) ||
			pGrElement->Attribute(getName(PARAM_GR_WEIGHT_CCL2_RECRUITMENT)) ||
			pGrElement->Attribute(getName(PARAM_GR_WEIGHT_CCL5_RECRUITMENT)) ||
			pGrElement->Attribute(getName(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT));

	res &= readParam(pGrElement, PARAM_GR_WEIGHT_TNF_RECRUITMENT, 1, false);
	res &= readParam(pGrElement, PARAM_GR_WEIGHT_CCL2_RECRUITMENT,
		getParam(PARAM_MAC_SEC_RATE_TNF) / getParam(PARAM_MAC_SEC_RATE_CCL2), false);
	res &= readParam(pGrElement, PARAM_GR_WEIGHT_CCL5_RECRUITMENT,
		getParam(PARAM_MAC_SEC_RATE_TNF) / getParam(PARAM_MAC_SEC_RATE_CCL5), false);
	res &= readParam(pGrElement, PARAM_GR_WEIGHT_CXCL9_RECRUITMENT,
		getParam(PARAM_MAC_SEC_RATE_TNF) / getParam(PARAM_MAC_SEC_RATE_CXCL9), false);

	if (_ode)
	{
		res &= readParam(pGrElement, PARAM_muMDC_LN, false);
		res &= readParam(pGrElement, PARAM_sn4, false);
		res &= readParam(pGrElement, PARAM_muN4, false);
		res &= readParam(pGrElement, PARAM_k13, false);
		res &= readParam(pGrElement, PARAM_hs13, false);
		res &= readParam(pGrElement, PARAM_k14, false);
		res &= readParam(pGrElement, PARAM_k15, false);
		res &= readParam(pGrElement, PARAM_rho2, false);
		res &= readParam(pGrElement, PARAM_k20a, false);
		res &= readParam(pGrElement, PARAM_hs20a, false);
		res &= readParam(pGrElement, PARAM_csi1, false);
		res &= readParam(pGrElement, PARAM_csi1a, false);
		res &= readParam(pGrElement, PARAM_sn8, false);
		res &= readParam(pGrElement, PARAM_muN8, false);
		res &= readParam(pGrElement, PARAM_wT80, false);
		res &= readParam(pGrElement, PARAM_k16, false);
		res &= readParam(pGrElement, PARAM_hs16, false);
		res &= readParam(pGrElement, PARAM_k17, false);
		res &= readParam(pGrElement, PARAM_hs17, false);
		res &= readParam(pGrElement, PARAM_k18, false);
		res &= readParam(pGrElement, PARAM_rho3, false);
		res &= readParam(pGrElement, PARAM_k24a, false);
		res &= readParam(pGrElement, PARAM_hs24a, false);
		res &= readParam(pGrElement, PARAM_csi2, false);
		res &= readParam(pGrElement, PARAM_csi2a, false);
		res &= readParam(pGrElement, PARAM_csi2b, false);
		res &= readParam(pGrElement, PARAM_scaling, false);
		res &= readParam(pGrElement, PARAM_m, false);
	}

	return res;
}

bool Params::readMacElement(const TiXmlElement* pMacElement)
{
	assert(pMacElement);

	bool res = true;

	res &= readParam(pMacElement, PARAM_MAC_AGE, true);
	res &= readParam(pMacElement, PARAM_MAC_INIT_NUMBER, true);
	res &= readParam(pMacElement, PARAM_MAC_SEC_RATE_TNF, false);
	res &= readParam(pMacElement, PARAM_MAC_SEC_RATE_CCL2, false);
	res &= readParam(pMacElement, PARAM_MAC_SEC_RATE_CCL5, false);
	res &= readParam(pMacElement, PARAM_MAC_SEC_RATE_CXCL9, false);
	res &= readParam(pMacElement, PARAM_MAC_THRESHOLD_NFKB_TNF, false);
	res &= readParam(pMacElement, PARAM_MAC_K_NFKB, false);
	res &= readParam(pMacElement, PARAM_MAC_NR_UPTAKE_RI_EXTMTB, false);
	res &= readParam(pMacElement, PARAM_MAC_PROB_KILL_R_EXTMTB, false);
	res &= readParam(pMacElement, PARAM_MAC_THRESHOLD_NFKB_EXTMTB, false);
	res &= readParam(pMacElement, PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB, false);
	res &= readParam(pMacElement, PARAM_MAC_PROB_STAT1_TGAM, true);
	res &= readParam(pMacElement, PARAM_MAC_THRESHOLD_BURST_CI_INTMTB, false);
	res &= readParam(pMacElement, PARAM_MAC_NR_UPTAKE_A_EXTMTB, false);
	res &= readParam(pMacElement, PARAM_MAC_A_AGE, true);
	res &= readParam(pMacElement, PARAM_MAC_THRESHOLD_RECRUITMENT, false);
	res &= readParam(pMacElement, PARAM_MAC_PROB_RECRUITMENT, true);
	res &= readParam(pMacElement, PARAM_MAC_MOVEMENT_RESTING, true);
	res &= readParam(pMacElement, PARAM_MAC_MOVEMENT_INFECTED, true);
	res &= readParam(pMacElement, PARAM_MAC_MOVEMENT_ACTIVE, true);
	res &= readParam(pMacElement, PARAM_MAC_TIMESPAN_REGULATED, true);

	return res;
}

bool Params::readTcellElement(const TiXmlElement* pTcellElement)
{
	assert(pTcellElement);

	const TiXmlElement* pTgamElement = pTcellElement->FirstChildElement("Tgam");
	if (!pTgamElement)
	{
		std::cerr << "Expected element '/GR/Tcell/Tgam'" << std::endl;
		return false;
	}

	const TiXmlElement* pTcytElement = pTcellElement->FirstChildElement("Tcyt");
	if (!pTcytElement)
	{
		std::cerr << "Expected element '/GR/Tcell/Tcyt'" << std::endl;
		return false;
	}

	const TiXmlElement* pTregElement = pTcellElement->FirstChildElement("Treg");
	if (!pTregElement)
	{
		std::cerr << "Expected element '/GR/Tcell/Treg'" << std::endl;
		return false;
	}

	bool res = true;

	res &= readParam(pTcellElement, PARAM_TCELL_AGE, true);
	res &= readParam(pTcellElement, PARAM_TCELL_TIME_RECRUITMENT_ENABLED, true);
	res &= readParam(pTcellElement, PARAM_TCELL_PROB_MOVE_TO_MAC, true);
	res &= readParam(pTcellElement, PARAM_TCELL_PROB_MOVE_TO_TCELL, true);
	res &= readParam(pTcellElement, PARAM_TCELL_PROB_RECRUITMENT, true);

	res &= readTgamElement(pTgamElement);
	res &= readTcytElement(pTcytElement);
	res &= readTregElement(pTregElement);

	return res;
}

bool Params::readTgamElement(const TiXmlElement* pTgamElement)
{
	assert(pTgamElement);

	bool res = true;

	res &= readParam(pTgamElement, PARAM_TGAM_PROB_APOPTOSIS_FAS_FASL, true);
	res &= readParam(pTgamElement, PARAM_TGAM_PROB_RECRUITMENT, true);
	res &= readParam(pTgamElement, PARAM_TGAM_THRESHOLD_RECRUITMENT, false);
	res &= readParam(pTgamElement, PARAM_TGAM_TIMESPAN_REGULATED, true);

	return res;
}

bool Params::readTcytElement(const TiXmlElement* pTcytElement)
{
	assert(pTcytElement);

	bool res = true;

	res &= readParam(pTcytElement, PARAM_TCYT_PROB_KILL_MAC, true);
	res &= readParam(pTcytElement, PARAM_TCYT_PROB_KILL_MAC_CLEANLY, true);
	res &= readParam(pTcytElement, PARAM_TCYT_PROB_RECRUITMENT, true);
	res &= readParam(pTcytElement, PARAM_TCYT_THRESHOLD_RECRUITMENT, false);
	res &= readParam(pTcytElement, PARAM_TCYT_TIMESPAN_REGULATED, true);

	return res;
}

bool Params::readTregElement(const TiXmlElement* pTregElement)
{
	assert(pTregElement);

	bool res = true;

	res &= readParam(pTregElement, PARAM_TREG_PROB_RECRUITMENT, true);
	res &= readParam(pTregElement, PARAM_TREG_THRESHOLD_RECRUITMENT, false);

	return res;
}

bool Params::readMtbElement(const TiXmlElement* pMtbElement)
{
	assert(pMtbElement);

	bool res = true;

	res &= readParam(pMtbElement, PARAM_INTMTB_GROWTH_RATE, false);
	res &= readParam(pMtbElement, PARAM_EXTMTB_GROWTH_RATE, false);
	res &= readParam(pMtbElement, PARAM_EXTMTB_UPPER_BOUND, false);

	return res;
}

bool Params::readInitElement(const TiXmlElement* pGrElement)
{
	assert(pGrElement);

	_initialMacs.clear();
	_initialExtMtb.clear();

	bool res = true;

	const TiXmlElement* pInitElement = pGrElement->FirstChildElement("Init");
	if (pInitElement)
	{
		for (const TiXmlElement* pInitChildElement = pInitElement->FirstChildElement();
			pInitChildElement;
			pInitChildElement = pInitChildElement->NextSiblingElement())
		{
			Pos pos;
			if (!strcmp(pInitChildElement->Value(), "Mac"))
			{
				res &= readParam(pInitChildElement, "row", &pos.first, true);
				res &= readParam(pInitChildElement, "col", &pos.second, true);

				if (!(0 <= pos.first && pos.first < NROWS) ||
					!(0 <= pos.second && pos.second < NCOLS))
				{
					res = false;
					std::cerr << "Position (" << pos.first << "," << pos.second
							<< ") is out of range" << std::endl;
				}
				else
				{
					_initialMacs.push_back(pos);
				}
			}
			else if (!strcmp(pInitChildElement->Value(), "ExtMtb"))
			{
				res &= readParam(pInitChildElement, "row", &pos.first, true);
				res &= readParam(pInitChildElement, "col", &pos.second, true);

				if (!(0 <= pos.first && pos.first < NROWS) ||
					!(0 <= pos.second && pos.second < NCOLS))
				{
					res = false;
					std::cerr << "Position (" << pos.first << "," << pos.second
							<< ") is out of range" << std::endl;
				}
				else
				{
					_initialExtMtb.push_back(pos);
				}
			}
			else
			{
				std::cerr << "Unexpected element with name '"
						<< pInitChildElement->Value() << "'" << std::endl;
				res = false;
			}
		}
	}
	else
	{
		Pos pos(NROWS/2, NCOLS/2);
		_initialMacs.push_back(pos);
	}

	return res;
}

bool Params::fromXml(const char* filename)
{
	TiXmlDocument xmlDoc(filename);

	if (!xmlDoc.LoadFile())
	{
		std::cerr << xmlDoc.ErrorDesc() << " '" << filename << "'" << std::endl;
		return false;
	}

	TiXmlElement* pGrElement = xmlDoc.RootElement();
	if (!pGrElement || strcmp(pGrElement->Value(), "GR"))
	{
		std::cerr << "Expected root element with name 'GR'" << std::endl;
		return false;
	}

	return readGRElement(pGrElement);
}

void Params::updateRecruitmentWeights()
{
	setParam(PARAM_GR_WEIGHT_CCL2_RECRUITMENT,
		getParam(PARAM_MAC_SEC_RATE_TNF) / getParam(PARAM_MAC_SEC_RATE_CCL2));
	setParam(PARAM_GR_WEIGHT_CCL5_RECRUITMENT,
		getParam(PARAM_MAC_SEC_RATE_TNF) / getParam(PARAM_MAC_SEC_RATE_CCL5));
	setParam(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT,
		getParam(PARAM_MAC_SEC_RATE_TNF) / getParam(PARAM_MAC_SEC_RATE_CXCL9));
}

void Params::printCSV(bool header) const
{
	if (header)
	{
		int n = (int) PARAM_DOUBLE_COUNT;
		for (int i = 0; i < n; i++)
		{
			if (i != 0)
				std::cout << ',';

			std::cout << '"' << getType((ParamDoubleType) i) << '@' << getName((ParamDoubleType) i) << '"';
		}

		n = (int) PARAM_INT_COUNT;
		for (int i = 0; i < n; i++)
		{
			if (i != 0)
				std::cout << ',';

			std::cout << '"' << getType((ParamIntType) i) << '@' << getName((ParamIntType) i) << '"';
		}

		std::cout << std::endl;
	}

	int n = (int) PARAM_DOUBLE_COUNT;
	for (int i = 0; i < n; i++)
	{
		if (i != 0)
			std::cout << ',';

		std::cout << '"' << getParam((ParamDoubleType) i) << '"';
	}

	n = (int) PARAM_INT_COUNT;
	for (int i = 0; i < n; i++)
	{
		if (i != 0)
			std::cout << ',';

		std::cout << '"' << getParam((ParamIntType) i) << '"';
	}

	std::cout << std::endl;
}
