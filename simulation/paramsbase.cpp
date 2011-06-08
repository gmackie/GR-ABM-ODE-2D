/*
 * params.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "paramsbase.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

// For each node, its: node type, parent, text name, whether or not it is required.
// The first element of this array is considered the root node and must be present in the parameter file as the root node.
const NodeDescription ParamsBase::_element[NODE_COUNT] = {
   // type			parent		name		required
	{ GR_NODE, 		NULL_NODE, 	"GR" },
	{ MAC_NODE, 	GR_NODE, 	"Mac" },
	{ TCELL_NODE, 	GR_NODE, 	"Tcell" },
	{ TGAM_NODE, 	TCELL_NODE, "Tgam" },
	{ TCYT_NODE, 	TCELL_NODE, "Tcyt" },
	{ TREG_NODE, 	TCELL_NODE, "Treg" },
	{ MTB_NODE, 	GR_NODE,	"Mtb" },
	{ INIT_NODE, 	GR_NODE,	"Init", }
};

const int ParamsBase::_PARAM_COUNT = PARAM_DOUBLE_COUNT + PARAM_INT_COUNT;

// The _description array defines the attributes for the parameters.
// The order of the items in this array must exactly match the order of the items in the ParamDoubleType and ParamIntType enums,
// with all the double parameters appearing first in _description followed by the int parameters, since the ParamDoubleType
// and ParamIntType items are used as indices into _description.

// To add a parameter, make a new entry in this array and also add a parameter enum item to the corresponding
// position of the ParamDoubleType enum, if the parameter is to have floating point values, or to the
// corresponding position of the ParamIntType enum, if the parameter is to have integer values

// Give it a name, which does not have to be unique, it only has to be unique within the particular XML element it belongs to.
// For example, there can be only one MAC_NODE probRec, but there can be a MAC_NODE probRec and a TCYT_NODE probRec.

// The XmlElement is one of the items in the XmlElement enum defined in paramsbase.h.
// It specifies where in the parameter XML file this parameter will appear as an attribute.
// If this parameter is for a new XML element that wasn't used before, then also add a new entry in the
// XmlElement enum and to the _element array.

// If probPos is true it specifies a constraint on the parameter value.
// For a double parameter it means it is a probability, and so must be in the range [0.0, 1.0], inclusive.
// For an int parameter is means it must be positive (really non-negative), i.e. > 0.

// useDefault specifies whether or not a default value should be used if the parameter does not appear in a parameter
// file being read.
// If false and the parameter does not appear in a parameter file that is treated as an error.
// If true and the parameter does not appear in a parameter file either the double or int default value is used,
// based on whether or not the parameter is a double parameter or an int parameter.
//
// useDefault should be set true if the parameter is to be calculated from another parameter, so it isn't flagged
// as a missing parameter when it is not present in a parameter file.
//
// double Default is the double value to use as a default.
// int Default is the int value to use as a default.
// Both need to be specified regardless of the value for useDefault.

const ParamDescription ParamsBase::_description[_PARAM_COUNT] =
{
	// 														     use	  default
	// 										XmlElement	probPos	Default	double	int	  unit				description
	{ "diffusivityTNF",						GR_NODE,	false,	false,	0.0,	0,	"cm^2/s",			"TNF diffusivity" },
	{ "diffusivityShedTNFR2",				GR_NODE,	false,	false,	0.0,	0,	"cm^2/s",			"Shed sTNF/TNFR2 comples diffusivity" },
	{ "diffusivityChemokines",				GR_NODE,	false,	false,	0.0,	0,	"cm^2/s",			"Chemokine diffusivity" },
	{ "degRateTNF",							GR_NODE,	false,	false,	0.0,	0,	"",					"TNF degradation rate (per 6s)" },
	{ "degRateChemokines",					GR_NODE,	false,	false,	0.0,	0,	"",					"Chemokine degradation rate (per 6s)" },
	{ "minChemotaxis",						GR_NODE,	false,	false,	0.0,	0,	"#molecules",		"Chemotaxis sensivity range (lower bound)" },
	{ "maxChemotaxis",						GR_NODE,	false,	false,	0.0,	0,	"#molecules",		"Chemotaxis sensivity range (upper bound)" },
	{ "dTNF_Tgam",							GR_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of TNF by Tgam" },
	{ "dTNF_Tcyt",							GR_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of TNF by Tgam" },
	// molecular TNF-associated parameters
	{ "kSynthMac",							GR_NODE,	false,	false,	0.0,	0,	"#/cell.sec",		"Rate of mTNF synthesis by a macrophage" },
	{ "kSynthTcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell.sec",		"Rate of mTNF synthesis by a T cell" },
	{ "kTaceMac",							GR_NODE,	false,	false,	0.0,	0,	"1/sec",			"Rate of mTNF release from a macrophage" },
	{ "kTaceTcell",							GR_NODE,	false,	false,	0.0,	0,	"1/sec",			"Rate of mTNF release from a T cell" },
	{ "KD1",								GR_NODE,	false,	false,	0.0,	0,	"M",				"sTNF/TNFR1 dissociation equilibrium constant" },
	{ "KD2",								GR_NODE,	false,	false,	0.0,	0,	"M",				"sTNF/TNFR2 dissociation equilibrium constant" },
	{ "kOn1",								GR_NODE,	false,	false,	0.0,	0,	"1/M.s",			"sTNF/TNFR1 association rate constant" },
	{ "kOn2",								GR_NODE,	false,	false,	0.0,	0,	"1/M.s",			"sTNF/TNFR2 association rate constant" },
	{ "kInt1",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"sTNF/TNFR1 complex internaliation rate constant" },
	{ "kInt2",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"sTNF/TNFR2 complex internaliation rate constant" },
	{ "kShed",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"sTNF/TNFR2 complex shedding rate constant" },
	{ "kRec1",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"TNFR1 recycling rate constant" },
	{ "kRec2",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"TNFR2 recycling rate constant" },
	{ "kT1",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"TNFR1 turn-over rate constant" },
	{ "kT2",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"TNFR2 turn-over rate constant" },
	{ "kDeg1",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"TNFR1 degradtion rate constant" },
	{ "kDeg2",								GR_NODE,	false,	false,	0.0,	0,	"1/s",				"TNFR2 degradtion rate constant" },
	{ "maxTNFR1Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Maximum density of TNFR1 on macrophages" },
	{ "minTNFR1Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Minimum density of TNFR1 on macrophages" },
	{ "maxTNFR2Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Maximum density of TNFR2 on macrophages" },
	{ "minTNFR2Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Minimum density of TNFR2 on macrophages" },
	{ "maxTNFR1Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Maximum density of TNFR1 on T cells" },
	{ "minTNFR1Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Minimum density of TNFR1 on T cells" },
	{ "maxTNFR2Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Maximum density of TNFR2 on T cells" },
	{ "minTNFR2Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Minimum density of TNFR2 on T cells" },
	// end of molecular TNF-associated parameters
	{ "thresholdApoptosisTNF",				GR_NODE,	false,	false,	0.0,	0,	"fraction",			"TNF threshold for TNF-induced apoptosis" },
	{ "kApoptosis",							GR_NODE,	false,	false,	0.0,	0,	"1/s",				"Rate of apoptosis happening" },
	{ "thresholdApoptosisTNF_Molecular",	GR_NODE,	false,	false,	0.0,	0,	"#molecules",		"TNF threshold for TNF-induced apoptosis" },
	{ "kApoptosis_Molecular",				GR_NODE,	false,	false,	0.0,	0,	"1/s",				"Rate of apoptosis happening" },
	{ "probApoptosisTNF",					GR_NODE,	true,	false,	0.0,	0,	"",					"Probability of TNF-induced apoptosis" },
	{ "effectRecTNF",						GR_NODE,	false,	true,	1.0,	1,	"",					"Effect of TNF on recruitment" },
	{ "effectRecCCL2",						GR_NODE,	false,	true,	0.0,	0,	"",					"Effect of CCL2 on recruitment" },
	{ "effectRecCCL5",						GR_NODE,	false,	true,	0.0,	0,	"",					"Effect of CCL5 on recruitment" },
	{ "effectRecCXCL9",						GR_NODE,	false,	true,	0.0,	0,	"",					"Effect of CXCL9 on recruitment" },
	{ "dAttractant",						GR_NODE,	false,	false,	0.0,	0,	"",					"Secretion rate of macrophage attractant secreted by debris" },
	{ "dCCL2",								MAC_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of CCL2" },
	{ "dCCL5",								MAC_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of CCL5" },
	{ "dCXCL9",								MAC_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of CXCL9" },
	{ "dTNF",								MAC_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of TNF" },
	{ "thresholdNFkBTNF",					MAC_NODE,	false,	false,	0.0,	0,	"fraction",			"TNF threshold for NFkB activation" },
	{ "kNFkB",								MAC_NODE,	false,	false,	0.0,	0,	"1/s",				"Rate of NFkB activation" },
	{ "thresholdNFkBTNF_Molecular",			MAC_NODE,	false,	false,	0.0,	0,	"#molecules",		"TNF threshold for NFkB activation" },
	{ "kNFkB_Molecular",					MAC_NODE,	false,	false,	0.0,	0,	"1/s",				"Rate of NFkB activation" },
	{ "nrExtMtbUptakeRest",					MAC_NODE,	false,	false,	0.0,	0,	"#bacteria",		"Number of extracellular bacteria a resting macrophage can uptake/kill" },
	{ "probKillExtMtbRest",					MAC_NODE,	false,	false,	0.0,	0,	"",					"Probability of a resting macrophage to kill some extracellular bacteria" },
	{ "nrExtMtbNFkB",						MAC_NODE,	false,	false,	0.0,	0,	"#bacteria",		"Number of extracellular bacteria for NFkB activation in an infected macrophage" },
	{ "nrIntMtbCInf",						MAC_NODE,	false,	false,	0.0,	0,	"#bacteria",		"Number of intracellular bacteria for an infected macrophage to become chronically infected" },
	{ "nrIntMtbBurstCInf",					MAC_NODE,	false,	false,	0.0,	0,	"#bacteria",		"Number of intracellular bacteria for an chronically infected macrophage to burst" },
	{ "probStat1Tgam",						MAC_NODE,	true,	false,	0.0,	0,	"",					"Probability of STAT1 activation due to a Tgam cell" },
	{ "nrExtMtbUptakeAct",					MAC_NODE,	false,	false,	0.0,	0,	"#bacteria",		"Number of extracellular bacteria an active macrophage can uptake and subsequently kill" },
	{ "thresholdRec",						MAC_NODE,	false,	false,	0.0,	0,	"",					"TNF/chemokine threshold for macrophage recruitment" },
	{ "probRec",							MAC_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a macrophage" },
	{ "movementBonusFactor",				MAC_NODE,	false,	true,	1.0,	1,	"",					"Bonus factor for neighbor compartment with highest chemical concentration" },
	{ "probRec",							TCELL_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a T cell" },
	{ "probMoveToMac",						TCELL_NODE,	true,	false,	0.0,	0,	"",					"Probability of a T cell moving onto a compartment already containing a macrophage" },
	{ "probMoveToTcell",					TCELL_NODE,	true,	false,	0.0,	0,	"",					"Probability of a T cell moving onto a compartment already containing another T cell" },
	{ "movementBonusFactor",				TCELL_NODE,	false,	true,	1.0,	1,	"",					"Bonus factor for neighbor compartment with highest chemical concentration" },
	{ "probApoptosisFasFasL",				TGAM_NODE,	true,	false,	0.0,	0,	"",					"Probability of Fas/FasL induced apoptosis by a Tgam cell" },
	{ "thresholdRec",						TGAM_NODE,	false,	false,	0.0,	0,	"",					"TNF/chemokine threshold for Tgam recruitment" },
	{ "probRec",							TGAM_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a Tgam cell" },
	{ "thresholdRec",						TCYT_NODE,	false,	false,	0.0,	0,	"",					"TNF/chemokine threshold for Tcyt recruitment" },
	{ "probRec",							TCYT_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a Tcyt cell" },
	{ "probKillMac",						TCYT_NODE,	true,	false,	0.0,	0,	"",					"Probability of a Tcyt cell killing a (chronically) infected macrophage" },
	{ "probKillMacCleanly",					TCYT_NODE,	true,	false,	0.0,	0,	"",					"Probability of a Tcyt cell killing a chronically infected macrophage cleanly" },
	{ "thresholdRec",						TREG_NODE,	false,	false,	0.0,	0,	"",					"TNF/chemokine threshold for Treg recruitment" },
	{ "probRec",							TREG_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a Treg cell" },
	{ "growthRateIntMtb",					MTB_NODE,	false,	false,	0.0,	0,	"",					"Growth rate of intracellular bacteria" },
	{ "growthRateExtMtb",					MTB_NODE,	false,	false,	0.0,	0,	"",					"Growth rate of extracellular bacteria" },
	{ "growthExtMtbBound",					MTB_NODE,	false,	false,	0.0,	0,	"#bacteria",		"Upper bound on the number of extracellular bacteria used in growth function" },
	{ "muMDC_LN",							GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "muN4",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k13",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "hs13",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k14",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k15",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "rho2",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k20a",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "hs20a",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "csi1",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "csi1a",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "muN8",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "wT80",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k16",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "hs16",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k17",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "hs17",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k18",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "rho3",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "k24a",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "hs24a",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "csi2",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "csi2a",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "csi2b",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "scaling",							GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "m",									GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff" },
	{ "scalingMDC",							GR_NODE,	false,	true,	1.0,	1,	"",					"ODE stuff" },
	{ "scalingLung",						GR_NODE,	false,	true,	1.0,	1,	"",					"ODE stuff - scaling factor for MDC proxy to account for GR size" },
	{ "scalingLN",							GR_NODE,	false,	true,	1.0,	1,	"",					"ODE stuff - scaling factor for T cell fluxes to account for GR size" },
	{ "initN4",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff - N4(0)" },
	{ "initN8",								GR_NODE,	false,	false,	0.0,	0,	"",					"ODE stuff - N8(0)" },
	/* INT */
	{ "nrSources",							GR_NODE,	true,	false,	0.0,	0,	"",					"Number of vascular sources on the grid" },
	{ "nrKillingsCaseation",				GR_NODE,	true,	false,	0.0,	0,	"",					"Number of killings for a compartment to become caseated" },
	{ "maxAge",								MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Maximal macrophage age" },
	{ "maxAgeAct",							MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Maximal activated macrophage age" },
	{ "initNumber",							MAC_NODE,	true,	false,	0.0,	0,	"",					"Initial number of resting macrophages on the grid" },
	{ "maxAge",								TCELL_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Maximal T cell age" },
	{ "timeRecEnabled",						TCELL_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time after which T cell recruitment is enabled" },
	{ "maxTimeReg",							TGAM_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time span during which a Tgam cell remains down-regulated" },
	{ "maxTimeReg",							TCYT_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time span during which a Tcyt cell remains down-regulated" },
	{ "movementRest",						MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time required for a resting macrophage to move one micro-compartment" },
	{ "movementInf",						MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time required for an infected macrophage to move one micro-compartment" },
	{ "movementAct",						MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time required for an active macrophage to move one micro-compartment" },
	{ "maxTimeReg",							MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time span during which a macrophage remains down-regulated" },
};

ParamsBase::ParamsBase(bool ode)
	: _ode(ode)
	, _doubleParam()
	, _intParam()
	, _initialMacs()
	, _initialExtMtb()
{
	memset(_doubleParam, 0, sizeof(double) * PARAM_DOUBLE_COUNT);
	memset(_intParam, 0, sizeof(int) * PARAM_INT_COUNT);


}

ParamsBase::~ParamsBase()
{
}

void ParamsBase::defineDefaults()
{
	// For those parameters that use defaults, set their parameter values to their defaults,
	// in case they are not specified in the parameter file.

	for(int i = 0; i < _PARAM_COUNT; i++)
	{
		if (_description[i].useDefault)
		{
			if (isDouble(i))
			{
				_doubleParam[i] = _description[i].doubleDefault;
			}
			else
			{
				_intParam[intIndex(i)] = _description[i].intDefault;
			}
		}
	}
}

// Look for an entry in the parameter description array with the specified parameter name and XmlElement.
// If found return its array index, if not return -1.
int ParamsBase::findParameterDescription(std::string name, const TiXmlElement* pElement) const
{
	XmlElement xmlElement = findElementDescription(pElement->Value());
	if (xmlElement == NULL_NODE)
	{
		return -1;
	}

	for(int i = 0; i < _PARAM_COUNT; i++)
	{
		if (_description[i].name == name && _description[i].xmlElement == xmlElement)
		{
			return i;
		}
	}

	return -1;
}

// Look for an entry in the node array with the specified node name.
// If found return its array index, if not return -1.
XmlElement ParamsBase::findElementDescription(std::string name) const
{
	for(int i = 0; i < NODE_COUNT; i++)
	{
		if (_element[i].name == name)
		{
			return (XmlElement) i;
		}
	}

	return NULL_NODE;
}

bool ParamsBase::readParam(const TiXmlElement* pElement, const std::string paramName, double* pVar, bool prob)
{
	switch (pElement->QueryDoubleAttribute(paramName.c_str(), pVar))
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

bool ParamsBase::readParam(const TiXmlElement* pElement, const std::string paramName, double* pVar, double defaultVal, bool prob)
{
	switch (pElement->QueryDoubleAttribute(paramName.c_str(), pVar))
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

bool ParamsBase::readParam(const TiXmlElement* pElement, const std::string paramName, int* pVar, bool pos)
{
	switch (pElement->QueryIntAttribute(paramName.c_str(), pVar))
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

bool ParamsBase::readParam(const TiXmlElement* pElement, const std::string paramName, int* pVar, int defaultVal, bool pos)
{
	switch (pElement->QueryIntAttribute(paramName.c_str(), pVar))
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

// Read the "Init" element, which specifies the simulation initial conditions.
// It is optional. If not specified hard coded initial conditions are used.
bool ParamsBase::readInitElement(const TiXmlElement* pRootElement)
{
	assert(pRootElement);

	_initialMacs.clear();
	_initialExtMtb.clear();

	bool res = true;

	const TiXmlElement* pInitElement = pRootElement->FirstChildElement("Init");
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
		// No "Init" element was specified.
		// Put one infected macrophage in the middle of the grid.
		Pos pos(NROWS/2, NCOLS/2);
		_initialMacs.push_back(pos);
	}

	return res;
}

bool ParamsBase::fromXml(const char* filename)
{
	bool res;

	if (!_xmlDoc.LoadFile(filename))
	{
		std::cerr << _xmlDoc.ErrorDesc() << " '" << filename << "'" << std::endl;
		return false;
	}

	const TiXmlElement* pRootElement = _xmlDoc.RootElement();
	if (!pRootElement || strcmp(pRootElement->Value(), _element[0].name.c_str()))
	{
		std::cerr << "Expected root element with name '" << _element[0].name << "'." << std::endl;
		return false;
	}

	defineDefaults();

	// Used to determine whether required parameters were specified or not.
	for (int i = 0; i < _PARAM_COUNT; i++)
	{
		_paramsRead[i] = false;
	}

	// Read the parameters starting at the root element.
	res = readElement(pRootElement, _paramsRead);

	// Check that required parameters were read.
	for (int i = 0; i < _PARAM_COUNT; i++)
	{
		if (!_description[i].useDefault && !_paramsRead[i])
		{
			res = false;
			std::cerr << "Required parameter " << _description[i].name << " is missing." << std::endl;
		}
	}

	// Read the "Init" element which specifies the simulation initial conditions.
	readInitElement(pRootElement);

	std::cerr << std::endl << "ParamsBase::fromXml, res: " << res << std::endl << std::endl;

	#if 0
	for (int i = 0; i < PARAM_DOUBLE_COUNT; i++)
	{
		std::cerr << "i: " << i << " " << _description[i].name << " " << _doubleParam[i] << std::endl;
	}
	std::cerr << std::endl;
	for (int i = PARAM_DOUBLE_COUNT; i < PARAM_DOUBLE_COUNT + PARAM_INT_COUNT; i++)
	{
		std::cerr << "i: " << i << " " << _description[i].name << " "  << _intParam[intIndex(i)] << std::endl;
	}
	#endif

	return res;
}

// Read the parameters for an XML node and for all its children.
bool ParamsBase::readElement(const TiXmlElement* pElement, bool paramsRead[])
{
	bool res = true;

	// Read parameters from the element attributes.
	// If an attribute name does not match any parameter name then ignore it.
	// If an attribute name does match a parameter name then read it.
	const TiXmlAttribute* pAttrib=pElement->FirstAttribute();
	while (pAttrib)
	{
		std::string attributeName = std::string(pAttrib->Name());

		int parameterIndex = findParameterDescription(attributeName, pElement);
		//std::cerr << "attributeName: " << attributeName << " parameterIndex: " << parameterIndex  << " isDouble: " << isDouble(parameterIndex) << " " << pAttrib->Value() << std::endl; //DBG

		if (parameterIndex >= 0)
		{
			paramsRead[parameterIndex] = true;

			if (isDouble(parameterIndex))
			{
				res &= readParam(pElement, pAttrib, (ParamDoubleType) parameterIndex);
			}
			else
			{
				res &= readParam(pElement, pAttrib, (ParamIntType) parameterIndex);
			}
		}

		pAttrib=pAttrib->Next();
	}

	//DBG
	#if 0
	std::cerr << std::endl << " res: " << res << std::endl;
	for (int i = 0; i < PARAM_DOUBLE_COUNT; i++)
	{
		std::cerr << "i: " << i << " " << _description[i].name << " " << _doubleParam[i] << std::endl;
	}
	std::cerr << std::endl;
	for (int i = PARAM_DOUBLE_COUNT; i < PARAM_DOUBLE_COUNT + PARAM_INT_COUNT; i++)
	{
		std::cerr << "i: " << i << " " << _description[i].name << " "  << _intParam[intIndex(i)] << std::endl;
	}
	#endif
	//DBG

	// Read the parameters from the children of this element that are themselves XML elements.
	const TiXmlNode* child = 0;
	while( (child = pElement->IterateChildren( child )) )
	{
//		std::cerr << "child: " << child->Value() << std::endl; //DBG

		// Ignore children that are not XML elements.
		if (child->ToElement())
		{
//			std::cerr << "element: " << child->Value() << std::endl; //DBG
			res &= readElement(child->ToElement(), paramsRead);
		}
	}

	return res;
}


// Define any computed parameters including any that are computed based on other parameters.
void ParamsBase::computeParams()
{
	std::cerr << " ParamsBase::computeParams: " << std::endl; //DBG
	defineRecruitmentWeight(PARAM_GR_WEIGHT_CCL2_RECRUITMENT, PARAM_MAC_SEC_RATE_CCL2);
	defineRecruitmentWeight(PARAM_GR_WEIGHT_CCL5_RECRUITMENT, PARAM_MAC_SEC_RATE_CCL5);
	defineRecruitmentWeight(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT, PARAM_MAC_SEC_RATE_CXCL9);
}

void ParamsBase::defineRecruitmentWeight(ParamDoubleType recruitmentWeightParam, ParamDoubleType secretionParam)
{
	std::cerr << " recruitmentWeightParam: " << recruitmentWeightParam << " read: " << _paramsRead[recruitmentWeightParam] << std::endl; //DBG
	if (!_paramsRead[recruitmentWeightParam])
	{
		std::cerr << " ParamsBase::defineRecruitmentWeight: "<< _description[recruitmentWeightParam].name << std::endl;//DBG
		std::cerr << " PARAM_MAC_SEC_RATE_TNF: " << getParam(PARAM_MAC_SEC_RATE_TNF) << " secretionParam: " << getParam(secretionParam) << std::endl; //DBG
		setParam(recruitmentWeightParam, getParam(PARAM_MAC_SEC_RATE_TNF) / getParam(secretionParam));
	}
}

// Perform any additional checks in addition to those performed by the readParam functions.
// These checks can't be included in the code that reads specific parameters because that
// code also gets called by the Lhs sub-class of ParamsBase. Class Lhs reads an LHS parameter
// file which has parameter ranges, not individual parameter values, so the double and int
// parameter value arrays are not defined so these checks will typically fail in that case.
// However class Lhs does call this function for each parameter file it generates.
bool ParamsBase::checkParams() const
{
	bool res = true;

	res &= movementBonusFactorCheck(PARAM_MAC_MOVEMENT_BONUSFACTOR, "macrophage");
	res &= movementBonusFactorCheck(PARAM_TCELL_MOVEMENT_BONUSFACTOR, "T cell");

	return res;
}

bool ParamsBase::movementBonusFactorCheck(ParamDoubleType param, std::string s) const
{
	if (getParam(param) < 1.0)
	{
		std::cerr << "The " << s << " movement bonus factor of " << getParam(param) << " is < 1.0." << std::endl;
		return false;
	}

	return true;
}

bool ParamsBase::toXml(const char* filename) const
{
	std::ofstream outFile(filename);
	if (!outFile.good())
	{
		return false;
	}

	const TiXmlElement* pRootElement = _xmlDoc.RootElement();
	if (!pRootElement || strcmp(pRootElement->Value(), _element[0].name.c_str()))
	{
		std::cerr << "Expected root element with name '" << _element[0].name << "'." << std::endl;
		return false;
	}

	// Write the parameters starting at the root element.
	// This only writes parameters that were present in the parameter file that was read.
	// When an LHS parameter file is read this writes a generated parameter file but only
	// for those parameters that had ranges specified in the LHS parameter file.
	writeElement(outFile, pRootElement, 0);

	return true;
}

void ParamsBase::writeElement(std::ostream& out, const TiXmlElement* pElement, int indent) const
{

	// Write the opening tag for the XML element.
	doIdentation(out, indent);
	out << "<" << pElement->Value() << std::endl;

	// Write parameters for the element attributes.
	// If an attribute name does not match any parameter name then ignore it.
	// If an attribute name does match a parameter name then write it.
	const TiXmlAttribute* pAttrib=pElement->FirstAttribute();
	int attribCount = 0;
	while (pAttrib)
	{
		std::string attributeName = std::string(pAttrib->Name());

		int parameterIndex = findParameterDescription(attributeName, pElement);

		if (parameterIndex >= 0)
		{
			writeParameter(out, parameterIndex, indent);
		}

		attribCount++;
		pAttrib=pAttrib->Next();
	}

	// Count the number of element children.
	const TiXmlNode* child = 0;
	int childCount = 0;
	while( (child = pElement->IterateChildren( child )) )
	{
		if (child->ToElement())
		{
			childCount++;
		}
	}

	// Close the XML element opening tag.
	if (attribCount > 0)
	{
		// The start tag close is indented the same as the attributes were.
		doIdentation(out, indent+1);
	}

	if (childCount > 0)
	{
		out << ">" << std::endl;
	}
	else
	{
		out << "/>" << std::endl;
	}

	// Write the parameters from the children of this element that are themselves XML elements.
	child = 0;
	while( (child = pElement->IterateChildren( child )) )
	{
//		std::cerr << "child: " << child->Value() << std::endl; //DBG

		// Ignore children that are not XML elements. Also ignore the Init element, since it
		// has no attributes and its children are read into the _initialMacs and _initialExtMtb lists,
		// which are processed later.
		if (child->ToElement() && strcmp(child->Value(), _element[INIT_NODE].name.c_str()))
		{
			writeElement(out, child->ToElement(), indent+1);
		}
	}

	// For the root element, write the "Init" element.
	if (pElement == _xmlDoc.RootElement())
	{
		writeInitNode(out, indent+1);
	}

	// Write the element closing tag, if needed.
	if (childCount > 0)
	{
		doIdentation(out, indent);
		out << "</" << pElement->Value() << ">" << std::endl;;
	}
}

void ParamsBase::writeParameter(std::ostream& out, int parameterIndex, int indent) const
{
	// Parameters are indented relative to the start tag for the XML element they are attributes of.
	doIdentation(out, indent+1);
	out << _description[parameterIndex].name << " = \"";

	if (isDouble(parameterIndex))
	{
		out << getParam((ParamDoubleType) parameterIndex);
	}
	else
	{
		out << getParam((ParamIntType) intIndex(parameterIndex));
	}

	out << "\""<< std::endl;
}

void ParamsBase::writeInitNode(std::ostream& out, int indent) const
{

	doIdentation(out, indent);
	out << "<Init>" << std::endl;
	for (PosVector::const_iterator it = _initialMacs.begin(); it != _initialMacs.end(); it++)
	{
		doIdentation(out, indent+1);
		out << "<Mac row = \"" << it->first << "\" col = \"" << it->second << "\"/>" << std::endl;
	}
	for (PosVector::const_iterator it = _initialExtMtb.begin(); it != _initialExtMtb.end(); it++)
	{
		doIdentation(out, indent+1);
		out << "<ExtMtb row = \"" << it->first << "\" col = \"" << it->second << "\"/>" << std::endl;
	}
	doIdentation(out, indent);
	out << "</Init>" << std::endl;
}

void ParamsBase::doIdentation(std::ostream& out, int indent) const
{
	for (int i = 0; i < indent; i++)
		out << "\t";
}

