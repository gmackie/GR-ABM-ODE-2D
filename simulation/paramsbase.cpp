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

// For each XML element, its XML element type and text name as it appears in a parameter file.
// The first element of this array is considered the root XML element and must be present in the parameter file as the root XML element.
const NodeDescription ParamsBase::_element[NODE_COUNT] = {
   // type			name		Description
	{ GR_NODE, 		"GR",		"" },
	{ MAC_NODE, 	"Mac",		"Macrophage specific parameters" },
	{ TCELL_NODE, 	"Tcell",	"Tcell specific parameters" },
	{ TGAM_NODE, 	"Tgam",		"Tgam specific parameters" },
	{ TCYT_NODE, 	"Tcyt",		"Tcyt specific parameters" },
	{ TREG_NODE, 	"Treg",		"Treg specific parameters" },
	{ MTB_NODE, 	"Mtb",		"Mtb specific parameters" },
	{ INIT_NODE, 	"Init",		""}
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
	// Name									XmlElement	probPos	Default	double	int	  unit				description
	{ "diffusivityTNF",						GR_NODE,	false,	false,	0.0,	0,	"cm^2/s",			"TNF diffusivity" },
	{ "diffusivityShedTNFR2",				GR_NODE,	false,	false,	0.0,	0,	"cm^2/s",			"Shed sTNF/TNFR2 comples diffusivity" },
	{ "diffusivityChemokines",				GR_NODE,	false,	false,	0.0,	0,	"cm^2/s",			"Chemokine diffusivity" },
	{ "diffusivityIL10",					GR_NODE,	false,	false,	0.0,	0,	"cm^2/s",			"IL10 diffusivity" },	
    { "degRateTNF",							GR_NODE,	false,	false,	0.0,	0,	"",					"TNF degradation rate (per 6s)" },
    { "degRateIL10",					    GR_NODE,	false,	false,	0.0,	0,	"",					"IL10 degradation rate (per 6s)" },
	{ "degRateChemokines",					GR_NODE,	false,	false,	0.0,	0,	"",					"Chemokine degradation rate (per 6s)" },
	{ "minChemotaxis",						GR_NODE,	false,	false,	0.0,	0,	"#molecules",		"Chemotaxis sensivity range (lower bound)" },
	{ "maxChemotaxis",						GR_NODE,	false,	false,	0.0,	0,	"#molecules",		"Chemotaxis sensivity range (upper bound)" },
	{ "dTNF_Tgam",							GR_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of TNF by Tgam" },
	{ "dTNF_Tcyt",							GR_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of TNF by Tcyt" },
    { "dIL10_Tgam",							GR_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of IL10 by Tgam" },
    { "dIL10_Tcyt",							GR_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of IL10 by Tcyt" },
    { "dIL10_Treg",							GR_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of IL10 by Treg" },
    

	// molecular TNF-associated parameters
	{ "kSynthMac",							GR_NODE,	false,	false,	0.0,	0,	"#/cell.sec",		"Basal rate of mTNF synthesis by a macrophage" },
	{ "kSynthTcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell.sec",		"Basal rate of mTNF synthesis by a T cell" },
	{ "kTrans",                             GR_NODE,	false,	false,	0.0,	0,	"1/sec",            "Rate of translation of TNF mRNA to mTNF" },
    { "kRNAMac",                            GR_NODE,	false,	false,	0.0,	0,	"1/sec",            "Basal rate of TNF mRNA synthesis by a macrophage" },
	{ "kRNATcell",                          GR_NODE,	false,	false,	0.0,	0,	"1/sec",            "Basal rate of TNF mRNA synthesis bt a T cell" },
    { "kTaceMac",							GR_NODE,	false,	false,	0.0,	0,	"1/sec",			"Basal rate of mTNF release from a macrophage" },
	{ "kTaceTcell",							GR_NODE,	false,	false,	0.0,	0,	"1/sec",			"Basal rate of mTNF release from a T cell" },
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
	{ "meanTNFR1Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Mean density of TNFR1 on macrophages" },
	{ "stdTNFR1Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Standard deviation of density of TNFR1 on macrophages" },
	{ "meanTNFR2Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Mean density of TNFR2 on macrophages" },
	{ "stdTNFR2Mac",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Standard deviation of density of TNFR2 on macrophages" },
	{ "meanTNFR1Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Mean density of TNFR1 on T cells" },
	{ "stdTNFR1Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Standard deviation of density of TNFR1 on T cells" },
	{ "meanTNFR2Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Mean density of TNFR2 on T cells" },
	{ "stdTNFR2Tcell",						GR_NODE,	false,	false,	0.0,	0,	"#/cell",			"Standard deviation of density of TNFR2 on T cells" },
	// end of molecular TNF-associated parameters

    // molecular IL10 associated parameters
    { "IkSynthMacInf",                      GR_NODE,	false,	false,	0.0,	0,	"#/cell.s",		    "sIL10 synthesis rate for infected Macs" },	
    { "IkSynthMacAct",                      GR_NODE,	false,	false,	0.0,	0,	"#/cell.s",		    "sIL10 synthesis rate for activated Macs" },	
    { "IkSynthTcell",                       GR_NODE,	false,	false,	0.0,	0,	"#/cell.s",		    "sIL10 synthesis rate for T cells" },	
    { "IkD",                                GR_NODE,	false,	false,	0.0,	0,	"M",		        "sIL10/IL10R equilibrium dissociation rate constant" },	
    { "IkOn",							    GR_NODE,	false,	false,	0.0,	0,	"1/M.s",		    "sIL10/IL10R association rate constant" },
	{ "IkOff",							    GR_NODE,	false,	false,	0.0,	0,	"1/s",		        "sIL10/IL10R dissociation rate constant" },
	{ "IkT",                                GR_NODE,	false,	false,	0.0,	0,	"1/s",		        "sIL10/IL10R association rate constant" },
	{ "IkInt",                              GR_NODE,	false,	false,	0.0,	0,	"1/s",		        "sIL10/IL10R association rate constant" },
	{ "meanIL10RMac",                       GR_NODE,	false,	false,	0.0,	0,	"#/cell",		    "Mean density of IL10R on macrophages" },
	{ "stdIL10RMac",                        GR_NODE,	false,	false,	0.0,	0,	"#/cell",		    "Standard deviation of mean density of IL10R on macrophages" },
	{ "meanIL10RTcell",                     GR_NODE,	false,	false,	0.0,	0,	"#/cell",		    "Mean density of IL10R on T cells" },
	{ "stdIL10RTcell",                      GR_NODE,	false,	false,	0.0,	0,	"#/cell",		    "Standard deviation of mean density of IL10R on T cells" },
    { "Imod",                               GR_NODE,    false,  false,  0.0,    0,  "",                 "Scaling factor for coarse grained IL10 dynamics" },
    // end of molecular IL10 associated parameters

    // TNF and IL10 linking parameters
	{ "LinkRNADelta",                       GR_NODE,	false,	false,	0.0,	0,	"",                 "Logistic function parameter value - changes slope" },
	{ "LinkRNAGamma",                       GR_NODE,	false,	false,	0.0,	0,	"",		            "Logistic function parameter value - changes where IC50 occurs" },
    { "LinkRNAMod",                         GR_NODE,	false,	false,	0.0,	0,	"",                 "Modification factor for when IL10 ODEs are not turned on" },
    { "LinkLogAlpha",                       GR_NODE,	false,	false,	0.0,	0,	"log10(ng/mL)",		"Coarse grained TNF/IL10 dose dependence parameter alpha" },
    { "LinkLogBeta",                        GR_NODE,	false,	false,	0.0,	0,	"log10(ng/mL)",     "Coarse grained TNF/IL10 dose dependence parameter beta" },
    // end of TNF and IL10 linking parameters
    
	// intracellular NFkB signaling pathway parameters
	{ "KN",									GR_NODE,	false,	true,	0.0,	0,	"#/cell",			"number of IKKK molecules" },
	{ "KNN",								GR_NODE,	false,	true,	0.0,	0,	"#/cell",			"number of IKK molecules" },
	{ "meanNFkB",							GR_NODE,	false,	true,	0.0,	0,	"#/cell",			"avearge number of total NFkB molecules" },
	{ "ka",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IKKK activation rate" },
	{ "ki",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IKKK inactivation rate" },
	{ "k1",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IKKn activation rate" },
	{ "kA20",								GR_NODE,	false,	true,	0.0,	0,	"#/cell",			"Michaelis coeff. in TNFR activity attenuation" },
	{ "k2",									GR_NODE,	false,	true,	0.0,	0,	"#/cell",			"Michaelis coeff. in IKKa inactivation" },
	{ "k3",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IKKn inactivation rate" },
	{ "k4",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IKKi->IKKii and IKKii->IKKn transormation" },
	{ "q1",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"NF-kB binding at A20 and IkB gene promoters" },
	{ "q2",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB inducible NF-kB detaching from A20 and IkB genes" },
	{ "c1",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"inducible A20 and IkB mRNA synthesis" },
	{ "c3",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"A20 and IkB mRNA degradation" },
	{ "c4",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"A20 and IkB traslation" },
	{ "c5",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"A20 degradation rate" },
	{ "a1",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB association NF-kB" },
	{ "a2",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB phosphorylation" },
	{ "a3",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB phosphorylation in IkB/NF-kB complexes" },
	{ "tp",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"degradation of phosphorylated IkB" },
	{ "c5a",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"spontaneous IkB degradation" },
	{ "c6a",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"spontaneous IkB degradation in IkB/NF-kB complexes" },
	{ "i1",									GR_NODE,	false,	true,	0.0,	0,	"1/s",				"NF-kB nuclear import" },
	{ "e2a",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB/NF-kB nuclear export" },
	{ "i1a",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB nuclear import" },
	{ "e1a",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB nuclear export" },
	{ "q1r",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"NF-kB binding at reporter gene promoter" },
	{ "q2r",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IkB inducible NF-kB detaching from reporter gene" },
	{ "q2rr",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"spontaneous NF-kB detaching from reporter gene" },
	{ "c1r",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"inducible reporter mRNA synthesis" },
	{ "epsilon1",							GR_NODE,	false,	true,	0.0,	0,	"-"	,				"ratio for mRNA synthesis: NF-kB independent / NF-kB dependent" },
	{ "epsilon2",							GR_NODE,	false,	true,	0.0,	0,	"-"	,				"ratio for NF-kB dependent mRNA synthesis: resting Mac / infected Mac" },
	{ "c3rChem",							GR_NODE,	false,	true,	0.0,	0,	"1/s",				"chemokine mRNA degradation rate" },
	{ "c4Chem",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"chemokine translation" },
	{ "c5Chem",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"intracellualr chemokine degradation rate" },
	{ "e3Chem",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"chemokine export rate" },
	{ "c3rTNF",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"TNF mRNA degradation rate" },
	{ "c4TNF",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"TNF translation" },
	{ "c5TNF",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"intracellualr TNF degradation rate" },
	{ "e3TNF",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"TNF export rate" },
	{ "c1rrACT",							GR_NODE,	false,	true,	0.0,	0,	"1/s",				"mac-activation molecule mRNA constitutive synthesis" },
	{ "c3rACT",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"mac-activation molecule mRNA degradation rate" },
	{ "c4ACT",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"mac-activation molecule translation" },
	{ "c5ACT",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"mac-activation molecule degradation rate" },
	{ "actThreshod",						GR_NODE,	false,	true,	0.0,	0,	"-"	,				"ACT threshold for Mac activation" },
	{ "activationRate",						GR_NODE,	false,	true,	0.0,	0,	"1/s",				"mac activation rate constant after ACT threshold is met" },
	{ "c1rrIAP",							GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IAP mRNA constitutive synthesis" },
	{ "c3rIAP",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IAP mRNA degradation rate" },
	{ "c4IAP",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IAP translation" },
	{ "c5IAP",								GR_NODE,	false,	true,	0.0,	0,	"1/s",				"IAP degradation rate" },
	{ "kIAP",								GR_NODE,	false,	true,	0.0,	0,	"-",				"IAP inhibitory effect on k_apoptosis" },
	// End of intracellular NFkB signaling pathway parameters

	{ "thresholdApoptosisTNF",				GR_NODE,	false,	false,	0.0,	0,	"fraction",			"TNF threshold for TNF-induced apoptosis" },
	{ "kApoptosis",							GR_NODE,	false,	false,	0.0,	0,	"1/s",				"Rate of apoptosis happening" },
	{ "thresholdApoptosisTNF_Molecular",	GR_NODE,	false,	false,	0.0,	0,	"#molecules",		"TNF threshold for TNF-induced apoptosis" },
	{ "kApoptosis_Molecular",				GR_NODE,	false,	false,	0.0,	0,	"1/s",				"Rate of apoptosis happening" },
	{ "kApoptosis_NFkB_Molecular",			GR_NODE,	false,	false,	0.0,	0,	"1/s",				"Rate of apoptosis happening when NFkB dynamics are on" },
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
    { "dIL10",								MAC_NODE,	false,	false,	0.0,	0,	"#molecules/6s",	"Secretion rate of IL10" },
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
	{ "initDensity",						MAC_NODE,	true,	true,	0.0,	0,	"",					"Initial density of resting macrophages on the grid" },
	{ "thresholdICOS",						MAC_NODE,	false,	false,	0.0,	0,	"",					"Threshold for ICOS-L inhibition on a Mac by IL10" },
	{ "mtbRecEnabled",						TCELL_NODE,	false,	false,	0.0,	0,	"#mtb",				"The tot mtb count above which T cell recruitment is enabled" },
	{ "probRec",							TCELL_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a T cell" },
	{ "probMoveToMac",						TCELL_NODE,	true,	false,	0.0,	0,	"",					"Probability of a T cell moving onto a compartment already containing a macrophage" },
	{ "probMoveToTcell",					TCELL_NODE,	true,	false,	0.0,	0,	"",					"Probability of a T cell moving onto a compartment already containing another T cell" },
	{ "movementBonusFactor",				TCELL_NODE,	false,	true,	1.0,	1,	"",					"Bonus factor for neighbor compartment with highest chemical concentration" },
	{ "probApoptosisFasFasL",				TGAM_NODE,	true,	false,	0.0,	0,	"",					"Probability of Fas/FasL induced apoptosis by a Tgam cell" },
	{ "thresholdRec",						TGAM_NODE,	false,	false,	0.0,	0,	"",					"TNF/chemokine threshold for Tgam recruitment" },
	{ "probRec",							TGAM_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a Tgam cell" },
	{ "thresholdExtMtb",					TGAM_NODE,	false,	false,	0.0,	0,	"",					"Ext. Mtb threshold for successfull antigen presentation by mac" },
	{ "probAntigenPresentation",            TGAM_NODE,	true,	false,	0.0,	0,	"",					"Probability of successfull antigen presentation by infected mac" },
	{ "probTGFB",                           TGAM_NODE,	false,	false,	0.0,	0,	"",					"Probability rate of TGF-B Tgam10 activation" },
	{ "probICOS",                           TGAM_NODE,	false,	false,	0.0,	0,	"",					"Probability rate of ICOS Tgam10 activation" },
	{ "probAgDegree",                       TGAM_NODE,	false,	false,	0.0,	0,	"",					"Probability rate of 2degree antigen stimulation" },
    { "thresholdRec",						TCYT_NODE,	false,	false,	0.0,	0,	"",					"TNF/chemokine threshold for Tcyt recruitment" },
	{ "probRec",							TCYT_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a Tcyt cell" },
	{ "probKillMac",						TCYT_NODE,	true,	false,	0.0,	0,	"",					"Probability of a Tcyt cell killing a (chronically) infected macrophage" },
	{ "probKillMacCleanly",					TCYT_NODE,	true,	false,	0.0,	0,	"",					"Probability of a Tcyt cell killing a chronically infected macrophage cleanly" },
	{ "thresholdRec",						TREG_NODE,	false,	false,	0.0,	0,	"",					"TNF/chemokine threshold for Treg recruitment" },
	{ "probRec",							TREG_NODE,	true,	false,	0.0,	0,	"",					"Probability of recruiting a Treg cell" },
	{ "growthRateIntMtb",					MTB_NODE,	false,	false,	0.0,	0,	"",					"Growth rate of intracellular bacteria" },
	{ "growthRateFactorPostAdaptiveIntMtb",	MTB_NODE,	false,	true,	1.0,	1,	"",					"Factor applied to growth rate of intracellular bacteria (parameter growthRateIntMtb) after adaptive immunity begins" },
	{ "growthRateFactorDelayIntMtb",		MTB_NODE,	false,	true,	0.0,	0,	"",					"Delay, in time steps, after adaptive immunity begins (parameter timeRecEnabled), before applying parameter growthRateFactorPostAdaptiveIntMtb" },
	{ "growthRateExtMtb",					MTB_NODE,	false,	false,	0.0,	0,	"",					"Growth rate of extracellular bacteria" },
	{ "growthExtMtbBound",					MTB_NODE,	false,	false,	0.0,	0,	"#bacteria",		"Upper bound on the number of extracellular bacteria used in growth function" },
	{ "mtbStoppingThreshold",				GR_NODE,	false,	true,	0.0,	0,	"#bacteria", 		"Stop simulation if mtb count above this threshold and simulation is at the time step specified by mtbStopppingTimeStep"},
	{ "mtbStoppingThreshold2",				GR_NODE,	false,	true,	0.0,	0,	"#bacteria", 		"Stop simulation if mtb count above this threshold and simulation is at the time step specified by mtbStopppingTimeStep"},
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
	{ "sourceDensity",      				GR_NODE,	false,	true,	0.2, 	0,	"",					"Density of vascular sources on the grid"},

	/* INT */
	{ "nrSources",							GR_NODE,	true,	true,	0.0,	0,	"",					"Number of vascular sources on the grid" },
	{ "nrKillingsCaseation",				GR_NODE,	true,	false,	0.0,	0,	"",					"Number of killings for a compartment to become caseated" },
	{ "NFkBTimeCoefficient",				GR_NODE,	true,	false,	0.0,	0,	"-",				"number of NF-kB time-steps within a single diffusion time-step" },
	{ "maxAge",								MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Maximal macrophage age" },
	{ "maxAgeAct",							MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Maximal activated macrophage age" },
	{ "initNumber",							MAC_NODE,	true,	true,	0.0,	0,	"",					"Initial number of resting macrophages on the grid" },
	{ "maxAge",								TCELL_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Maximal T cell age" },
	{ "timeRecEnabled",						TCELL_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time after which T cell recruitment is enabled" },
	{ "maxTimeReg",							TGAM_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time span during which a Tgam cell remains down-regulated" },
    { "maxTimeReg",							TCYT_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time span during which a Tcyt cell remains down-regulated" },
	{ "movementRest",						MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time required for a resting macrophage to move one micro-compartment" },
	{ "movementInf",						MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time required for an infected macrophage to move one micro-compartment" },
	{ "movementAct",						MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time required for an active macrophage to move one micro-compartment" },
	{ "maxTimeReg",							MAC_NODE,	true,	false,	0.0,	0,	"#timesteps",		"Time span during which a macrophage remains down-regulated" },
	{ "mtbStopppingTimeStep",				GR_NODE,	true,	true,	0.0,	0,	"",					"Time step at which to check mtbStoppingThreshold" },
	{ "mtbStopppingTimeStep2",				GR_NODE,	true,	true,	0.0,	0,	"",					"Time step at which to check mtbStoppingThreshold" },
	{ "areaCellDensityStoppingThreshold",	GR_NODE,	false,	true,	0.0,	0,	"um^2",				"Stop simulation if area by cell density above this threshold and simulation is at the time step specified by areaCellDensityStopppingTimeStep"},
	{ "areaCellDensityStoppingThreshold2",	GR_NODE,	false,	true,	0.0,	0,	"um^2",				"Stop simulation if area by cell density above this threshold and simulation is at the time step specified by areaCellDensityStopppingTimeStep"},
	{ "areaCellDensityStopppingTimeStep",	GR_NODE,	true,	true,	0.0,	0,	"",					"Time step at which to check areaCellDensityStoppingThreshold" },
	{ "areaCellDensityStopppingTimeStep2",	GR_NODE,	true,	true,	0.0,	0,	"",					"Time step at which to check areaCellDensityStoppingThreshold" },
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
bool ParamsBase::readParam(const TiXmlElement* pElement, const std::string paramName, Pos& p) {
  std::string val;
  switch(pElement->QueryValueAttribute(paramName, &val))
  {
  case TIXML_NO_ATTRIBUTE: break;
  case TIXML_WRONG_TYPE:
		std::cerr << "Value of attribute '" << pElement->Value() << "/@"
			<< paramName << "' must be a list of integers or 'center'";
    break;
  default:
    if(val == std::string("center"))
      p = Pos(NROWS/2, NCOLS/2);
    else {
      std::stringstream ss(val);
      ss>>p; //throws parse error
    }
    return true;
  }
  return false;
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
        res = readParam(pInitChildElement, "pos", pos);
        if(res) {
          Pos offset;
          res |= readParam(pInitChildElement, "offset", offset);
          if(res) {
            pos.first += offset.first;
            pos.second += offset.second;
          }
        }
        else {
  				res = readParam(pInitChildElement, "row", &pos.first, true);
  				res &= readParam(pInitChildElement, "col", &pos.second, true);
        }

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
        res = readParam(pInitChildElement, "pos", pos);
        if(res) {
          Pos offset;
          res |= readParam(pInitChildElement, "offset", offset);
          if(res) {
            pos.first += offset.first;
            pos.second += offset.second;
          }
        }
        else {
  				res = readParam(pInitChildElement, "row", &pos.first, true);
  				res &= readParam(pInitChildElement, "col", &pos.second, true);
        }

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
	res &= readInitElement(pRootElement);

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

		if (parameterIndex >= 0)
		{
			paramsRead[parameterIndex] = true;

			if (isDouble(parameterIndex))
			{
				res &= readParam(pElement, pAttrib, (ParamDoubleType) parameterIndex);
			}
			else
			{
				res &= readParam(pElement, pAttrib, (ParamIntType) intIndex(parameterIndex));
			}
		}

		pAttrib=pAttrib->Next();
	}

	// Read the parameters from the children of this element that are themselves XML elements.
	const TiXmlNode* child = 0;
	while( (child = pElement->IterateChildren( child )) )
	{
		// Ignore children that are not XML elements.
		if (child->ToElement())
		{
			res &= readElement(child->ToElement(), paramsRead);
		}
	}

	return res;
}


// Define any computed parameters including any that are computed based on other parameters.
void ParamsBase::computeParams()
{
  if(_paramsRead[PARAM_MAC_INIT_DENSITY])
    setParam(PARAM_MAC_INIT_NUMBER, getParam(PARAM_MAC_INIT_DENSITY) * FLOAT_TYPE(NROWS*NCOLS));
  else if(_paramsRead[PARAM_MAC_INIT_NUMBER])
    setParam(PARAM_MAC_INIT_DENSITY, getParam(PARAM_MAC_INIT_NUMBER) / FLOAT_TYPE(NROWS*NCOLS));
  else
    throw std::runtime_error("Initial resting macs not specified in parameter file");

  if(_paramsRead[PARAM_SOURCE_DENSITY])
    setParam(PARAM_GR_NR_SOURCES, getParam(PARAM_SOURCE_DENSITY) * FLOAT_TYPE(NROWS*NCOLS));
  else if(_paramsRead[PARAM_GR_NR_SOURCES])
    setParam(PARAM_SOURCE_DENSITY, getParam(PARAM_GR_NR_SOURCES) / FLOAT_TYPE(NROWS*NCOLS));
  else
    throw std::runtime_error("Sources not specified in parameter file");

	defineRecruitmentWeight(PARAM_GR_WEIGHT_CCL2_RECRUITMENT, PARAM_MAC_SEC_RATE_CCL2);
	defineRecruitmentWeight(PARAM_GR_WEIGHT_CCL5_RECRUITMENT, PARAM_MAC_SEC_RATE_CCL5);
	defineRecruitmentWeight(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT, PARAM_MAC_SEC_RATE_CXCL9);
}

void ParamsBase::defineRecruitmentWeight(ParamDoubleType recruitmentWeightParam, ParamDoubleType secretionParam)
{
	if (!_paramsRead[recruitmentWeightParam])
	{
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
  const TiXmlElement* pInitNode = _xmlDoc.RootElement()->FirstChildElement("Init");
	out << "<Init>" << std::endl;
  for(const TiXmlElement* pElem = pInitNode->FirstChildElement(); pElem; pElem = pElem->NextSiblingElement()) {
    doIdentation(out, indent+1);
    out << '<' << pElem->Value();
    for(const TiXmlAttribute* pAtt = pElem->FirstAttribute(); pAtt; pAtt = pAtt->Next()) {
      out<< ' ' << pAtt->Name() << "=\"" << pAtt->Value() << "\" ";
    }
    out<<"/>"<<std::endl;
  }
	doIdentation(out, indent);
	out << "</Init>" << std::endl;
}

void ParamsBase::doIdentation(std::ostream& out, int indent) const
{
	for (int i = 0; i < indent; i++)
		out << "\t";
}


