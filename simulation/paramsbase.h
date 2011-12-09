/*
 * params.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef PARAMS_H
#define PARAMS_H

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif //TIXML_USE_STL

#include "tinyxml/tinyxml.h"
#include "gr.h"

typedef enum {
	  NULL_NODE = -1
	, GR_NODE
	, MAC_NODE
	, TCELL_NODE
	, TGAM_NODE
	, TCYT_NODE
	, TREG_NODE
	, MTB_NODE
	, INIT_NODE
	, NODE_COUNT
} XmlElement;

struct NodeDescription
{
	XmlElement type;
	std::string name;
	std::string description;
};

// These are used as indices into the _doubleParam array and the _description array.
typedef enum {
	PARAM_GR_D_TNF,
	PARAM_GR_D_SHED_TNFR2,
	PARAM_GR_D_CHEMOKINES,
	PARAM_GR_DEG_TNF,
	PARAM_GR_DEG_CHEMOKINES,
	PARAM_GR_MIN_CHEMOTAXIS,
	PARAM_GR_MAX_CHEMOTAXIS,
	PARAM_TGAM_SEC_RATE_TNF,
	PARAM_TCYT_SEC_RATE_TNF,

	// molecular TNF-associated parameters
	PARAM_GR_K_SYNTH_MAC,
	PARAM_GR_K_SYNTH_TCELL,
	PARAM_GR_K_TACE_MAC,
	PARAM_GR_K_TACE_TCELL,
	PARAM_GR_KD1,
	PARAM_GR_KD2,
	PARAM_GR_K_ON1,
	PARAM_GR_K_ON2,
	PARAM_GR_K_INT1,
	PARAM_GR_K_INT2,
	PARAM_GR_K_SHED,
	PARAM_GR_K_REC1,
	PARAM_GR_K_REC2,
	PARAM_GR_K_T1,
	PARAM_GR_K_T2,
	PARAM_GR_K_DEG1,
	PARAM_GR_K_DEG2,
	PARAM_GR_MEAN_TNFR1_MAC,
	PARAM_GR_STD_TNFR1_MAC,
	PARAM_GR_MEAN_TNFR2_MAC,
	PARAM_GR_STD_TNFR2_MAC,
	PARAM_GR_MEAN_TNFR1_TCELL,
	PARAM_GR_STD_TNFR1_TCELL,
	PARAM_GR_MEAN_TNFR2_TCELL,
	PARAM_GR_STD_TNFR2_TCELL,
	// end of molecular TNF-associated parameters

	// intracellular NFkB signaling pathway parameters
	PARAM_GR_KN,
	PARAM_GR_KNN,
	PARAM_GR_MEAN_NFKB,
	PARAM_GR_ka,
	PARAM_GR_ki,
	PARAM_GR_k1,
	PARAM_GR_kA20,
	PARAM_GR_k2,
	PARAM_GR_k3,
	PARAM_GR_k4,
	PARAM_GR_q1,
	PARAM_GR_q2,
	PARAM_GR_c1,
	PARAM_GR_c3,
	PARAM_GR_c4,
	PARAM_GR_c5,
	PARAM_GR_a1,
	PARAM_GR_a2,
	PARAM_GR_a3,
	PARAM_GR_tp,
	PARAM_GR_c5a,
	PARAM_GR_c6a,
	PARAM_GR_i1,
	PARAM_GR_e2a,
	PARAM_GR_i1a,
	PARAM_GR_e1a,
	PARAM_GR_q1r,
	PARAM_GR_q2r,
	PARAM_GR_q2rr,
	PARAM_GR_c1r,
	PARAM_GR_epsilon1,
	PARAM_GR_epsilon2,
	PARAM_GR_c3rChem,
	PARAM_GR_c4Chem,
	PARAM_GR_c5Chem,
	PARAM_GR_e3Chem,
	PARAM_GR_c3rTNF,
	PARAM_GR_c4TNF,
	PARAM_GR_c5TNF,
	PARAM_GR_e3TNF,
	PARAM_GR_c1rrACT,
	PARAM_GR_c3rACT,
	PARAM_GR_c4ACT,
	PARAM_GR_c5ACT,
	PARAM_GR_ACT_THRESHOLD,
	PARAM_GR_ACT_K,
	PARAM_GR_c1rrIAP,
	PARAM_GR_c3rIAP,
	PARAM_GR_c4IAP,
	PARAM_GR_c5IAP,
	PARAM_GR_K_IAP,
	// End of intracellular NFkB signaling pathway parameters

	PARAM_GR_THRESHOLD_APOPTOSIS_TNF,
	PARAM_GR_K_APOPTOSIS,
	PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR,
	PARAM_GR_K_APOPTOSIS_MOLECULAR,
	PARAM_GR_K_APOPTOSIS_NFKB_MOLECULAR,
	PARAM_GR_PROB_APOPTOSIS_TNF,
	PARAM_GR_WEIGHT_TNF_RECRUITMENT,
	PARAM_GR_WEIGHT_CCL2_RECRUITMENT,
	PARAM_GR_WEIGHT_CCL5_RECRUITMENT,
	PARAM_GR_WEIGHT_CXCL9_RECRUITMENT,
	PARAM_GR_SEC_RATE_ATTRACTANT,
	PARAM_MAC_SEC_RATE_CCL2,
	PARAM_MAC_SEC_RATE_CCL5,
	PARAM_MAC_SEC_RATE_CXCL9,
	PARAM_MAC_SEC_RATE_TNF,
	PARAM_MAC_THRESHOLD_NFKB_TNF,
	PARAM_MAC_K_NFKB,
	PARAM_MAC_THRESHOLD_NFKB_TNF_MOLECULAR,
	PARAM_MAC_K_NFKB_MOLECULAR,
	PARAM_MAC_NR_UPTAKE_RI_EXTMTB,
	PARAM_MAC_PROB_KILL_R_EXTMTB,
	PARAM_MAC_THRESHOLD_NFKB_EXTMTB,
	PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB,
	PARAM_MAC_THRESHOLD_BURST_CI_INTMTB,
	PARAM_MAC_PROB_STAT1_TGAM,
	PARAM_MAC_NR_UPTAKE_A_EXTMTB,
	PARAM_MAC_THRESHOLD_RECRUITMENT,
	PARAM_MAC_PROB_RECRUITMENT,
	PARAM_MAC_MOVEMENT_BONUSFACTOR,
	PARAM_MAC_INIT_DENSITY,
	PARAM_TCELL_PROB_RECRUITMENT,
	PARAM_TCELL_PROB_MOVE_TO_MAC,
	PARAM_TCELL_PROB_MOVE_TO_TCELL,
	PARAM_TCELL_MOVEMENT_BONUSFACTOR,
	PARAM_TGAM_PROB_APOPTOSIS_FAS_FASL,
	PARAM_TGAM_THRESHOLD_RECRUITMENT,
	PARAM_TGAM_PROB_RECRUITMENT,
	PARAM_TCYT_THRESHOLD_RECRUITMENT,
	PARAM_TCYT_PROB_RECRUITMENT,
	PARAM_TCYT_PROB_KILL_MAC,
	PARAM_TCYT_PROB_KILL_MAC_CLEANLY,
	PARAM_TREG_THRESHOLD_RECRUITMENT,
	PARAM_TREG_PROB_RECRUITMENT,
	PARAM_INTMTB_GROWTH_RATE,
	PARAM_INTMTB_GROWTH_RATE_FACTOR_POST_ADAPTIVE,
	PARAM_INTMTB_GROWTH_RATE_FACTOR_DELAY,
	PARAM_EXTMTB_GROWTH_RATE,
	PARAM_EXTMTB_UPPER_BOUND,
	PARAM_MTB_STOPPING_THRESHOLD,
	PARAM_muMDC_LN,
	PARAM_muN4,
	PARAM_k13,
	PARAM_hs13,
	PARAM_k14,
	PARAM_k15,
	PARAM_rho2,
	PARAM_k20a,
	PARAM_hs20a,
	PARAM_csi1,
	PARAM_csi1a,
	PARAM_muN8,
	PARAM_wT80,
	PARAM_k16,
	PARAM_hs16,
	PARAM_k17,
	PARAM_hs17,
	PARAM_k18,
	PARAM_rho3,
	PARAM_k24a,
	PARAM_hs24a,
	PARAM_csi2,
	PARAM_csi2a,
	PARAM_csi2b,
	PARAM_scaling,
	PARAM_m,
	PARAM_scaling_MDC,
	PARAM_scaling_LUNG,
	PARAM_scaling_LN,
	PARAM_initN4,
	PARAM_initN8,
	PARAM_SOURCE_DENSITY,
	PARAM_DOUBLE_COUNT // dummy for the count
} ParamDoubleType;

// These are used as indices into the _intParam array.
// These enum values + PARAM_DOUBLE_COUNT are used as indices into the _description array.
typedef enum {
	PARAM_GR_NR_SOURCES,
	PARAM_GR_NR_KILLINGS_FOR_CASEATION,
	PARAM_GR_NF_KB_TIME_COEFF,
	PARAM_MAC_AGE,
	PARAM_MAC_A_AGE,
	PARAM_MAC_INIT_NUMBER,
	PARAM_TCELL_AGE,
	PARAM_TCELL_TIME_RECRUITMENT_ENABLED,
	PARAM_TGAM_TIMESPAN_REGULATED,
	PARAM_TCYT_TIMESPAN_REGULATED,
	PARAM_MAC_MOVEMENT_RESTING,
	PARAM_MAC_MOVEMENT_INFECTED,
	PARAM_MAC_MOVEMENT_ACTIVE,
	PARAM_MAC_TIMESPAN_REGULATED,
	PARAM_MTB_STOPPING_TIME_STEP,
	PARAM_AREA_CELL_DENSITY_STOPPING_THRESHOLD,
	PARAM_AREA_CELL_DENSITY_STOPPING_TIME_STEP,
	PARAM_INT_COUNT // dummy for the count
} ParamIntType;

struct ParamDescription
{
	std::string name;

	XmlElement xmlElement;

	// For a double parameter whether or not it should be a probability (between 0 and 1).
	// For an int parameter whether or not it should be a positive value (greater than 0).
	bool probPos;

	// Whether or not the parameter should be read with a default value or not.
	// If so, then if not specified in a parameter file the corresponding default value is used -
	// doubleDefault for a double parameter, intDefault for an int parameter.
	// If not, then it is considered a required parameter and if not present in the parameter file
	// it is considered an error.
	bool useDefault;
	double doubleDefault;
	int intDefault;

	std::string unit;
	std::string description;
};

class ParamsBase
{

protected:
	static const NodeDescription _element[NODE_COUNT];
	static const int _PARAM_COUNT;
	static const ParamDescription _description[];

	TiXmlDocument _xmlDoc;
	bool _ode;

	bool _paramsRead[PARAM_DOUBLE_COUNT + PARAM_INT_COUNT];
	double _doubleParam[PARAM_DOUBLE_COUNT];
	int _intParam[PARAM_INT_COUNT];
	PosVector _initialMacs;
	PosVector _initialExtMtb;

	ParamsBase(bool ode);
	virtual ~ParamsBase();

	void defineDefaults();

	bool readElement(const TiXmlElement* pElement, bool paramsRead[]);
	virtual bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamDoubleType param) = 0;
	virtual bool readParam(const TiXmlElement* pElement, const TiXmlAttribute* pAttrib,  ParamIntType param) = 0;

	bool readParam(const TiXmlElement* pElement, const std::string paramName, double* pVar, bool prob);
	bool readParam(const TiXmlElement* pElement, const std::string paramName, double* pVar, double defaultVal, bool prob);
	bool readParam(const TiXmlElement* pElement, const std::string paramName, int* pVar, bool pos);
	bool readParam(const TiXmlElement* pElement, const std::string paramName, int* pVar, int defaultVal, bool pos);
  bool readParam(const TiXmlElement*, std::string, Pos&);

	bool readInitElement(const TiXmlElement* pRootElement);
	void computeParams();
	void defineRecruitmentWeight(ParamDoubleType recruitmentWeightParam, ParamDoubleType secretionParam);
	virtual bool checkParams() const;
	bool movementBonusFactorCheck(ParamDoubleType param, std::string s) const;

	void writeElement(std::ostream& out, const TiXmlElement* pElement, int indent) const;
	void writeParameter(std::ostream& out, int parameterIndex, int indent) const;
	void writeInitNode(std::ostream& out, int indent) const;
	void doIdentation(std::ostream& out, int indent) const;

public:

	static const std::string& getXmlElementName(XmlElement element);
	static const std::string& getXmlElementDescription(XmlElement element);

	bool getUseOde() const;
	bool fromXml(const char* filename);
	bool toXml(const char* filename) const;
	const TiXmlDocument& getXmlDoc() const;
	double getParam(ParamDoubleType param) const;
	int getParam(ParamIntType param) const;
	void setParam(ParamDoubleType param, double val);
	void setParam(ParamIntType param, int val);
	PosVector getInitialMacs() const;
	PosVector getInitialExtMtb() const;
	int findParameterDescription(std::string name, const TiXmlElement* pElement) const;
	XmlElement findElementDescription(std::string name) const;
	bool isDouble(int param) const;
	ParamIntType intIndex(int i) const;
	int paramIndex(ParamIntType intIndex) const;
	const std::string& getName(ParamDoubleType param) const;
	const std::string& getName(ParamIntType param) const;
	const XmlElement& getXmlElement(ParamDoubleType param) const;
	const XmlElement& getXmlElement(ParamIntType param) const;
	const std::string& getUnit(ParamDoubleType param) const;
	const std::string& getUnit(ParamIntType param) const;
	const std::string& getDescription(ParamDoubleType param) const;
	const std::string& getDescription(ParamIntType param) const;
	const TiXmlElement* getRootElement() const;
};


inline const TiXmlDocument& ParamsBase::getXmlDoc() const
{
	return _xmlDoc;
}

// Given an index into the _description array for an integer parameter,
// return the index into the _intParam array for that parameter.
inline 	ParamIntType ParamsBase::intIndex(int paramIndex) const
{
	return (ParamIntType) (paramIndex - PARAM_DOUBLE_COUNT);
}

// Given an index into the _intParam array for an integer parameter,
// return the index into the _description array for that parameter.
inline int ParamsBase:: paramIndex(ParamIntType intIndex) const
{
	return intIndex + PARAM_DOUBLE_COUNT;
}

inline const std::string& ParamsBase::getXmlElementName(XmlElement element)
{
	return _element[element].name;
}

inline const std::string& ParamsBase::getXmlElementDescription(XmlElement element)
{
	return _element[element].description;
}

inline bool ParamsBase::getUseOde() const
{
	return _ode;
}

inline PosVector ParamsBase::getInitialMacs() const
{
	return _initialMacs;
}

inline PosVector ParamsBase::getInitialExtMtb() const
{
	return _initialExtMtb;
}

inline double ParamsBase::getParam(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _doubleParam[param];
}

inline int ParamsBase::getParam(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _intParam[param];
}

inline void ParamsBase::setParam(ParamDoubleType param, double val)
{
	assert(param != PARAM_DOUBLE_COUNT);
	_doubleParam[param] = val;
}

inline void ParamsBase::setParam(ParamIntType param, int val)
{
	assert(param != PARAM_INT_COUNT);
	_intParam[param] = val;
}

inline const XmlElement& ParamsBase::getXmlElement(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param].xmlElement;
}

inline const XmlElement& ParamsBase::getXmlElement(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[paramIndex(param)].xmlElement;
}

inline 	bool ParamsBase::isDouble(int param) const
{
	if (param < PARAM_DOUBLE_COUNT)
	{
		return true;
	}
	else
	{
		return false;
	}
}

inline const std::string& ParamsBase::getName(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param].name;
}

inline const std::string& ParamsBase::getName(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[paramIndex(param)].name;
}

inline const std::string& ParamsBase::getUnit(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param].unit;
}

inline const std::string& ParamsBase::getUnit(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[paramIndex(param)].unit;
}

inline const std::string& ParamsBase::getDescription(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param].description;
}

inline const std::string& ParamsBase::getDescription(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[paramIndex(param)].description;
}

inline const TiXmlElement* ParamsBase::getRootElement() const
{
	return _xmlDoc.RootElement();
}

#endif /* PARAMS_H */
