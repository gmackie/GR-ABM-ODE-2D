/*
 * params.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef PARAMS_H
#define PARAMS_H

#include "tinyxml/tinyxml.h"
#include "gr.h"

// introduce a convenient shorthand
#define _PARAM(a) Params::getInstance()->getParam((a))

typedef enum {
	PARAM_GR_D_TNF,
	PARAM_GR_D_CHEMOKINES,
	PARAM_GR_DEG_TNF,
	PARAM_GR_DEG_CHEMOKINES,
	PARAM_GR_MIN_CHEMOTAXIS,
	PARAM_GR_MAX_CHEMOTAXIS,
	PARAM_GR_THRESHOLD_APOPTOSIS_TNF,
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
	PARAM_MAC_NR_UPTAKE_RI_EXTMTB,
	PARAM_MAC_PROB_KILL_R_EXTMTB,
	PARAM_MAC_THRESHOLD_NFKB_EXTMTB,
	PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB,
	PARAM_MAC_THRESHOLD_BURST_CI_INTMTB,
	PARAM_MAC_PROB_STAT1_TGAM,
	PARAM_MAC_NR_UPTAKE_A_EXTMTB,
	PARAM_MAC_THRESHOLD_RECRUITMENT,
	PARAM_MAC_PROB_RECRUITMENT,
	PARAM_TCELL_PROB_RECRUITMENT,
	PARAM_TCELL_PROB_MOVE_TO_MAC,
	PARAM_TCELL_PROB_MOVE_TO_TCELL,
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
	PARAM_EXTMTB_GROWTH_RATE,
	PARAM_EXTMTB_UPPER_BOUND,
	PARAM_muMDC_LN,
	PARAM_sn4,
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
	PARAM_sn8,
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
	PARAM_DOUBLE_COUNT // dummy for the count
} ParamDoubleType;

typedef enum {
	PARAM_GR_NR_SOURCES,
	PARAM_GR_NR_KILLINGS_FOR_CASEATION,
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
	PARAM_INT_COUNT // dummy for the count
} ParamIntType;

class Params
{
private:
	static const char* _description[][4];
	static Params* _pInstance;
	bool _useRecruitmentWeights;
	bool _ode;
	double _doubleParam[PARAM_DOUBLE_COUNT];
	int _intParam[PARAM_INT_COUNT];
	PosVector _initialMacs;
	PosVector _initialExtMtb;

	enum ClosingTagType { CLOSE_START_TAG, CLOSE_END_TAG, CLOSE_NONE };

protected:
	Params(bool ode = false);
	virtual ~Params();
	bool readParam(const TiXmlElement* pElement, const char* paramName, double* pVar, bool prob);
	bool readParam(const TiXmlElement* pElement, const char* paramName, double* pVar, double defaultVal, bool prob);
	bool readParam(const TiXmlElement* pElement, const char* paramName, int* pVar, bool pos);
	bool readParam(const TiXmlElement* pElement, const char* paramName, int* pVar, int defaultVal, bool pos);
	virtual bool readParam(const TiXmlElement* pElement, ParamDoubleType param, bool prob);
	virtual bool readParam(const TiXmlElement* pElement, ParamDoubleType param, double defaultVal, bool prob);
	virtual bool readParam(const TiXmlElement* pElement, ParamIntType param, bool pos);
	virtual bool readParam(const TiXmlElement* pElement, ParamIntType param, int defaultVal, bool pos);
	bool readGRElement(const TiXmlElement* pGrElement);
	bool readMacElement(const TiXmlElement* pMacElement);
	bool readTcellElement(const TiXmlElement* pTcellElement);
	bool readTgamElement(const TiXmlElement* pTgamElement);
	bool readTcytElement(const TiXmlElement* pTcytElement);
	bool readTregElement(const TiXmlElement* pTregElement);
	bool readMtbElement(const TiXmlElement* pMtbElement);
	bool readInitElement(const TiXmlElement* pGrElement);
	void writeParam(std::ostream& out, ParamDoubleType param, int indentation, ClosingTagType close) const;
	void writeParam(std::ostream& out, ParamIntType param, int indentation, ClosingTagType close) const;
	void updateRecruitmentWeights();

public:
	void printCSV(bool header) const;
	bool getUseOde() const;
	bool fromXml(const char* filename);
	bool toXml(const char* filename) const;
	double getParam(ParamDoubleType param) const;
	int getParam(ParamIntType param) const;
	void setParam(ParamDoubleType param, double val);
	void setParam(ParamIntType param, int val);
	bool getUseRecruitmentWeights() const;
	PosVector getInitialMacs() const;
	PosVector getInitialExtMtb() const;
	const char* getName(ParamDoubleType param) const;
	const char* getName(ParamIntType param) const;
	const char* getType(ParamDoubleType param) const;
	const char* getType(ParamIntType param) const;
	const char* getUnit(ParamDoubleType param) const;
	const char* getUnit(ParamIntType param) const;
	const char* getDescription(ParamDoubleType param) const;
	const char* getDescription(ParamIntType param) const;
	static Params* getInstance(bool ode = false);
	static bool reinit(const char* filename);
};

inline bool Params::getUseOde() const
{
	return _ode;
}

inline Params* Params::getInstance(bool ode)
{
	if (!_pInstance)
		_pInstance = new Params(ode);

	return _pInstance;
}

inline bool Params::getUseRecruitmentWeights() const
{
	return _useRecruitmentWeights;
}

inline PosVector Params::getInitialMacs() const
{
	return _initialMacs;
}

inline PosVector Params::getInitialExtMtb() const
{
	return _initialExtMtb;
}

inline double Params::getParam(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _doubleParam[param];
}

inline int Params::getParam(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _intParam[param];
}

inline void Params::setParam(ParamDoubleType param, double val)
{
	assert(param != PARAM_DOUBLE_COUNT);
	_doubleParam[param] = val;
}

inline void Params::setParam(ParamIntType param, int val)
{
	assert(param != PARAM_INT_COUNT);
	_intParam[param] = val;
}

inline const char* Params::getType(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param][3];
}

inline const char* Params::getType(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[param + PARAM_DOUBLE_COUNT][3];
}

inline const char* Params::getName(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param][0];
}

inline const char* Params::getName(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[param + PARAM_DOUBLE_COUNT][0];
}

inline const char* Params::getUnit(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param][2];
}

inline const char* Params::getUnit(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[param + PARAM_DOUBLE_COUNT][2];
}

inline const char* Params::getDescription(ParamDoubleType param) const
{
	assert(param != PARAM_DOUBLE_COUNT);
	return _description[param][1];
}

inline const char* Params::getDescription(ParamIntType param) const
{
	assert(param != PARAM_INT_COUNT);
	return _description[param + PARAM_DOUBLE_COUNT][1];
}

#endif /* PARAMS_H */
