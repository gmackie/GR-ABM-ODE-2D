/*
 * recruitmentprob.h
 *
 *  Created on: Jun 21, 2010
 *      Author: mohammed
 */

#ifndef RECRUITMENTPROB_H_
#define RECRUITMENTPROB_H_

#include "gr.h"
#include "grgrid.h"
#include "stat.h"
#include "recruitmentbase.h"

class RecruitmentProb : public RecruitmentBase
{
private:
	void recruitMac(GrSimulation& sim, const Pos& pSource);
	void recruitTcell(GrSimulation& sim, const Pos& pSource);

public:
	RecruitmentProb();
	virtual ~RecruitmentProb();

	RecruitmentMethod getMethod() const;
	void serialize(std::ostream&) const;
	void deserialize(std::istream&);

	void recruit(GrSimulation& sim);

  /*virtual*/ RecruitmentBase* clone() const { return new RecruitmentProb(*this); }
};

inline RecruitmentMethod RecruitmentProb::getMethod() const
{
	return RECR_PROB;
}

inline void RecruitmentProb::serialize(std::ostream&) const
{
}

inline void RecruitmentProb::deserialize(std::istream&)
{
}

#endif /* RECRUITMENTPROB_H_ */
