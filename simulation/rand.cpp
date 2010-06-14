/*
 * rand.cpp
 *
 *  Created on: Feb 16, 2010
 *      Author: mohammed
 */

#include "rand.h"

Rand::Rand(unsigned int	seed)
	: _eng()
	, _realEng()
	, _uniformDist01(0, 1)
	, _seed(seed)
{
	_eng.seed(seed);
	_realEng.seed(seed);
}

Rand::~Rand()
{
}

void Rand::serialize(std::ostream& out) const
{
	assert(out.good());

	out << _seed << std::endl;
	out << _realEng << std::endl;
	out << _eng << std::endl;
}

void Rand::deserialize(std::istream& in)
{
	assert(in.good());

	in >> _seed;
	in >> _realEng;
	in >> _eng;
}
