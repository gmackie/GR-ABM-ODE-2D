/*
 * rand.cpp
 *
 *  Created on: Feb 16, 2010
 *      Author: mohammed
 */

#include "rand.h"
#include "serialization.h"

const std::string Rand::_ClassName = "Rand";

Rand::Rand(unsigned int	seed)
	: _seed(seed)
	, _eng()
	, _realEng()
	, _uniformDist01(0, 1)
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

	Serialization::writeHeader(out, Rand::_ClassName);

	out << _seed << std::endl;
	out << _eng << std::endl;
	out << _realEng << std::endl;

	Serialization::writeFooter(out, Rand::_ClassName);
}

void Rand::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Rand::_ClassName))
	{
		exit(1);
	}

	in >> _seed;
	in >> _eng;
	in >> _realEng;

	if (!Serialization::readFooter(in, Rand::_ClassName))
	{
		exit(1);
	}
}

// For debugging problems with random number generation being out of synch between
// supposedly identical runs (same model, same parameter file, same seed).
void Rand::test(int time)
{
	std::cerr  << "ts: " << time << " rf: " << getReal() << " ri: " << getInt(101, 1) << std::endl;
}
