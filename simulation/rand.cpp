/*
 * rand.cpp
 *
 *  Created on: Feb 16, 2010
 *      Author: mohammed
 */

#include "rand.h"

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

	out << getSerialSize() << std::endl;
	out << _seed << std::endl;
	out << _eng << std::endl;
	out << _realEng << std::endl;

}

void Rand::deserialize(std::istream& in)
{
	assert(in.good());

	// This check isn't fool proof because of data alignment in memory.
	// An object's size can be bigger than the sum of the size's of its members.
	// When a new member is added or an existing one deleted the object size can remain unchanged.
	std::size_t currentSerialSize = getSerialSize();
	std::size_t savedSerialSize;
	in >> savedSerialSize;
	if (savedSerialSize != currentSerialSize)
	{
		std::cerr << "Error deserializing Rand object."<< std::endl;
		std::cerr << "The saved serial size of " << savedSerialSize << " does not match the current serial size of " << currentSerialSize << std::endl;
		exit(1);
	}

	in >> _seed;
	in >> _eng;
	in >> _realEng;

}

// For debugging problems with random number generation being out of synch between
// supposedly identical runs (same model, same parameter file, same seed).
void Rand::test(int time)
{
	std::cerr  << "ts: " << time << " rf: " << getReal() << " ri: " << getInt(101, 1) << std::endl;
}
