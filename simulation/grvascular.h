/*
* grvascular.h
*
*
* Created by ncilfone on 9/25/12.
*/


#ifndef GRVASCULAR_H
#define GRVASCULAR_H

#include "grgrid.h"
#include "gr.h"

using namespace std;

class Vascular
{
private:
    //double _bloodConcentrationIL10;
    //double _bloodConcentrationTNF;
    double _bloodConcentrationINH;
    double _bloodConcentrationRIF;

    void updateBloodODE(double& bloodConcentration, double clearance, double dt);

    void addDose(const int time, double dosage, const int doseInterval, double& bloodConcentration);
    bool isDoseStartTime(const int time);
    void modulateBlood(double& bloodconc, double change, Scalar MW);
    void solveMolecule();
    double calculateFluxChange(double localconc, double bloodconc, double dt, double permeability, Scalar MW);
public:
    Vascular();
    virtual ~Vascular();
    void solveVascularSources(GrGrid& grid, double dt, const int time, int DiffStep);
    double getAverageOfNeighbours(Scalar* grid, const Pos& pos, const Pos& dim, const Scalar MW);
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version);
    double getBloodConcentrationINH();
    double getBloodConcentrationRIF();

    Vascular* clone() const
    {
      return new Vascular(*this);
    }
};

inline void Vascular::addDose(const int time, double dosage, const int doseInterval, double& bloodConcentration)
{
    int timeSinceStart = time - _PARAM(_dosageStartTime)*600;
    bloodConcentration = (((timeSinceStart + doseInterval*600) % (doseInterval*600)) == 0) ? (((bloodConcentration * NAV * _PARAM(_bloodVolume)) + dosage)/(NAV * _PARAM(_bloodVolume))) : bloodConcentration;
    //std::cout << "at time:  "<< time << "blood conc is: "<< bloodConcentration << std::endl;
    //std::cout << "time since start dose is  "<< timeSinceStart << "check is: " << ((timeSinceStart + doseInterval*600) % (doseInterval*600))<< std::endl;
}

inline bool Vascular::isDoseStartTime(const int time)
{
    return time >= _PARAM(_dosageStartTime)*600;
}

inline double Vascular::calculateFluxChange(double localconc, double bloodconc, double dt, double permeability, Scalar MW)
{
    const double vascularSurfaceArea = (20e-4*20e-4*4.0);
    double change = permeability * vascularSurfaceArea * (bloodconc - localconc) * dt * MW;
    assert(abs(change)<=abs((bloodconc-localconc))*MW);
    return change;   // mg/L
}

inline void Vascular::modulateBlood(double &bloodconc, double change, Scalar MW)
{
    double change_molecules = change * VOL * NAV / MW;
    bloodconc = ((bloodconc * NAV * _PARAM(_bloodVolume)) - change_molecules)/(NAV*_PARAM(_bloodVolume));
}

inline double Vascular::getBloodConcentrationINH()
{
    return _bloodConcentrationINH;
}

inline double Vascular::getBloodConcentrationRIF()
{
    return _bloodConcentrationRIF;
}

template<class Archive>
void Vascular::serialize(Archive& ar, const unsigned int /*version*/)
{
    //ar & BOOST_SERIALIZATION_NVP(_bloodConcentrationIL10);
    //ar & BOOST_SERIALIZATION_NVP(_bloodConcentrationTNF);
    ar & BOOST_SERIALIZATION_NVP(_bloodConcentrationINH);
    ar & BOOST_SERIALIZATION_NVP(_bloodConcentrationRIF);
}

#endif /* GRVASCULAR_H */
