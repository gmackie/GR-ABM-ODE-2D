/*
 * grsimulation.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#include <boost/archive/codecvt_null.hpp>
#include <fstream>
#include <float.h>
#include "grsimulation.h"
#include "grdiffusion.h"
#include "grdiffusionbtcs.h"
#include "grdiffusionwrongbtcs.h"
#include "grdiffusionadeswap.h"
#include "grdiffusionftcsswap.h"
#include "grdiffusionadeswap.h"
#include "grvascular.h"
#include "areatest.h"
#include "mtbtest.h"
#include "recruitmentprob.h"
#include "recruitmentlnodepure.h"

using namespace std;

GrSimulation::GrSimulation(const Pos& dim)
    : _time(0)
    , _grid(dim)
    , _macList()
    , _tgamList()
    , _tcytList()
    , _tregList()
    , _statsPrevious()
    , _stats()
    , _areaThreshold(0.5f)
    , _areaThresholdCellDensity(0.5f)
    , _pDiffusion(new GrDiffusionADE_Swap())
    , _pTTest()
    , _pRecruitment(NULL)
    , _pVascular(new Vascular())
    , _tnfrDynamics(false)
    , _nfkbDynamics(false)
    , _il10rDynamics(false)
    , _tgammatransition(false)
    , _odeSolver(ODESolvers::NMethods)
    , _tnfDepletionTimeStep(-1)
    , _il10DepletionTimeStep(-1)
    , _tcellRecruitmentBegun(false)
    , _numMolecularPerDiffusion(0)
    , _numDiffusionPerAgent(0)
{
    for (int i = 0; i < NOUTCOMES; i++)
    {
        _pTTest[i] = NULL;
    }

    _numMolecularPerDiffusion = (int)(_PARAM(_timestepDiffusion)/_PARAM(_timestepMolecular)); // Number of molecular iterations per diffusion iteration
    _numDiffusionPerAgent = (int)(AGENT_TIME_STEP / _PARAM(_timestepDiffusion)); // Number of diffusion iterations per agent iteration
}

GrSimulation* GrSimulation::clone(GrSimulation* sim) const {
    assert(sim != this);  //Don't copy to yourself... please...
    GrSimulation* ret = NULL;
    if(sim == NULL)
        ret = new GrSimulation(*this);
    else {
        sim->~GrSimulation(); //Call the destructor *ONLY* in the case of placement new
        ret = new (sim) GrSimulation(*this);
    }
    for(MacList::iterator begin = ret->_macList.begin(); begin!=ret->_macList.end(); begin++)
    {
        Mac* pAgent = static_cast<Mac*>((*begin)->clone());
        for(size_t i=0;i<GrGrid::MAX_AGENTS_PER_CELL;i++) //Place it proper in the grid
            if(*begin == getGrid().agent(pAgent->getPosition(), i)) {
                ret->getGrid().agent(pAgent->getPosition(), i) = pAgent;
                break;
            }
        *begin = pAgent;
    }
    for(TgamList::iterator begin = ret->_tgamList.begin(); begin!=ret->_tgamList.end(); begin++)
    {
        Tgam* pAgent = static_cast<Tgam*>((*begin)->clone());
        for(size_t i=0;i<GrGrid::MAX_AGENTS_PER_CELL;i++) //Place it proper in the grid
            if(*begin == getGrid().agent(pAgent->getPosition(), i)) {
                ret->getGrid().agent(pAgent->getPosition(), i) = pAgent;
                break;
            }
        *begin = pAgent;
    }
    for(TregList::iterator begin = ret->_tregList.begin(); begin!=ret->_tregList.end(); begin++)
    {
        Treg* pAgent = static_cast<Treg*>((*begin)->clone());
        for(size_t i=0;i<GrGrid::MAX_AGENTS_PER_CELL;i++) //Place it proper in the grid
            if(*begin == getGrid().agent(pAgent->getPosition(), i)) {
                ret->getGrid().agent(pAgent->getPosition(), i) = pAgent;
                break;
            }
        *begin = pAgent;
    }
    for(TcytList::iterator begin = ret->_tcytList.begin(); begin!=ret->_tcytList.end(); begin++)
    {
        Tcyt* pAgent = static_cast<Tcyt*>((*begin)->clone());
        for(size_t i=0;i<GrGrid::MAX_AGENTS_PER_CELL;i++) //Place it proper in the grid
            if(*begin == getGrid().agent(pAgent->getPosition(), i)) {
                ret->getGrid().agent(pAgent->getPosition(), i) = pAgent;
                break;
            }
        *begin = pAgent;
    }
    if(ret->_pDiffusion)
        ret->_pDiffusion = ret->_pDiffusion->clone();
    if(ret->_pRecruitment)
        ret->_pRecruitment = ret->_pRecruitment->clone();
    if (ret->_pVascular)
        ret->_pVascular = ret->_pVascular->clone();
    for(size_t i=0;i<NOUTCOMES;i++)
        if(ret->_pTTest[i])
            ret->_pTTest[i] = ret->_pTTest[i]->clone();
    return ret;
}


template<typename T>
inline void clearPtrList(T& container) {
    typename T::iterator end = container.end();
    for(typename T::iterator begin = container.begin(); begin!=end; begin++)
        delete *begin;
    container.clear();
}

GrSimulation::~GrSimulation()
{
    if(_pRecruitment) delete _pRecruitment;
    if (_pVascular) delete _pVascular;
    if(_pDiffusion)   delete _pDiffusion;
    for (int i = 0; i < NOUTCOMES; i++)
        if(_pTTest[i]) delete _pTTest[i];
    clearPtrList(_macList);
    clearPtrList(_tgamList);
    clearPtrList(_tregList);
    clearPtrList(_tcytList);
}

void GrSimulation::save(const char* fname) const
{
    std::ofstream ofs(fname);
    save(ofs);
}

void GrSimulation::save(std::ostream& out) const
{
    assert(out.good());
    std::locale loc = out.getloc();
    std::locale default_locale(out.getloc(), new boost::archive::codecvt_null<char>);
    std::locale _locale(default_locale, new boost::math::nonfinite_num_put<char>);
    out.imbue(_locale);
    boost::archive::xml_oarchive oa(out, boost::archive::no_codecvt);
    oa << boost::serialization::make_nvp("GR", *this);
    assert(out.good());
    out.imbue(loc);
}

void GrSimulation::load(const char* fname)
{
    std::ifstream ifs(fname);
    load(ifs);
}

void GrSimulation::load(std::istream& in)
{
    assert(in.good());
    std::locale loc = in.getloc();
    std::locale default_locale(in.getloc(), new boost::archive::codecvt_null<char>);
    std::locale _locale(default_locale, new boost::math::nonfinite_num_get<char>);
    in.imbue(_locale);
    boost::archive::xml_iarchive oa(in, boost::archive::no_codecvt);
    oa >> boost::serialization::make_nvp("GR", *this);
    assert(in.good());
    in.imbue(loc);
    timestepSync();
}

#if 0
void GrSimulation::deserialize(std::istream& in)
{

    timestepSync();


    assert(in.good());

    if (!Serialization::readHeader(in, GrSimulation::_ClassName))
    {
        exit(1);
    }

    std::string svnVersion;
    in >> svnVersion;

    // deserialize time
    in >> _time;

    // deserialize area threshold
    in >> _areaThreshold;
    in >> _areaThresholdCellDensity;

    // deserialize diffusion method
    int intVal;
    in >> intVal;
    setDiffusionMethod((DiffusionMethod) intVal);

    // deserialize recruitment method and object.
    int intRecrMtehod;
    in >> intRecrMtehod;
    RecruitmentMethod recrmethod = (RecruitmentMethod) intRecrMtehod;
    deserializeRecruitmentMethod(recrmethod, in);
    assert(in.good());
    _pRecruitment->deserialize(in); //LnODE recruitment will deserialize twice.
    assert(in.good());

    assert(in.good());
    in >> _tnfrDynamics;
    in >> _nfkbDynamics;
    in >> _il10rDynamics;
    in >> _tgammatransition;
    in >> (int&)(_odeSolver);
    in >> _tnfDepletionTimeStep;
    in >> _il10DepletionTimeStep;
    in >> _tcellRecruitmentBegun;
    in >> _numMolecularPerDiffusion;
    in >> _numDiffusionPerAgent;
    assert(in.good());

    // deserialize grid
    _grid.deserialize(in);

    // Deserialize the agent class.
    Agent::classDeserialize(in);

    // Since the cell lists are cleared here, the grid must not contain
    // any cell pointers at this point. If it does they are dangling pointers
    // that will be referenced during cell deserialization, when determining
    // if the grid compartment for the deserialized cell has room for the
    // deserialized cell. This will cause undefined (bad) behavior - seg faults, etc.

    // deserialize macs
    _macList.clear();
    in >> intVal;
    for (int i = 0; i < intVal; i++)
    {
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

        // create a dummy mac
        _macList.push_back(new Mac());
        Mac* pMac = _macList.back();

        // update attributes of dummy mac and add to grid
        pMac->deserialize(in);

        if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
        {
            cerr << "Error deserializing mac " << i << " at grid position " << pMac->getPosition() << ": invalid grid agent array index " << index << endl;
            exit(1);
        }

        bool addedAgent = _grid.getGrid().addAgent(pMac, index);
        if (!addedAgent)
        {
            cerr << "Error deserializing mac " << i << " at grid position " << pMac->getPosition() << ": addAgent failed."<< endl;
            exit(1);
        }
    }

    // deserialize tgam cells
    _tgamList.clear();
    in >> intVal;
    for (int i = 0; i < intVal; i++)
    {
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

        // create a dummy tgam cell
        _tgamList.push_back(new Tgam());
        Tgam* pTgam = _tgamList.back();

        if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
        {
            cerr << "Error deserializing tgam " << i << " at grid position " << pTgam->getPosition() << ": invalid grid agent array index " << index << endl;
            exit(1);
        }

        // update attributes of dummy tgam and add to grid
        pTgam->deserialize(in);
        bool addedAgent = _grid.getGrid().addAgent(pTgam, index);
        if (!addedAgent)
        {
            cerr << "Error deserializing tgam " << i << " at grid position " << pTgam->getPosition() << ": addAgent failed."<< endl;
            exit(1);
        }
    }

    // deserialize tcyt cells
    _tcytList.clear();
    in >> intVal;
    for (int i = 0; i < intVal; i++)
    {
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

        // create a dummy tcyt cell
        _tcytList.push_back(new Tcyt());
        Tcyt* pTcyt = _tcytList.back();

        if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
        {
            cerr << "Error deserializing tcyt " << i << " at grid position " << pTcyt->getPosition() << ": invalid grid agent array index " << index << endl;
            exit(1);
        }

        // update attributes of dummy tcyt and add to grid
        pTcyt->deserialize(in);
        bool addedAgent = _grid.getGrid().addAgent(pTcyt, index);
        if (!addedAgent)
        {
            cerr << "Error deserializing tcyt " << i << " at grid position " << pTcyt->getPosition() << ": addAgent failed."<< endl;
            exit(1);
        }

    }

    // deserialize treg cells
    _tregList.clear();
    in >> intVal;
    for (int i = 0; i < intVal; i++)
    {
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

        // create a dummy treg cell
        _tregList.push_back(new Treg());
        Treg* pTreg = _tregList.back();

        if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
        {
            cerr << "Error deserializing treg " << i << " at grid position " << pTreg->getPosition() << ": invalid grid agent array index " << index << endl;
            exit(1);
        }

        // update attributes of dummy mac and add to grid
        pTreg->deserialize(in);
        bool addedAgent = _grid.getGrid().addAgent(pTreg, index);
        if (!addedAgent)
        {
            cerr << "Error deserializing treg " << i << " at grid position " << pTreg->getPosition() << ":addAgent failed."<< endl;
            exit(1);
        }

    }

    // deserialize statistics
    _statsPrevious.serialize(in);
    _stats.serialize(in);

    // deserialize random number generator
    g_Rand.deserialize(in);

    if (!Serialization::readFooter(in, GrSimulation::_ClassName))
    {
        exit(1);
    }
}
#endif

static void initGrowthRateGrid(GrGrid& g) {
    for(int i=0;i<g.getRange().x;i++)
        for(int j=0;j<g.getRange().y;j++)
            g.growthRate(i,j) = g_Rand.getReal(_PARAM(Mtb_growthRateExtMtbMin),_PARAM(Mtb_growthRateExtMtbMax));
}

void GrSimulation::init(const boost::property_tree::ptree& pt, bool xml)
{
    // Before initializing anything check to see if the time step criteria are met
    timestepSync();
    
    // Used to initialize the valarray length for solving ODEs using RK4
    // Must be run before any agents are created so they have the correct
    // valarray length

    // initialize the sources
    _grid.initSources();
    const boost::property_tree::ptree& initpt = pt.get_child("GR.Init");
    PosVector initMacs, initExtMtb;

    std::pair<boost::property_tree::ptree::const_assoc_iterator, boost::property_tree::ptree::const_assoc_iterator>
            range = initpt.equal_range("Mac");
    while(range.first != range.second) {
        const boost::property_tree::ptree& p = (range.first++)->second;
        initMacs.push_back(Pos((xml ? p.get_child("<xmlattr>") : p), _grid.getCenter()));
    }

    range = initpt.equal_range("ExtMtb");
    while(range.first != range.second) {
        const boost::property_tree::ptree& p = (range.first++)->second;
        initExtMtb.push_back(Pos((xml ? p.get_child("<xmlattr>") : p), _grid.getCenter()));
    }

    const int maxMacAge = _PARAM(Mac_maxAge);

    // Place initial infected macrophages on the grid
    for (PosVector::const_iterator it = initMacs.begin(); it != initMacs.end(); it++)
    {
        Mac* pMac = createMac(it->x, it->y, g_Rand.getInt(-1, -maxMacAge), Mac::MAC_INFECTED, false, false);

        if (_nfkbDynamics)
        {
            // initialize NF-kB signaling from steady-state (12 hours equilibrium)
            for (int i = 0; i < 7200*_PARAM(_NFkBTimeCoefficient); ++i)
                pMac->solveNFkBODEsEquilibrium(6.00/_PARAM(_NFkBTimeCoefficient));
        }

        pMac->setIntMtb(1);
        ++_statsPrevious.getTotIntMtb();
        ++_stats.getTotIntMtb();

        //        // Place Treg in Mac Compartment for ODE Testing
        //        //DBG
        //        createTreg(it->x, it->y, g_Rand.getInt(-1, _PARAM(Tcell_maxAge), Treg::TREG_ACTIVE);
        //        //DBG

    }

    // Place initial extracellular bacteria on the grid
    for (PosVector::const_iterator it = initExtMtb.begin(); it != initExtMtb.end(); it++)
    {
        _grid.getGrid().extMTB(*it) += 1;
        ++_statsPrevious.getTotExtMtb();
        ++_stats.getTotExtMtb();
    }

    // Add a fixed number (PARAM_MAC_INIT_NUMBER) of resting
    // macrophages to the grid
    int count = _PARAM(Mac_initNumber);

    int nrGridCompartments = (_grid.getSize());
    if (count > nrGridCompartments)
    {
        std::cerr << "The number of initial resting macrophages to place on the grid, " << count << ", is > the number of grid compartments, " << nrGridCompartments << "." << std::endl;
        exit(1);
    }

    while (count > 0)
    {
        int row = g_Rand.getInt(_grid.getRange().x);
        int col = g_Rand.getInt(_grid.getRange().y);
        if (!_grid.getGrid().hasAgentType(MAC, Pos(row, col)))
        {
            Mac* newMac = createMac(row, col, g_Rand.getInt(-1, -maxMacAge), Mac::MAC_RESTING, false, false);
            if (_nfkbDynamics)
            {
                // initialize NF-kB signaling from steady-state (12 hours equilibrium)
                for (int i = 0; i < 7200*_PARAM(_NFkBTimeCoefficient); ++i)
                    newMac->solveNFkBODEsEquilibrium(6.00/_PARAM(_NFkBTimeCoefficient));
            }
            count--;
        }
    }
    
    initGrowthRateGrid(getGrid());

}

// Mark each cell within the molecular tracking radius to be tracked.
// Clear a cell's tracking status otherwise, so that the cells to track for a loaded
// state doesn't include those that were marked to be tracked when the state was saved.
// This allows defining the cells to be tracked based on the simulation state at the
// time the saved state was loaded, rather than based in part on the simulation state
// at time 0, when the simulation was initialized, and which would be carried forward to
// the saved and then loaded state.
void GrSimulation::initMolecularTracking(Scalar molecularTrackingRadius)
{
    if (molecularTrackingRadius == 0.0)
    {
        return;
    }

    const Pos center = _grid.getCenter();

    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    {
        if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
        {
            (*it)->setTrackMolecularDynamics(true);
        }
        else
        {
            (*it)->setTrackMolecularDynamics(false);
        }

    }

    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    {
        if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
        {
            (*it)->setTrackMolecularDynamics(true);
        }
        else
        {
            (*it)->setTrackMolecularDynamics(false);
        }
    }

    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    {
        if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
        {
            (*it)->setTrackMolecularDynamics(true);
        }
        else
        {
            (*it)->setTrackMolecularDynamics(false);
        }
    }

    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    {
        if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
        {
            (*it)->setTrackMolecularDynamics(true);
        }
        else
        {
            (*it)->setTrackMolecularDynamics(false);
        }
    }
}
void GrSimulation::initMolecularTracking(const std::vector<size_t>& ids) {

    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
        for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
            if(*i == (*it)->getID())
                (*it)->setTrackMolecularDynamics(true);
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
        for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
            if(*i == (*it)->getID())
                (*it)->setTrackMolecularDynamics(true);
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
        for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
            if(*i == (*it)->getID())
                (*it)->setTrackMolecularDynamics(true);
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
        for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
            if(*i == (*it)->getID())
                (*it)->setTrackMolecularDynamics(true);
}

void GrSimulation::solve()
{
    using namespace boost::posix_time;
    // here we perform a timestep, which is 10 minutes
    _time++;
    ptime startWallTime(microsec_clock::local_time());

#if 0
    //DBG
    // This is useful for determining whether or not 2 versions of the code are
    // using the same random number sequence, if they are using the same seed
    // and run on the same system.
    g_Rand.test(_time);
    //DBG
#endif

    // Ensure that both grids have the same state.
    _grid.updateNextGrid();

    // reset statistics
    _stats.reset();

    bool tnfDepletion = false;
    bool il10Depletion = false;

    if (_tnfDepletionTimeStep >= 0 && _time >= _tnfDepletionTimeStep)  {
        tnfDepletion = true;
        _tnfrDynamics = false;
    }

    if (_il10DepletionTimeStep >= 0 && _time >= _il10DepletionTimeStep)  {
        il10Depletion = true;
        _il10rDynamics = false;
    }

    // Get cell density information to use for adjusting diffusion (for Drug Diffusion)
    if (_PARAM(_DrugDynamics) && (_time >= _PARAM(_dosageStartTime)))  {
        GrGrid& g = _grid.getGrid();
        Pos p;
        Pos dim = g.getRange();
        for (p.x = 0; p.x < dim.x; p.x++)
        {
            for (p.y = 0; p.y < dim.y; p.y++)
            {
                g.nCells(p) = g.getCellDensity(p);
            }
        }
    }
    
    // Calculate Diffusion and Molecular events for a 10 min timestep
    // Agent time step is 600 s (10 min)
    for (int DiffStep = 0; DiffStep < _numDiffusionPerAgent; DiffStep++)
    {
        // Calculate signle diffusion time step
        _pDiffusion->diffuse(_grid, _time);

        // Consume drugs on the grid if they are turned on
        if (_PARAM(_DrugDynamics) && (_time >= _PARAM(_dosageStartTime)))  {
            consumeDrugs(_grid.getGrid(),_PARAM(_timestepDiffusion));
        }

        // Secrete from cells every diffusion time step (independent of molecular dynamics)
        secreteFromMacrophages(tnfDepletion, il10Depletion, _PARAM(_timestepDiffusion));
        secreteFromTcells(tnfDepletion, il10Depletion, _PARAM(_timestepDiffusion));
        secreteFromCaseations(_PARAM(_timestepDiffusion));

        // If molecular scale dynamics and an adaptive solver are chosen
        // Solve the ODEs using the adaptive Molecular Scale
        if ((_nfkbDynamics || _tnfrDynamics || _il10rDynamics) && _adaptive) {
            // Solve molecular scale dynamics with an adaptive solver
            solveMolecularScaleAdaptive(_PARAM(_timestepDiffusion));
            // If not using all molecular scale dynamics then adjust the CORE ABM consumption
            if (!_il10rDynamics || !_tnfrDynamics)  {
                for (int MolStep = 0; MolStep < _numMolecularPerDiffusion; MolStep++)  {
                    adjustFauxDegradation(_PARAM(_timestepMolecular));
                }
            }
        }
        // If molecular scale dynamics and a non-adaptive solver are chosen
        // Solve the ODEs using the non-adaptive solver loop
        else
            for (int MolStep = 0; MolStep < _numMolecularPerDiffusion; MolStep++)  {
                // Solve molecular scale dynamics with a non-adaptive solver
                if (_nfkbDynamics || _tnfrDynamics || _il10rDynamics)  {
                    solveMolecularScale(_PARAM(_timestepMolecular));
                }
                // If not using all molecular scale dynamics then adjust the CORE ABM consumption
                if (!_il10rDynamics || !_tnfrDynamics)  {
                    adjustFauxDegradation(_PARAM(_timestepMolecular));
                }
            }

        // Solve vascular dynamics for drugs if turned on
        if (_PARAM(_DrugDynamics) && (_time >= _PARAM(_dosageStartTime)))  {
            _pVascular->solveVascularSources(_grid.getGrid(), _PARAM(_timestepDiffusion), _time, DiffStep);
            _stats.getBloodConcINH() = _pVascular->getBloodConcentrationINH();
            _stats.getBloodConcRIF() = _pVascular->getBloodConcentrationRIF();
        }

    }

    // move macrophages
    moveMacrophages();

    // move T cells every 10 minutes
    moveTcells();

    // Shuffle the sources and recruit agents from vascular sources every 10 minutes
    getGrid().shuffleSources();
    _pRecruitment->recruit(*this, _time);

    // reset statistics
    _stats.resetAgentStats();

    // compute next state every 10 minutes
    computeNextStates();

    // update extracellular Mtb
    growExtMtb();

	updateDrugConcentrationStats();

    // update states and remove dead agents from lists and grid
    updateStates();

    // This must be after growExtMtb and updateStates, since updateStates updates the stats with intMtb count
    // and growExtMtb updates the stats with extMtb counts.
    // The mtb counts are set to zero at the start of function solve above.
    checkTCellRecruitmentStart();

    // shuffle order of cells every hour for randomization in movement
    if (_time % 6 == 0)  {
        shuffleCells();
    }

    // perform t-test
    updateT_Test();

    // Copy statistics, so we can use them on the next time step.
    _statsPrevious = _stats;

    ptime endWallTime(microsec_clock::local_time());
    time_duration dur = (endWallTime - startWallTime);
    _stats.setTimePerStep(dur.total_microseconds()*1.0e-6);
}

void GrSimulation::updateStates()
{
    // update states and remove dead macrophages from lists
    for (MacList::iterator it = _macList.begin(); it != _macList.end();)
    {
        Mac& mac = **it;
        mac.updateState();
        _stats.updateAgentStatistics(&mac);

        if (mac.isDead())
        {
            assert_res(_grid.getGrid().removeAgent(&mac));
            delete *it;
            it = _macList.erase(it);
        }
        else
        {
            it++;
        }
    }

    // update states and remove dead Tgams from lists
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end();)
    {
        Tgam& tgam = **it;
        tgam.updateState();
        _stats.updateAgentStatistics(&tgam);

        if (tgam.isDead())
        {
            assert_res(_grid.getGrid().removeAgent(&tgam));
            delete *it;
            it = _tgamList.erase(it);
        }
        else
        {
            it++;
        }
    }

    // update states and remove dead Tcyts from lists
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end();)
    {
        Tcyt& tcyt = **it;
        tcyt.updateState();
        _stats.updateAgentStatistics(&tcyt);

        if (tcyt.isDead())
        {
            assert_res(_grid.getGrid().removeAgent(&tcyt));
            delete *it;
            it = _tcytList.erase(it);
        }
        else
        {
            it++;
        }
    }

    // update states and remove dead Tregs from lists
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end();)
    {
        Treg& treg = **it;
        treg.updateState();
        _stats.updateAgentStatistics(&treg);

        if ((*it)->isDead())
        {
            assert_res(_grid.getGrid().removeAgent(&treg));
            delete *it;
            it = _tregList.erase(it);
        }
        else
        {
            it++;
        }
    }
}

void GrSimulation::computeNextStates()
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    {
        (*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
        (*it)->updateM1M2Ratio(); //Only for macs... should refactor this better...
    }
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    {
        (*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
    }
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    {
        (*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
    }
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    {
        (*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
    }
}

void GrSimulation::secreteFromMacrophages(bool tnfDepletion, bool il10Depletion, double mdt)
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    {
        (*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
    }
}

void GrSimulation::secreteFromTcells(bool tnfDepletion, bool il10Depletion, double mdt)
{
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    {
        (*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
    }
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    {
        (*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
    }
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    {
        (*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
    }
}

void GrSimulation::secreteFromCaseations(double mdt)
{
    // secrete chemokines from caseated compartments only if infection is not cleared
    Pos p;
    GrGrid& g = _grid.getGrid();
    const Pos& dim = g.getRange();
    for (p.x = 0; p.x < dim.x; p.x++)
    {
        for (p.y = 0; p.y < dim.y; p.y++)
        {
            if (g.isCaseated(p))
            {
                g.incCCL2(p, (0.25 * _PARAM(Mac_dCCL2) * mdt));
                g.incCCL5(p, (0.25 * _PARAM(Mac_dCCL5) * mdt));
                g.incCXCL9(p,  (0.25 * _PARAM(Mac_dCXCL9) * mdt));
                g.incmacAttractant(p, (_PARAM(_dAttractant) * mdt));
            }
        }
    }
}

void GrSimulation::solveMolecularScale(double dt)
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);   //t is not correct, should be getTime()*600+SEC_SINCE_LAST
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);

}
void GrSimulation::solveMolecularScaleAdaptive(double dt)
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
        (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);   //t is not correct, should be getTime()*600+SEC_SINCE_LAST
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
        if(!(*it)->getODEstatus())
            (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
        if(!(*it)->getODEstatus())
            (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
        if(!(*it)->getODEstatus())
            (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);

    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
        (*it)->setODEstatus(false);
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
        (*it)->setODEstatus(false);
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
        (*it)->setODEstatus(false);
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
        (*it)->setODEstatus(false);
}

void GrSimulation::moveMacrophages()
{
    const int timeResting = _PARAM(Mac_movementRest);
    const int timeInfected = _PARAM(Mac_movementInf);
    const int timeActive = _PARAM(Mac_movementAct);

    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    {
        // move resting macrophages
        if (_time % timeResting == 0 && (*it)->getState() == Mac::MAC_RESTING)
        {
            (*it)->move(_grid.getGrid());
        }
        // move active macrophages
        else if (_time % timeActive == 0 && (*it)->getState() == Mac::MAC_ACTIVE)
        {
            (*it)->move(_grid.getGrid());
        }
        // move infected macrophages
        else if (_time % timeInfected == 0 && (*it)->getState() == Mac::MAC_INFECTED)
        {
            (*it)->move(_grid.getGrid());
        }
    }
}

void GrSimulation::moveTcells()
{

    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    {
        (*it)->move(_grid.getGrid());
    }
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    {
        (*it)->move(_grid.getGrid());
    }
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    {
        (*it)->move(_grid.getGrid());
    }
}

void GrSimulation::growExtMtb()
{
    const double upperBound = _PARAM(Mtb_growthExtMtbBound);

    GrGrid& g = _grid.getGrid();
    const Pos& dim = g.getRange();
    Pos p;
    for (p.x = 0; p.x < dim.x; p.x++)
    {
        for (p.y = 0; p.y < dim.y; p.y++)
        {
            double growthRate = _PARAM(Mtb_growthRateExtMtb) - 1;
            if(_PARAM(_RandomizeGrowthRate))
                growthRate = getGrid().growthRate(p) - 1;

            double adj_kill_rate = 0.0;
            if (_PARAM(_DrugDynamics) && (_time >= _PARAM(_dosageStartTime)))
            {
                double INH_Emax = _PARAM(_extraActivityINH);
                double INH_C50 = _PARAM(_C50_INH); // milligram/L
                double INH_C = g.INH(p);        // get concentration from grid in mg/L
                adj_kill_rate = INH_Emax * ((INH_C) / ( INH_C50 + INH_C ));
            }

            Scalar& extMtb = g.extMTB(p);

            if (g.isCaseated(p))  {
                // Bacteria don't gp.x in caseated compartments

                // Ext Mtb in caseation has to die somehow or it will increase with the nrCaseated compartments

                double killProb = (extMtb/upperBound);
                double dExtMtb = 0.0;

                if (g_Rand.getReal() <= killProb)
                    dExtMtb = _PARAM(Mtb_deathRateExtMtbCaseated) * extMtb * (killProb);
                if (dExtMtb > extMtb)
                    dExtMtb = extMtb;

                extMtb += (-1.0 * dExtMtb);

                // Bacteria that survive can die from antibiotics in caseated compartments given that the local concentration is high enough. We can also
                // later adjust kill rate for caseum (i.e. if local conditions such as oxygen and pH affect drug activity).
                dExtMtb = -adj_kill_rate * extMtb;
                assert(-dExtMtb<=extMtb);
                extMtb += dExtMtb;

                ++_stats.getNrCaseated();
                _stats.getTotNonRepExtMtb()+=(extMtb);
                _stats.getTotExtMtb()+=(extMtb);
            }
            else  {
                // Function to evaluate whether the extMtb is 'trapped' or not... (very basic)
                // Currently only identifies cells whose Moore neighborhood (minus its own compartment)
                // is completley caseated
                int caseationCount;

                if (extMtb > 0.0)
                    caseationCount = g.isTrapped(p);
                else
                    caseationCount = 0;

                // Scale growth rate based on local caseation
                // This mimicks the hypoxic environment in granulomas which causes Mtb to decrease its growth rate
                // This should also prevent the extMtb that are not caught by the trapped function to stop growing as fast
                // Ext Mtb with caseationCount >= 5 DO NOT GROW

                double dExtMtb = 0.0;
                if (caseationCount >= 2 && caseationCount < 5)
                    dExtMtb = (growthRate/caseationCount) * extMtb * (1.0 - extMtb / upperBound) - adj_kill_rate * extMtb;
                else if (caseationCount >= 5)
                    dExtMtb = 0.0 - adj_kill_rate * extMtb;
                else
                    dExtMtb = growthRate * extMtb * (1.0 - extMtb / upperBound) - adj_kill_rate * extMtb;
                assert(-dExtMtb<=extMtb);

                extMtb += dExtMtb;
                _stats.getTotExtMtb() += (extMtb);
            }

            _stats.getTotMacAttractant() += (g.macAttractant(p));
            _stats.getTotTNF() += (g.TNF(p));
            _stats.getTotIL10() += (g.il10(p));
            _stats.getTotCCL2() += (g.CCL2(p));
            _stats.getTotCCL5() += (g.CCL5(p));
            _stats.getTotCXCL9() += (g.CXCL9(p));

            if (g.TNF(p) >= _areaThreshold)
                ++_stats.getAreaTNF();

            if (g.getCellDensity(p) >= _areaThresholdCellDensity)
            {
                ++_stats.getAreaCellDensity();
            }

            if (_PARAM(_DrugDynamics) && (_time >= _PARAM(_dosageStartTime)))
            {
                _stats.getTotINH() += (g.INH(p));
                _stats.getTotRIF() += (g.RIF(p));
                if (g.getCellDensity(p) >= _areaThresholdCellDensity)
                {
                    _stats.getTotINHGran() += g.INH(p);
                    _stats.getTotRIFGran() += g.RIF(p);
                }
                else
                {
                    _stats.getTotINHNorm() += g.INH(p);
                    _stats.getTotRIFNorm() += g.RIF(p);
                }
            }
        }
    }
}

// Update the area under the curve of average drug concentration for the entire grid vs. time,
// for an interval of consecutive time steps, specified by the _drugConcentrationStatInterval
// parameter file parameter.
//
// Stat DrugConcentration contains the total concentration of the drug on the grid for the last 
// _drugConcentrationStatInterval time steps. However, they are not in order by time step since
// that would require extra moving of vector elements (deleting the first element in the vector
// by moving up all the subsequent elements and then entering the total grid drug concentration
// for the current time step to the end of the vector). Instead we use timestep % interval to
// choose the index into the vector for the value for the oldest timestep and replace that with 
// the newest version. For example, if we are tracking for an interval of 3 time steps the drug
// concentration vector would be the following to time steps 1 to 7:
// ts	drugconc	vector
// 0	0			[ 0.0, 0.0, 0.0] 0 elements defined
// 1    3.2			[ 0.0, 3.2, 0.0] 1 element  defined
// 2	2.1			[ 0.0, 3.2, 2.1] 2 elements defined
// 3	4.3			[ 4.3, 3.2, 2.1] 3 elements defined
// 4	1.7			[ 4.3, 1.7, 2.1] 3 elements defined
// 5	8.4			[ 4.3, 1.7, 8.4] 3 elements defined
//
// The number of elements defined will be less than the _drugConcentrationStatInterval
// parameter, so the drugConcentrationArea value will not be correct until 
// the simulation time step >= _drugConcentrationStatInterval. This is not
// considered a problem since the drugConcentrationArea value is needed for the end of
// a simulation, and the _drugConcentrationStatInterval parameter will be very small
// compared to the length of a simulation.
//
// For now we do this specifically for INH. We may add support for other drugs later.

void GrSimulation::updateDrugConcentrationStats()
{
	_stats.setDrugConcentrationINH(_time % _PARAM(_drugConcentrationStatInterval), _stats.getTotINH());

	std::vector<Scalar> drugConc =_stats.getDrugConcentrationINH();
	Scalar drugConcentrationArea = 0.0;

	for (size_t i = 0; i < drugConc.size(); ++i)
	{
		drugConcentrationArea += drugConc[i];
	}

	// Convert from total grid concentration to grid concentration average per
	// grid compartment.
	_stats.setDrugConcentrationAreaINH(drugConcentrationArea / _grid.getSize());
}

void GrSimulation::adjustFauxDegradation(double dt)
{
    
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    {
        (*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
    }
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    {
        (*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
    }
    for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    {
        (*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
    }
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    {
        (*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
    }
}

void GrSimulation::consumeDrugs(GrGrid grid, double dt)
{
    Pos pos;
    for (pos.x = 0; pos.x < grid.getRange().x; pos.x++)
    {
          for (pos.y = 0; pos.y < grid.getRange().y; pos.y++)
          {
              // simulate the consumption of drugs by extracellular bacteria
              double dINH;
              Scalar INH = grid.INH(pos);   //  mg/L
              double mtb = grid.extMTB(pos);

              if (grid.extMTB(pos)!=0)
              {
                  dINH = -_PARAM(_ConsumptionRate_INH) * (INH / (INH + _PARAM(_ConsumptionRate_INH))) * dt * mtb;  // M
                  dINH = dINH * MW_INH;    // mg/L
                  grid.incINH(pos, dINH);
              }
          }
    }
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    {
        if (Mac::isMac(*it,Mac::MAC_INFECTED)||Mac::isMac(*it,Mac::MAC_CINFECTED))
        {
            // simulate the consumption of drugs by intracellular bacteria - for now from the same pool as extracellular
            double dINH;
            pos.x = (*it)->Agent::getRow();
            pos.y = (*it)->Agent::getCol();
            Scalar INH = grid.INH(pos);   //  mg/L
            double mtb = (*it)->Mac::getIntMtb();

            dINH = -_PARAM(_ConsumptionRate_INH) * (INH / (INH + _PARAM(_ConsumptionRate_INH))) * dt * mtb;  // M
            dINH = dINH * MW_INH;    // mg/L
            assert(abs(dINH)<=INH);
            grid.incINH(pos, dINH);
            double INHafter = grid.INH(pos);
            assert(INHafter>=0);
        }
    }
}


void GrSimulation::updateT_Test()
{
    for (int i = 0; i < NOUTCOMES; i++)
    {
        if (_pTTest[i])
        {
            _pTTest[i]->update(_time, i, _stats);
        }
    }
}

void GrSimulation::performT_Test()
{
    for (int i = 0; i < NOUTCOMES; i++)
    {
        if (_pTTest[i])
        {
            _pTTest[i]->perform(i, _stats);
        }
    }
}

void GrSimulation::setDiffusionMethod(DiffusionMethod method)
{
    if (_pDiffusion) {
        if (_pDiffusion->getMethod() == method)
            return;
        delete _pDiffusion;
    }

    switch (method)
    {
    case DIFF_SOR_CORRECT:
        _pDiffusion = new GrDiffusionBTCS(); break;
    case DIFF_SOR_WRONG:
        _pDiffusion = new GrDiffusionWrongBTCS(); break;
    case DIFF_REC_EQ:
        std::cerr<<"*** WARNING: DIFF_REC_EQ disabled, defaulting to DIFF_REC_EQ_SWAP"<<std::endl;
    case DIFF_REC_EQ_SWAP:
        _pDiffusion = new GrDiffusionFTCS_Swap(); break;
    case DIFF_ADE_SWAP:
        if(_PARAM(_timestepDiffusion) > 30)
            std::cerr<<("*** WARNING: This diffusion method is inaccurate for timesteps greater than 30 seconds")<<std::endl;
        _pDiffusion = new GrDiffusionADE_Swap(); break;
    default:
        assert(!"Invalid diffusion method"); break;
    }
    if(method != DIFF_ADE_SWAP && _PARAM(_timestepDiffusion) > 12)
        std::cerr<<("*** WARNING: This diffusion method is unstable for timesteps greater than 12 seconds")<<std::endl;
}

RecruitmentMethod GrSimulation::getRecruitmentMethod()
{
    return _pRecruitment->getMethod();
}

void GrSimulation::setRecruitmentMethod(RecruitmentMethod method)
{
    if (_pRecruitment)
    {
        if (method == _pRecruitment->getMethod())
        {
            // Already have the correct type of recruitment object.
            // Don't want to replace it since it may have been deserialized from a saved simulation state.
            return;
        }

        delete _pRecruitment;
    }

    switch (method)
    {
    case RECR_PROB:
        _pRecruitment = new RecruitmentProb();
        break;

    case RECR_LN_ODE_PURE:
        _pRecruitment = new RecruitmentLnODEPure();
        break;

    default:
        std::cerr << "Invalid recruitment method: " << method << std::endl;
        exit(1);
    }
}

void GrSimulation::setOutcomeMethod(int index, OutcomeMethod method, double alpha,
                                    int testPeriod, int samplePeriod)
{
    assert(0 <= index && index < NOUTCOMES);

    if (_pTTest[index] && _pTTest[index]->getMethod() == method)
    {
        _pTTest[index]->setAlpha(alpha);
        _pTTest[index]->setTestPeriod(testPeriod);
        _pTTest[index]->setSamplePeriod(samplePeriod);
    }
    else
    {
        delete _pTTest[index];

        switch (method)
        {
        case OUTCOME_AREA:
            _pTTest[index] = new AreaTest(alpha, testPeriod, samplePeriod);
            break;
        case OUTCOME_MTB:
            _pTTest[index] = new MtbTest(alpha, testPeriod, samplePeriod);
            break;
        case OUTCOME_NONE:
            _pTTest[index] = NULL;
            break;
        default:
            assert(!"Invalid outcome method");
        }
    }
}

void GrSimulation::shuffleCells()
{
#if 1
    std::random_shuffle(_macList.begin(), _macList.end(), g_Rand);
    std::random_shuffle(_tgamList.begin(), _tgamList.end(), g_Rand);
    std::random_shuffle(_tcytList.begin(), _tcytList.end(), g_Rand);
    std::random_shuffle(_tregList.begin(), _tregList.end(), g_Rand);
#else
    //TODO: replace with random_shuffle
    MacList::iterator mac1 = _macList.begin();
    MacList::iterator mac2 = _macList.end();
    TgamList::iterator tgam1 = _tgamList.begin();
    TgamList::iterator tgam2 = _tgamList.end();
    TcytList::iterator tcyt1 = _tcytList.begin();
    TcytList::iterator tcyt2 = _tcytList.end();
    TregList::iterator treg1 = _tregList.begin();
    TregList::iterator treg2 = _tregList.end();
    MacList::size_type m = 1;
    MacList::size_type n = 1;

    // Don't shuffle an empty list.
    // The list size is unsigned so the shuffle loop limit, _list size - 1, is also unsigned.
    // When list size is 0, the limit is a large positive number instead of -1.

    if (_macList.size() > 0)
    {
        while (m + n < _macList.size() - 1)
        {
            if (g_Rand.getReal() < 0.5)
            {
                mac1++;
                m++;
            }
            if (g_Rand.getReal() < 0.5)
            {
                mac2--;
                n++;
            }
            std::swap(mac1, mac2);
        }
    }
    if (_tgamList.size() > 0)
    {
        m = 1;
        n = 1;
        while (m + n < _tgamList.size() - 1)
        {
            if (g_Rand.getReal() < 0.5)
            {
                tgam1++;
                m++;
            }
            if (g_Rand.getReal() < 0.5)
            {
                tgam2--;
                n++;
            }
            std::swap(tgam1, tgam2);
        }
    }

    if (_tcytList.size() > 0)
    {
        m = 1;
        n = 1;
        while (m + n < _tcytList.size() - 1)
        {
            if (g_Rand.getReal() < 0.5)
            {
                tcyt1++;
                m++;
            }
            if (g_Rand.getReal() < 0.5)
            {
                tcyt2--;
                n++;
            }
            std::swap(tcyt1, tcyt2);
        }
    }

    if (_tregList.size() > 0)
    {
        m = 1;
        n = 1;
        while (m + n < _tregList.size() - 1)
        {
            if (g_Rand.getReal() < 0.5)
            {
                treg1++;
                m++;
            }
            if (g_Rand.getReal() < 0.5)
            {
                treg2--;
                n++;
            }
            std::swap(treg1, treg2);
        }
    }
#endif
}
