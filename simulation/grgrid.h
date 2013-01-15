/*
 * grgrid.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef GRGRID_H
#define GRGRID_H

#include "gr.h"
#include "lungparams.h"
#include "pos.h"


// Sets up the world geometry
struct Indexer2D
{
  static inline size_t ind(const Pos& dim, const Pos& idx)
  {
    assert(idx.x >= 0 && idx.y >= 0);
    return ind(dim, idx.x, idx.y);
  }
  static inline size_t ind(const Pos& dim, size_t x, size_t y)
  {
    assert(x < (size_t)dim.x && y < (size_t)dim.y);
    return x * dim.y + y;
  }
  static inline size_t wrap(int d, int x)
  {
    assert((x + d) > 0);
    return (x + d) % d;
  }
  static inline size_t padInd(const Pos& dim, const Pos& idx)
  {
    return padInd(dim, idx.x, idx.y);
  }
  static inline size_t padInd(const Pos& dim, int x, int y)
  {
    assert(x > -2 && x < dim.x+2 && y > -2 && y < dim.y+2);
    return (x+1) * (dim.y+2) + (y+1);

  }
};
typedef Indexer2D Indexer;

#define GRIDS_DEFS  \
  GRID       (int, nKillings)         \
  GRID       (int, trappedCaseation)  \
  GRID       (int, nRecruitments)     \
  GRID       (int, nRecruitmentsMac)  \
  GRID       (int, nRecruitmentsTgam) \
  GRID       (int, nRecruitmentsTcyt) \
  GRID       (int, nRecruitmentsTreg) \
  GRID       (int, nSecretions)       \
  GRID       (int, nCells)            \
  PADDED_GRID(Scalar, TNF)            \
  PADDED_GRID(Scalar, macAttractant)  \
  PADDED_GRID(Scalar, INH)            \
  PADDED_GRID(Scalar, RIF)            \
  PADDED_GRID(Scalar, CCL2)           \
  PADDED_GRID(Scalar, CCL5)           \
  PADDED_GRID(Scalar, CXCL9)          \
  PADDED_GRID(Scalar, shedTNFR2)      \
  PADDED_GRID(Scalar, il10)           \
  GRID       (Scalar, extMTB)         \
  GRID       (Scalar, growthRate)     \
 

class GrGrid
{
public:
  static const unsigned MAX_AGENTS_PER_CELL = 2;
private:
  // Make this private so no one can accidently access it
  // (just for serialization)
  friend class boost::serialization::access;  
  GrGrid(){} //For serialization

  /*
   * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
   */
  Pos _dim;
  int _nCaseation;
  std::vector<Pos> _sources;
  std::vector<Agent*> _agents;

#define GRID(type, name) \
    std::vector<type> _##name;

#define PADDED_GRID(type, name) \
    std::vector<type> _##name;	\
    std::vector<type> _u_##name;	\
    std::vector<type> _v_##name;	\
 
  GRIDS_DEFS
#undef GRID
#undef PADDED_GRID

public:
  #define GRID(type, name) IDX_##name,
  #define PADDED_GRID(type, name) IDX_##name,
  enum GRID_IDX {
    __INIT_GRID_IDX__ = -1, //Don't use, just to ensure standard compliance
   GRIDS_DEFS
   IDX_NGRIDS
  };
#undef GRID
#undef PADDED_GRID

  GrGrid(const Pos& dim);
  ~GrGrid();

  bool inRange(const Pos& pos) const;
  bool inRange(int row, int col) const;
  const Pos& getRange() const
  {
    return _dim;
  }
  int getSize() const
  {
    return _dim.x*_dim.y;
  }
  const Pos getCenter() const
  {
    return _dim / 2;
  }

  const PosVector& getSources() const;
  void initSources();
  void shuffleSources();
  /**
  * @copydoc GrSimulation::serialize
  */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
  bool isSource(const Pos& p) const;

  int getOccupiedNeighborCount(int row, int col) const;
  int getOccupiedNeighborCount(const Pos& p) const;
  Scalar getCellDensity(int row, int col) const;
  Scalar getCellDensity(const Pos& p) const;
  bool isOccupied(int row, int col) const;
  bool isOccupied(const Pos& p) const;
  bool isCaseated(int row, int col) const;
  bool isCaseated(const Pos& p) const;
  int getNumberOfAgents(int row, int col) const;
  int getNumberOfAgents(const Pos& p) const;
  int hasTcell(const Pos& p) const;
  int hasAgentType(AgentType TYPE, const Pos& p) const;
  int hasAgentType(AgentType TYPE, int s, const Pos& p) const;
  bool addAgent(Agent* a, const Pos& p);
  bool addAgent(Agent* a);
  bool addAgent(Agent* a, int row, int col);
  bool addAgent(Agent* a, const Pos& p, int index);
  bool addAgent(Agent* a, int index);
  bool addAgent(Agent* a, int row, int col, int index);
  int agentIndex (const Agent* pAgent) const;
  bool removeAgent(const Agent* a);
  inline int mod_row(int row) const
  {
    return (row + _dim.x) % _dim.x;
  }
  inline int mod_col(int col) const
  {
    return (col + _dim.y) % _dim.y;
  }

  /*** Accessors ***/
#define GRID(type, name) \
  const type& name (const Pos& p) const;  \
  const type& name (int x, int y) const;  \
  const type* name () const;  \
        type& name (const Pos& p);  \
        type& name (int x, int y);  \
        type* name ();
#define PADDED_GRID(type, name) \
  const type& name (const Pos& p) const;  \
  const type& name (int x, int y) const;  \
        void set##name (const Pos& p, type v);  \
        void set##name (int x, int y, type v);  \
        void inc##name (const Pos& p, type v);  \
        void inc##name (int x, int y, type v);  \
  const type* name () const;  \
        type* name ();	\
  const type& u_##name (const Pos& p) const;  \
  const type& u_##name (int x, int y) const;  \
  const type* u_##name () const;  \
        type* u_##name ();	\
  const type& v_##name (const Pos& p) const;  \
  const type& v_##name (int x, int y) const;  \
  const type* v_##name () const;  \
        type* v_##name ();\
 
  GRIDS_DEFS

#undef GRID
#undef PADDED_GRID

  const Agent*  agent         (const Pos& p, size_t a) const;
  Agent*&       agent         (const Pos& p, size_t a);
  const Agent*  agent         (int x, int y, size_t a) const;
  Agent*&       agent         (int x, int y, size_t a);
  bool incKillings(const Pos& p);
  int isTrapped(const Pos& p);

// --- Visitor Methods ---

  /*! \brief visit each grid
  * \tparam Visitor Visitor type that models the Visitor concept for each stat
  * \param v Visitor that visits each grid
  */
  template<typename Visitor>
  void visit(Visitor& v)
  {
#define GRID(type, name)                v.visit(#name, _ ## name, "");
#define PADDED_GRID(type, name)         v.visit(#name, _ ## name, "");
    GRIDS_DEFS
#undef GRID
#undef PADDED_GRID
  }

  /*! \brief visit each grid
  * \tparam Visitor Visitor type that models the Visitor concept for each stat
  * \param v Visitor that visits each grid
  */
  template<typename Visitor>
  void visit(Visitor& v) const
  {
#define GRID(type, name)                v.visit(#name, _ ## name, "");
#define PADDED_GRID(type, name)         v.visit(#name, _ ## name, "");
    GRIDS_DEFS
#undef GRID
#undef PADDED_GRID
  }

  /**
   * @brief Given an templated index this sets val to the result of calling
   * the named grid's indexer.
   * @note This is probably not the best way to do this, but we're relying on
   * compiler optimization to remove the switch case entirely.  The reason for
   * implementing this is so we can index each grid, padded or not, in the same
   * manner, knowing only the index of the grid we want.  Only should be used
   * in the gui.
   * @details example:
   * getIndexedValue(IDX_TNF,0,0,val) is the same as
   * val = TNF(0,0);
   *
   * @tparam T type to cast to.  Must be copyable and assignable
   * @param i Index of the grid to get the value of
   * @param x Row
   * @param y Column
   * @param val Value to store result in
   */
  template<typename T>
  void getIndexedValue(GRID_IDX i, int x, int y, T& val) const {
    switch(i) {
#define GRID(type, name)        case IDX_##name : val = T(name(x,y)); break;
#define PADDED_GRID(type, name) case IDX_##name : val = T(name(x,y)); break;
        GRIDS_DEFS
#undef GRID
#undef PADDED_GRID
      default: assert(!"Invalid index"); break;
    }
  }
};

inline bool GrGrid::isSource(const Pos& p) const
{
  return std::find(_sources.begin(), _sources.end(), p) != _sources.end();
}

inline bool GrGrid::inRange(const Pos& pos) const
{
  return inRange(pos.x, pos.y);
}
inline bool GrGrid::inRange(int row, int col) const
{
  return row < _dim.x && col < _dim.y;
}

inline int GrGrid::getNumberOfAgents(int row, int col) const
{
  int sum = 0;
  for(unsigned i=0; i<MAX_AGENTS_PER_CELL; i++)
    sum += agent(Pos(row, col), i) != NULL ? 1 : 0;
  return sum;
}
inline int GrGrid::getNumberOfAgents(const Pos& p) const
{
  return getNumberOfAgents(GETROW(p), GETCOL(p));
}

inline bool GrGrid::isOccupied(int row, int col) const
{
  assert(row < _dim.x && col < _dim.y);
  return isCaseated(row, col) || (getNumberOfAgents(row, col) > 0) || (extMTB(Pos(row, col)) > 0.0);
}
inline bool GrGrid::isOccupied(const Pos& p) const
{
  return isOccupied(GETROW(p), GETCOL(p));
}

inline bool GrGrid::isCaseated(int row, int col) const
{
  return ((nKillings(Pos(row, col)) >= _nCaseation) || trappedCaseation(Pos(row,col)) >= 1);
}
inline bool GrGrid::isCaseated(const Pos& p) const
{
  return isCaseated(GETROW(p), GETCOL(p));
}

inline const std::vector<Pos>& GrGrid::getSources() const
{
  return _sources;
}

#define GRID(type, name) \
  inline const type& GrGrid:: name (const Pos& p) const  \
    { return _##name [Indexer::ind(_dim, p)]; } \
  inline const type& GrGrid:: name (int x, int y) const  \
    { return _##name [Indexer::ind(_dim, x, y)]; }  \
  inline const type* GrGrid:: name () const  \
      { return _##name.data(); }  \
  inline type& GrGrid:: name (const Pos& p) \
    { return _##name [Indexer::ind(_dim, p)]; } \
  inline type& GrGrid:: name (int x, int y)  \
    { return _##name [Indexer::ind(_dim, x, y)]; } \
  inline type* GrGrid:: name ()  \
      { return _##name.data(); }
#define PADDED_GRID(type, name) \
  inline const type& GrGrid:: name (const Pos& p) const  \
    { return _##name [Indexer::padInd(_dim, p)]; } \
  inline const type& GrGrid:: name (int x, int y) const  \
    { return _##name [Indexer::padInd(_dim, x, y)]; } \
  inline const type* GrGrid:: name () const  \
    { return _##name.data(); } \
  inline       type* GrGrid:: name () \
    { return _##name.data(); }  \
  inline const type& GrGrid:: u_##name (const Pos& p) const  \
    { return _u_##name [Indexer::padInd(_dim, p)]; } \
  inline const type& GrGrid:: u_##name (int x, int y) const  \
    { return _u_##name [Indexer::padInd(_dim, x, y)]; } \
  inline const type* GrGrid::u_##name () const  \
    { return _u_##name.data(); } \
  inline       type* GrGrid::u_##name () \
    { return _u_##name.data(); }  \
  inline const type& GrGrid:: v_##name (const Pos& p) const  \
    { return _v_##name [Indexer::padInd(_dim, p)]; } \
  inline const type& GrGrid:: v_##name (int x, int y) const  \
    { return _v_##name [Indexer::padInd(_dim, x, y)]; } \
  inline const type* GrGrid::v_##name () const  \
    { return _v_##name.data(); } \
  inline       type* GrGrid::v_##name () \
    { return _v_##name.data(); }  \
  inline void GrGrid::set##name (const Pos& p, type v) \
    { set##name(p.x, p.y, v); } \
  inline void GrGrid::set##name (int x, int y, type v)  \
    { _u_##name[Indexer::padInd(_dim, x, y)] += v - _##name[Indexer::padInd(_dim, x, y)];					\
      _v_##name[Indexer::padInd(_dim, x, y)] += v - _##name[Indexer::padInd(_dim, x, y)];					\
      _##name [Indexer::padInd(_dim, x, y)] = v; \
    }	\
  inline void GrGrid::inc##name (const Pos& p, type v) \
    { inc##name(p.x, p.y, v); } \
  inline void GrGrid::inc##name (int x, int y, type v)  \
    { _u_##name[Indexer::padInd(_dim, x, y)] += v;					\
      _v_##name[Indexer::padInd(_dim, x, y)] += v;					\
      _##name [Indexer::padInd(_dim, x, y)] += v; \
    }

GRIDS_DEFS

#undef GRID
#undef PADDED_GRID

inline const Agent*  GrGrid::agent         (const Pos& p, size_t a) const
{
  return _agents[Indexer::ind(_dim, p) * MAX_AGENTS_PER_CELL + a];
}
inline Agent*& GrGrid::agent(const Pos& p, size_t a)
{
  return _agents[Indexer::ind(_dim, p) * MAX_AGENTS_PER_CELL + a];
}
inline const Agent*  GrGrid::agent         (int x, int y, size_t a) const
{
  return _agents[Indexer::ind(_dim, x, y) * MAX_AGENTS_PER_CELL + a];
}
inline Agent*& GrGrid::agent(int x, int y, size_t a)
{
  return _agents[Indexer::ind(_dim, x, y) * MAX_AGENTS_PER_CELL + a];
}

template<class Archive>
void GrGrid::serialize(Archive& ar, const unsigned int /*version*/)
{
  ar & BOOST_SERIALIZATION_NVP(_dim);
  #define GRID(type, name) ar & BOOST_SERIALIZATION_NVP(_##name);
  #define PADDED_GRID(type, name) \
    ar & BOOST_SERIALIZATION_NVP(_##name);  \
    ar & BOOST_SERIALIZATION_NVP(_u_##name);  \
    ar & BOOST_SERIALIZATION_NVP(_v_##name);
   GRIDS_DEFS
  #undef GRID
  #undef PADDED_GRID
  ar & BOOST_SERIALIZATION_NVP(_agents);
  ar & BOOST_SERIALIZATION_NVP(_nCaseation);
  ar & BOOST_SERIALIZATION_NVP(_sources);
}

inline std::ostream& operator<<(std::ostream& s, GrGrid::GRID_IDX e) {
#define GRID(type, name)        case GrGrid::IDX_##name : s<<#name; break;
#define PADDED_GRID(type, name) case GrGrid::IDX_##name : s<<#name; break;
  switch(e) {
    GRIDS_DEFS
  default:
    assert(!"Invalid grid index."); break;
  }
#undef GRID
#undef PADDED_GRID
  return s;
}

#endif /* GRID_H */
