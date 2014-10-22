#ifndef DYNAMICSMANAGER_H
#define DYNAMICSMANAGER_H

template<int FT>
class FluidT<FT>;

enum __attribute__ ((__packed__)) FluidType
{
  DEFAULT = 0,
  MCG03 = 1,
  BT07 = 2,
  ICS13 = 3,
  NUMTYPES = 4;
};

// When adding a new type of fluid, add a line to the macro below.
#define FOREACH_FT_EVAL_MACRO(FUNC) \
  FUNC(DEFAULT); \
  FUNC(MCG03); \
  FUNC(BT07); \
  FUNC(ICS13);

class DynamicsManager
{
public:
  DynamicsManager();
  ~DynamicsManager();

  // templated "typedef"
  template<int FT>
  using FluidVecT = std::vector< FluidT<FT> >;

  typedef std::vector< Boundary > BoundaryVec;
  
// The following macros generate the fluid vector members and their respective 
// getters. 
#define FLUIDVEC_GETTER(FT) \
  FluidVecT<FT> &get_fluids<FT>() { return m_fluids_FT; }

  FOREACH_FT_EVAL_MACRO(FLUIDVEC_GETTER)

  BoundaryVec &get_bounds() { return m_boundaries; }

private:
#define FLUIDVEC_MEMEBER(FT) \
  FluidVecT<FT> m_fluids_FT;

  FOREACH_FT_EVAL_MACRO(FLUIDVEC_MEMBER)

  BoundaryVec m_boundaries;
};

#endif // DYNAMICSMANAGER_H
