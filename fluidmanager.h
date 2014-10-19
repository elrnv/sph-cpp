#ifndef FLUIDMANAGER_H
#define FLUIDMANAGER_H

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

class FluidManager
{
public:
  FluidManager();
  ~FluidManager();

  template<int FT> // templated "typedef"
  using FluidVecT = std::vector< FluidT<FT> >;
  
#define FLUIDVEC_GETTER(FT) \
  FluidVecT<FT> &get_fluids<FT>() { return m_fluids_FT; }

  FLUIDVEC_GETTER(DEFAULT);
  FLUIDVEC_GETTER(MCG03);
  FLUIDVEC_GETTER(BT07);
  FLUIDVEC_GETTER(ICS13);

private:
#define FLUIDVEC_MEMEBER(FT) \
  FluidVecT<FT> m_fluids_FT;

  FLUIDVEC_MEMBER(DEFAULT);
  FLUIDVEC_MEMBER(MCG03);
  FLUIDVEC_MEMBER(BT07);
  FLUIDVEC_MEMBER(ICS13);
};

#endif // FLUIDMANAGER_H
