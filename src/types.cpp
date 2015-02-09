#include "types.h"
AlignedBox3f UnitBox(Vector3f(-1.f,-1.f,-1.f), Vector3f(1.f,1.f,1.f));

#define SPH_TYPE_TO_STRING(unused,data,elem) BOOST_PP_STRINGIZE(elem) ,
const char * SPHParticleTypeString[] = {
  BOOST_PP_SEQ_FOR_EACH(SPH_TYPE_TO_STRING,_,SPH_TYPES)
};
#undef SPH_TYPE_TO_STRING
