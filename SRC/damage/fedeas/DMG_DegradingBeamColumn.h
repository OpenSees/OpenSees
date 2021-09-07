/* 
 * libdmg
 * TODO: Copyright, authors
 */

// Claudio Perez

#include <string.h>
#include <map>

#include <tcl.h>
#include <g30/matrix.h>

template<int N>
struct BaseState{
/* This struct defines the interface between
 * the public beam column wrapper and it's
 * external implementation.
 */
  g23::vector<N>   q;
  g23::vector<N>   se;
  g23::vector<N>   vp;
  g23::vector<N>   vpP;
  g23::matrix<N,N> k;
  g23::matrix<N,N> ke;
  g23::matrix<N,N> kt;
};


template<typename BaseElem>
class DMG_DegradingBeamColumn:
  public BaseElem
{
public:

  typedef int (*degrade_f)(BaseState<3>* past, BaseState<3>* pres);
  
  int setLibDmgWrapper(Tcl_Interp *interp, char *tag){
    auto newDamageWrap = (degrade_f* (*)(const char *))Tcl_GetAssocData(interp, "elle::libdmg::newDamageWrap", NULL);
    this->degrade = *(newDamageWrap(tag));
    return 1;
  }

  int update(void){
    BaseElem::update();
    /*
    BaseState pres = {
      .se = BaseElem::getResistingForce(),
      .vp = BaseElem::vppres,
      .kt = BaseElem::getTangentStiff(),
      .ke = BaseElem::getInitialStiff(),
    };
    */

    // Vector &se = BaseElem::getResistingForce();
    // wrapper->degrade(past, pres);
    return 1;
  }

private:
  degrade_f degrade = NULL;
}; // class DMG_DegradingBeamColumn


