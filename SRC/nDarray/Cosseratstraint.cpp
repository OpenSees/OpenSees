//################################################################################
//# COPY-YES  (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           strain tensor with all necessery functions                #
//# CLASS:             Cosseratstraintensor                                     #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              July 25 '93                                               #
//# UPDATE HISTORY:    December 15 '93 replaced polinomial root solver for       #
//#                    principal straines with explicit formula                  #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//
#ifndef COSSERATSTRAINTENSOR_CC
#define COSSERATSTRAINTENSOR_CC

#include "Cosseratstraint.h"

//##############################################################################
Cosseratstraintensor::Cosseratstraintensor (int rank_of_tensor, double initval):
  tensor(rank_of_tensor, Cosserat_def_dim_2, initval)
    {   } // default constructor

//##############################################################################
Cosseratstraintensor::Cosseratstraintensor ( double *values ):
  tensor( 2, Cosserat_def_dim_2, values)
    {  }

//##############################################################################
Cosseratstraintensor::Cosseratstraintensor ( double initvalue ):
  tensor( 2, Cosserat_def_dim_2, initvalue) {  }

//##############################################################################
Cosseratstraintensor::Cosseratstraintensor( const Cosseratstraintensor & x ):
  tensor("NO")
    {
      x.pc_nDarray_rep->n++;  // tell the rval it has another reference
//      x.reference_count(+1);              // we're adding another reference.
      pc_nDarray_rep = x.pc_nDarray_rep;  // point to the new tensor_rep.
// add the indices
      indices1 = x.indices1;
      indices2 = x.indices2;
    }


//##############################################################################
Cosseratstraintensor::Cosseratstraintensor(const tensor & x):
  tensor( x ) {  } // copy-initializer

//##############################################################################
Cosseratstraintensor::Cosseratstraintensor(const nDarray & x):
  tensor( x ) {  }  // copy-initializer


//#//##############################################################################
//#Cosseratstraintensor::Cosseratstraintensor(Cosseratstraintensor & x)
//#  {
//#    x.reference_count(+1);              // we're adding another reference.
//#    pc_nDarray_rep = x.pc_nDarray_rep;  // point to the new tensor_rep.
//#// add the indices
//#    indices1 = x.indices1;
//#    indices2 = x.indices2;
//#  }


 // IT IS NOT INHERITED so must be defined in all derived classes
 // See ARM page 277.
 //##############################################################################
// Cosseratstraintensor::~Cosseratstraintensor()
// {
//   if (reference_count(-1) == 0)  // if reference count  goes to 0
//     {
// // DEallocate memory of the actual nDarray
// //    delete [pc_nDarray_rep->pc_nDarray_rep->total_numb] pc_nDarray_rep->pd_nDdata;
// // nema potrebe za brojem clanova koji se brisu## see ELLIS & STROUSTRUP $18.3
// //                                                and note on the p.65($5.3.4)
//     delete [] data();
//     delete [] dim();
//     delete pc_nDarray_rep;
//   }
// }
//

// IT IS NOT INHERITED so must be defined in all derived classes
// See ARM page 306.
//##############################################################################
Cosseratstraintensor Cosseratstraintensor::operator=( const Cosseratstraintensor & rval)
{
    rval.pc_nDarray_rep->n++;  // tell the rval it has another reference
//    rval.reference_count(+1);  // tell the rval it has another reference
//   /*  It is important to increment the reference_counter in the new
//       tensor before decrementing the reference_counter in the
//       old tensor_rep to ensure proper operation when assigning a
//       tensor_rep to itself ( after ARKoenig JOOP May/June '90 )  */

 // clean up current value;
//    if(--pc_nDarray_rep->n == 0)  // if nobody else is referencing us.
    if( reference_count(-1) == 0)  // if nobody else is referencing us.
      {
// DEallocate memory of the actual tensor
//      delete [pc_tensor_rep->pc_tensor_rep->total_numb] pc_tensor_rep->pd_nDdata;
// nema potrebe za brojem clanova koji se brisu## see ELLIS & STROUSTRUP $18.3
//                                                and note on the p.65($5.3.4)
//        delete  pc_nDarray_rep->pd_nDdata;
        delete [] data();
        delete [] dim();
// ovo ne smem da brisem jer nije dinamicki alocirano
//        delete pc_tensor_rep->indices;
        delete pc_nDarray_rep;
      }

 // connect to new value
    pc_nDarray_rep = rval.pc_nDarray_rep;  // point at the rval tensor_rep
// Temporary out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// null indices in the rval
//    rval.indices1 = NULL;
//    rval.indices2 = NULL;
//    rval.null_indices();
// Temporary out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    this->null_indices();
    return *this;
}




// IT IS NOT INHERITED so must be defined in all derived classes
// See ARM page 306.
//##############################################################################
Cosseratstraintensor Cosseratstraintensor::operator=( const tensor & rval)
{
    rval.pc_nDarray_rep->n++;  // tell the rval it has another reference
//    rval.reference_count(+1);  // tell the rval it has another reference
//   /*  It is important to increment the reference_counter in the new
//       tensor before decrementing the reference_counter in the
//       old tensor_rep to ensure proper operation when assigning a
//       tensor_rep to itself ( after ARKoenig JOOP May/June '90 )  */

 // clean up current value;
    if( reference_count(-1) == 0)  // if nobody else is referencing us.
      {
// DEallocate memory of the actual tensor
//      delete [pc_tensor_rep->pc_tensor_rep->total_numb] pc_tensor_rep->pd_nDdata;
// nema potrebe za brojem clanova koji se brisu## see ELLIS & STROUSTRUP $18.3
//                                                and note on the p.65($5.3.4)
//        delete  pc_nDarray_rep->pd_nDdata;
        delete [] data();
        delete [] dim();
// ovo ne smem da brisem jer nije dinamicki alocirano
//        delete pc_tensor_rep->indices;
        delete pc_nDarray_rep;
      }

 // connect to new value
    pc_nDarray_rep = rval.pc_nDarray_rep;  // point at the rval tensor_rep
// Temporary out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// null indices in the rval
//    rval.indices1 = NULL;
//    rval.indices2 = NULL;
//    rval.null_indices();
// Temporary out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    this->null_indices();
    return *this;
}

// IT IS NOT INHERITED so must be defined in all derived classes
// See ARM page 306.
//##############################################################################
Cosseratstraintensor Cosseratstraintensor::operator=( const nDarray & rval)
{
    rval.pc_nDarray_rep->n++;  // tell the rval it has another reference
//    rval.reference_count(+1);  // tell the rval it has another reference
//   /*  It is important to increment the reference_counter in the new
//       tensor before decrementing the reference_counter in the
//       old tensor_rep to ensure proper operation when assigning a
//       tensor_rep to itself ( after ARKoenig JOOP May/June '90 )  */

    if( reference_count(-1) == 0)  // if nobody else is referencing us.
      {
        delete [] data();
        delete [] dim();
        delete pc_nDarray_rep;
      }

 // connect to new value
    pc_nDarray_rep = rval.pc_nDarray_rep;  // point at the rval tensor_rep
    return *this;
}

//##############################################################################
// makes a complete new copy of Cosseratstraintensor!!
Cosseratstraintensor Cosseratstraintensor::deep_copy(void)
  {
    return Cosseratstraintensor(this->data()); // call constructor and return it !
  }
//..//##############################################################################
//..// returns a pointer to this for a deep copy
//..Cosseratstraintensor Cosseratstraintensor::p_deep_copy(void)
//..  {
//..    return &this->deep_copy(); // call constructor and return it !
//..  }

//ini  //##############################################################################
//ini  // use "from" and initialize already allocated strain tensor from "from" values
//ini  void Cosseratstraintensor::Initialize( const Cosseratstraintensor & from )
//ini    {
//ini  // copy onlu data because everything else is default
//ini      for ( int i=0 ; i<pc_nDarray_rep->total_numb ; i++ )
//ini        this->pc_nDarray_rep->pd_nDdata[i] = from.pc_nDarray_rep->pd_nDdata[i] ;
//ini    }

//___//##############################################################################
//___//##############################################################################
//___//##############################################################################
//___Cosseratstraintensor & Cosseratstraintensor::operator()(short ir, short is, short it,
//___                                        short tr, short ts, short tt  )
//___// another overloading of operator() . . .  // overloaded for THREE arguments
//___  {
//___    short where = ir - 1;
//___          where = where*ts + is - 1;
//___          where = where*tt + it - 1;
//___
//___//::printf(" w=%ld ",where);
//___    Cosseratstraintensor *p_value = this + where;
//___    return (*p_value);
//___  }


//##############################################################################
// invariants of the strain tensor              // Chen W.F. "plasticity for
double Cosseratstraintensor::Iinvariant1()const       // Structural Engineers"
  {
    return (cval(1,1)+cval(2,2)+cval(3,3));
  }

//##############################################################################
double Cosseratstraintensor::Iinvariant2() const
  {
    return (cval(2,2)*cval(3,3)-cval(3,2)*cval(2,3)+
            cval(1,1)*cval(3,3)-cval(3,1)*cval(1,3)+
            cval(1,1)*cval(2,2)-cval(2,1)*cval(1,2));
  }

//##############################################################################
double Cosseratstraintensor::Iinvariant3() const
  {

    double I3 = cval(1,1)*cval(2,2)*cval(3,3) +
                cval(1,2)*cval(2,3)*cval(3,1) +
                cval(1,3)*cval(2,1)*cval(3,2) -
                cval(1,3)*cval(2,2)*cval(3,1) -
                cval(1,2)*cval(2,1)*cval(3,3) -
                cval(1,1)*cval(2,3)*cval(3,2) ;

    return I3;

//    return ( this->determinant());
  }



//##############################################################################
// invariants of the deviatoric strain tensor
double Cosseratstraintensor::Jinvariant1() const
  {
    return (0.0);
  }

//##############################################################################
double Cosseratstraintensor::Jinvariant2() const
  {
    double EPS = sqrt(d_macheps());
    double temp1 = (Iinvariant1()*Iinvariant1()-3.0*Iinvariant2())/3.0;
    if ( temp1 < 0.0 || fabs(temp1) < EPS )
      {                    // this is because it might be close
        temp1 = 0.0;       // to zero ( -1e-19 ) numericaly
      }                    // but sqrt() does not accept it
    return ( temp1 );      // as (-) and theoreticaly J2d>0
  }

//##############################################################################
double Cosseratstraintensor::Jinvariant3() const
  {
    return ( (2.0*Iinvariant1()*Iinvariant1()*Iinvariant1()-
              9.0*Iinvariant1()*Iinvariant2() +
              27.0*Iinvariant3())/27.0 );
  }



//##############################################################################
double Cosseratstraintensor::equivalent( ) const	  //Zhaohui added 09-02-2000
{   
    // Evaluating e_eq = sqrt( 2.0 * epsilon_ij * epsilon_ij / 3.0)
    Cosseratstraintensor pstrain =  *this;
    tensor temp  = pstrain("ij") * pstrain("ij");
    double tempd = temp.trace();
    double e_eq  = pow( 2.0 * tempd / 3.0, 0.5 );
    //cout << "e_eq = " << e_eq << endlnn;
    return e_eq;

}


//##############################################################################
Cosseratstraintensor Cosseratstraintensor::principal() const
  {

    Cosseratstraintensor ret;

    double p_     = this->p_hydrostatic();
    double q_     = this->q_deviatoric();
    double theta_ = this->theta();
//old    while ( theta_ >= 2.0*PI )
//old      theta_ = theta_ - 2.0*PI; // if bigger than full cycle
//old    while ( theta_ >= 4.0*PI/3.0 )
//old      theta_ = theta_ - 4.0*PI/3.0; // if bigger than four thirds of half cycle
//old    while ( theta_ >= 2.0*PI/3.0 )
//old      theta_ = theta_ - 2.0*PI/3.0; // if bigger than two third of half cycle
//old    while ( theta_ >= PI/3.0 )
//old      theta_ = 2.0*PI/3.0 - theta_; // if bigger than one third of half cycle


//    ret.report("ret");
//    double sqrtJ2D = q/3.0;
//    double 2osqrt3 = 2.0/sqrt(3.0);
    double temp = (2.0*q_)/3.0;

    double ct  = cos( theta_  );
    double ctm = cos( theta_ - 2.0*PI/3.0 );
    double ctp = cos( theta_ + 2.0*PI/3.0 );

    ret.val(1,1) = p_ + temp*ct;
    ret.val(2,2) = p_ + temp*ctm;
    ret.val(3,3) = p_ + temp*ctp;

//    ret.report("ret");
    return ret;
//..    static complex ac[4];               // Chen W.F. "plasticity for
//..    static complex roots[4];            // Structural Engineers"
//..    int polish = 1 ;                    // page 53
//..    int m = 3;
//..
//..    ac[0] = complex( -(this->Iinvariant3()), 0.0 );
//..    ac[1] = complex(  (this->Iinvariant2()), 0.0);
//..    ac[2] = complex( -(this->Iinvariant1()), 0.0);
//..    ac[3] = complex(  1.0, 0.0);
//..
//..// what was obtained for coefficients
//..//DEBUGprint ::printf("ac[0].r = %lf  ac[0].i = %lf\n", real(ac[0]), imag(ac[0]));
//..//DEBUGprint ::printf("ac[1].r = %lf  ac[1].i = %lf\n", real(ac[1]), imag(ac[1]));
//..//DEBUGprint ::printf("ac[2].r = %lf  ac[2].i = %lf\n", real(ac[2]), imag(ac[2]));
//..//DEBUGprint ::printf("ac[3].r = %lf  ac[3].i = %lf\n", real(ac[3]), imag(ac[3]));
//..//DEBUGprint ::printf("m  = %d\n",m);
//..
//..    zroots( ac, m, roots, polish);
//..
//..//DEBUGprint ::printf("\nroots[1].r=%lf  roots[1].i=%lf\n",real(roots[1]),imag(roots[1]));
//..//DEBUGprint ::printf("roots[2].r=%lf  roots[2].i=%lf\n",real(roots[2]),imag(roots[2]));
//..//DEBUGprint ::printf("roots[3].r=%lf  roots[3].i=%lf\n",real(roots[3]),imag(roots[3]));
//..
//..
//..    Cosseratstraintensor principal(0.0);
//..
//..    principal.val(1,1) = real(roots[3]); // since they are sorted by
//..    principal.val(2,2) = real(roots[2]); // the zroot function in ascending
//..    principal.val(3,3) = real(roots[1]); // order . . .
//..                                         // sig1>sig2>sig3
//..    return principal;
//..
//..
  }

//##############################################################################
Cosseratstraintensor Cosseratstraintensor::deviator() const
  {
    tensor I2("I", 2, Cosserat_def_dim_2);
    Cosseratstraintensor st_vol = I2 * (trace()*(1./3.));
    Cosseratstraintensor st_dev = (*this) - st_vol;
    return st_dev;
  }



//##############################################################################
double Cosseratstraintensor::sigma_octahedral() const  // Chen W.F. "plasticity for
  {                                            // Structural Engineers"
    return ( this->Iinvariant1()/3.0 );        // page 59-60
  }

//##############################################################################
double Cosseratstraintensor::tau_octahedral() const    // Chen W.F. "plasticity for
  {                                             // Structural Engineers"
    return(sqrt(2.0/3.0*(this->Jinvariant2())));// page 59-60
  }



//##############################################################################
double Cosseratstraintensor::ksi()  const                     // Chen W.F. "plasticity for
  {                                            // Structural Engineers"
    return( (this->Iinvariant1())/sqrt(3.0) ); // page 66
  }


//##############################################################################
double Cosseratstraintensor::ro() const                        // Chen W.F. "plasticity for
  {                                             // Structural Engineers"
    double temp1 = this->Jinvariant2();         // page 68
    double EPS = pow(d_macheps(),(1./2.));
    if ( temp1 < 0.0 || fabs(temp1) < EPS )
      {
        temp1 = 0.0;
      }
    return( sqrt(2.0*(temp1)));
  }

//##############################################################################
double Cosseratstraintensor::p_hydrostatic() const         // Desai "Constitutive Laws
  {                                         // for Engineering Materials"
    // Joey modified to make it consistent with stress tensor
    return( - (this->Iinvariant1())*ONEOVERTHREE );  // page 283
    
    //return( (this->Iinvariant1()) );    // page 283
  }


//##############################################################################
double Cosseratstraintensor::q_deviatoric() const        // Desai "Constitutive Laws
  {                                       // for Engineering Materials"
//     double temp1 = this->Jinvariant2();   // page 283
//     return( sqrt(4.0/3.0*temp1) );

    double tempsqrt = 2./3.*(deviator()("ij")*deviator()("ij")).trace();
    double EPS = d_macheps();
// this is because it might be close
// to zero ( -1e-19 ) numericaly
// but sqrt() does not accept it
    if ( tempsqrt < 0.0 )
      {
::fprintf(stdout,"tempsqrt < 0.0 || fabs(tempsqrt) < EPS in ");
::fprintf(stdout," double Cosseratstraintensor::q_deviatoric() const\a\a\n");
::fprintf(stderr,"tempsqrt < 0.0 || fabs(tempsqrt) < EPS in ");
::fprintf(stderr," double Cosseratstraintensor::q_deviatoric() const \a\a\n");
        ::exit(1);
      }
    if ( fabs(tempsqrt) < EPS )
      {
        tempsqrt = 0.0;
      }
     double temp1 = sqrt(tempsqrt);
     return(temp1);

  }


//##############################################################################
double Cosseratstraintensor::theta()  const             // Chen W.F. "plasticity for
  {                                             // Structural Engineers"
                                                // page 70
//    double MULT = 1000000.0;
    double EPS = pow(d_macheps(),(1./2.));
    double temp1 = (3.0*sqrt(3.0)/2.0);
    double temp2 = (this->Jinvariant3());
    double temp3 = (this->Jinvariant2());
//....    double temp2 = (this->Jinvariant3())*MULT*MULT*MULT;
//....    double temp3 = (this->Jinvariant2())*MULT*MULT;
    if ( temp3 < 0.0 && fabs(temp3) < EPS )
      {
        temp3 = 0.0;
      }
    double temp4 = temp3 * temp3 * temp3;
    double temp5 = temp1 * temp2;
    double temp6 = sqrt(temp4);
    if ( (fabs(temp6)) <= fabs(temp5 * EPS) ) return ( 0.000001 );// slight perturbation because of 1/0 KERU06sep96
    double temp7 = temp5 / temp6;
    double tempabs1 = (fabs(temp7-1.0));
    double tempabs2 = (fabs(temp7+1.0));
    if ( tempabs1 < 0.001 ) return ( 0.000001 / 3.0 );// slight perturbation because of 1/0 KERU06sep96
//    if ( tempabs2 < 0.001 ) return ( PI / 3.0 );
    if ( tempabs2 < 0.001 ) return ( 3.14159/3.0 );// slight perturbation because of 1/0 KERU06sep96
    if ( temp7>1.0 || temp7<-1.0 )
      {
        ::fprintf(stderr,"\a\n something is wrong in Cosseratstraintensor::theta() (temp7>1.0||temp7<-1.0)\n");

        ::fprintf(stderr,"temp1 = %.20e\n", temp1);
        ::fprintf(stderr,"temp2 = %.20e\n", temp2);
        ::fprintf(stderr,"temp3 = %.20e\n", temp3);
        ::fprintf(stderr,"temp4 = %.20e\n", temp4);
        ::fprintf(stderr,"temp5 = %.20e\n", temp5);
        ::fprintf(stderr,"temp6 = %.20e\n", temp6);
        ::fprintf(stderr,"temp7 = %.20e\n", temp7);

        ::fprintf(stdout,"\a\n something is wrong in Cosseratstraintensor::theta() (temp7>1.0||temp7<-1.0)\n");
        this->print("s","Cosseratstraintensor s");
        ::fprintf(stdout,"temp1 = %.20e\n", temp1);
        ::fprintf(stdout,"temp2 = %.20e\n", temp2);
        ::fprintf(stdout,"temp3 = %.20e\n", temp3);
        ::fprintf(stdout,"temp4 = %.20e\n", temp4);
        ::fprintf(stdout,"temp5 = %.20e\n", temp5);
        ::fprintf(stdout,"temp6 = %.20e\n", temp6);
        ::fprintf(stdout,"temp7 = %.20e\n", temp7);

        exit (1);
      }
    double temp8 = acos(temp7);
    double temp9 = temp8 / 3.0;

    return ( temp9 );
  }

//##############################################################################
double Cosseratstraintensor::thetaPI() const
  {
    double thetaPI = theta() / PI;
    return thetaPI;
  }



//##############################################################################
Cosseratstraintensor Cosseratstraintensor::pqtheta2strain( double p, double q, double theta)
  {
    Cosseratstraintensor ret;

//    double sqrtJ2D = q/3.0;
//    double 2osqrt3 = 2.0/sqrt(3.0);
    double temp = (2.0*q)/(3.0);

    while ( theta >= 2.0*PI )
      theta = theta - 2.0*PI; // if bigger than full cycle
    while ( theta >= 4.0*PI/3.0 )
      theta = theta - 4.0*PI/3.0; // if bigger than four thirds of half cycle
    while ( theta >= 2.0*PI/3.0 )
      theta = theta - 2.0*PI/3.0; // if bigger than two third of half cycle
    while ( theta >= PI/3.0 )
      theta = 2.0*PI/3.0 - theta; // if bigger than one third of half cycle

    double ct  = cos( theta  );
    double ctm = cos( theta - 2.0*PI/3.0 );
    double ctp = cos( theta + 2.0*PI/3.0 );

    ret.val(1,1) = p + temp*ct;
    ret.val(2,2) = p + temp*ctm;
    ret.val(3,3) = p + temp*ctp;

    return ret;

  }


//##############################################################################
Cosseratstraintensor Cosseratstraintensor::evoleq2strain( double evol, double eq )
  {
    Cosseratstraintensor ret;

    ret.val(1,1) = 1./3.*evol + 1./2.*eq;
    ret.val(2,2) = 1./3.*evol + 1./2.*eq;
    ret.val(3,3) = 1./3.*evol - eq;

    return ret;

  }



//##############################################################################
void Cosseratstraintensor::report(char * msg) const
  {
    ::printf("\n****************  strain tensor report ****************\n");
    if ( msg ) ::printf("%s",msg);

    this->print("st","Cosseratstraintensor st");

    ::printf("I1 = %.8e ; I2 = %.8e ; I3 = %.8e \n",
              Iinvariant1(),Iinvariant2(),Iinvariant3());

    printf("st_trace = %.8e,  mean pressure p = %.8e\n",
             trace(),  trace()/3.0);

    tensor I2("I", 2, Cosserat_def_dim_2);

    Cosseratstraintensor st_vol = I2 * trace() * (1./3.);
    st_vol.print("st_v","tensor st_vol (volumetric part of the st tensor)");

    Cosseratstraintensor st_dev = this->deviator();  // - st_vol;
    st_dev.print("st_d","tensor st_dev (strain deviator)");

    ::printf("J1 = %.8e ; J2 = %.8e ; J3 = %.8e ;\n",
              Jinvariant1(),Jinvariant2(),Jinvariant3());
//...  Cosseratstraintensor Jinv2 = st_dev("ij")*st_dev("ji")*0.5;

    Cosseratstraintensor st_principal = principal();
    st_principal.print("st_p","principal strain tensor");


    ::printf("sig_oct = %.8e  , tau_oct = %.8e\n",
              sigma_octahedral(), tau_octahedral());


    ::printf("ksi=%.6e, ro=%.6e, theta=%.6e=%.6e*PI \n",
              ksi(),    ro(),    theta(),  thetaPI());

    if ( msg ) ::printf("%s",msg);
    ::printf("\n############  end of strain tensor report ############\n");
  }


//##############################################################################
void Cosseratstraintensor::reportshort(char * msg) const
  {
//    ::printf("\n         ****************** short strain tensor report ***\n");
//    if ( msg ) ::printf("         %s",msg);

    this->print("st"," ");

//    ::printf("ksi = %.8e ,ro = %.8e ,theta = %.8e\n",
//              ksi(),       ro(),      theta());

//    ::printf("p=%.12e , q=%.12e , theta=%.12e*PI\n",
//              p_hydrostatic(), q_deviatoric(), thetaPI());

  }

#endif
