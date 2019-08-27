/*
 
Copyright (C) 2007 Lucas K. Wagner

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
*/

#ifndef RELPSEUDOPOTENTIAL_H_INCLUDED
#define RELPSEUDOPOTENTIAL_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Wavefunction.h"
#include "Force_fitter.h"
#include "Basis_function.h"
#include "Pseudopotential.h"
class Program_options;
class System;
class Wavefunction;
class Sample_point;
class Wavefunction_data;
class Wavefunction_storage;


/*!

*/
class RelPseudopotential: public Pseudopotential
{
public:
  RelPseudopotential():maxaip(85)
  {deterministic=0;}
  

  void read(vector < vector <string> > & pseudotext,
            System * sys);
  /*!
    \brief
    Standard calculation of the pseudopotential.

    Uses random evaluation of the pseudopotential automatically.
   */

  void calcNonlocTmove(Wavefunction_data * wfdata, System *,
                       Sample_point * sample,
                       Wavefunction * wf,
                       Array1 <doublevar> & totalv,  //total p.e. from the psp
                       vector <Tmove> & tmoves  //variables for T-moves of Casula
                       );
  /*!
    \brief
    Provide your own random numbers for random evaluation of psp

    If all of them are set to one, everything is evaluated.  If it's set
    to zero, nothing is evaluated.
  */
  void calcNonlocWithTest(Wavefunction_data *, System *,Sample_point *, Wavefunction *,
                          const Array1 <doublevar> & accept_var,
                          Array1 <doublevar> & totalv);

  /*!
    \brief
    Uses file created by initializeStatic to calculated the nonlocal energy

    Uses a cutoff radius for evaluation, completely nonrandom
   */
  void calcNonlocWithFile(Wavefunction_data *, System *,Sample_point *, Wavefunction *,
                          Array1 <doublevar> &, Pseudo_buffer & input);


  
  void calcNonlocParmDeriv(Wavefunction_data * wfdata, System *,
                                            Sample_point * sample,
                                            Wavefunction * wf,
                                            const Array1 <doublevar> & accept_var,
                                            Array1 <doublevar> & totalv, Array1 <doublevar> & parm_deriv);
    
  /*!
    The worker function; the rest just provide simple defaults when functions don't need everything
   
    */
  void calcNonlocWithAllvariables(Wavefunction_data * wfdata, System *,
                                  Sample_point * sample,
                                  Wavefunction * wf,
                                  const Array1 <doublevar> & accept_var, //random variables for stochastic evaluation
                                  Array1 <doublevar> & totalv,  //total p.e. from the psp
                                  bool do_tmoves,vector <Tmove> & tmoves,  //variables for T-moves of Casula
                                  bool parm_derivatives, Array1 <doublevar> & parm_deriv //derivatives wrt wf parameters
                                  );

  ~RelPseudopotential();


private:
  int deterministic;

  void getRadial(int at, int spin, Sample_point * sample,
                                Array1 <doublevar> & r, Array1 <doublevar> & v_l);
  void getRadial(int at, int spin, Sample_point * sample,
                                  Array1 <doublevar> & r, 
                                  Array2 <doublevar> & v_l);
  const int maxaip; //!< Maximum number of atomic integration points
  int nelectrons;
  Array1 <int> aip;
  Array3 <doublevar> integralpt;
  Array3 <doublevar> integralpt_orig;
  Array2 <doublevar> integralweight;
  Array1 <doublevar> cutoff;
  Array1 <int> numL;
  vector <string> atomnames;
  Array1 <bool> addzeff; //!< whether or not to add Z_eff/r to the local function

  Storage_container wfStore;

  Array2 <Basis_function *> radial_basis;

};

#endif //PSEUDOPOTENTIAL_H_INCLUDED
//------------------------------------------------------------------------
