/*
 )
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


#include "RelPseudopotential.h"
#include "qmc_io.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include <iomanip>
#include "ulec.h"
#include "Wavefunction_data.h"
#include "Basis_function.h"
#include "System.h"
#include "Array45.h"
using namespace std;
/*!
Some legendre polynomials..
*/


dcomplex SphericalHarmonic(int l,int m, const Array1<doublevar> & r) {
  assert(r.GetDim(0) == 5);
  switch (l) {
  case 0:
    return 1;
    break;

  case 1:
    switch(m) {
    case -1:
      return sqrt(3./2.)*dcomplex(r(2)/r(0),-r(3)/r(0));
      break;
    case 0:
      return sqrt(3.)*dcomplex(r(4)/r(0),0.0);
      break;
    case 1:
      return -sqrt(3./2.)*dcomplex(r(2)/r(0),r(3)/r(0));
      break;
    default:
      return 0.0;
      break;
    }
    break;

  case 2:
    switch(m) {
    case -2:
      return sqrt(15./8.)*(dcomplex(r(2),-r(3))*dcomplex(r(2),-r(3)))/r(1);
      break;
    case -1:
      return sqrt(15./2.)*dcomplex(r(2)*r(4),-r(3)*r(4))/r(1);
      break;
    case 0:
      return sqrt(5./4.)*dcomplex(2*r(4)*r(4)-r(2)*r(2)-r(3)*r(3),0.0)/r(1);
      break;
    case 1:
      return -sqrt(15./2.)*dcomplex(r(2)*r(4),r(3)*r(4))/r(1);
      break;
    case 2:
      return sqrt(15./8.)*(dcomplex(r(2),r(3))*dcomplex(r(2),r(3)))/r(1);
      break;
    default:
      return 0.0;
      break;
    }
    break;

  case 3:
    switch(m) {
    case -3:
      return sqrt(35./16.)*pow(dcomplex(r(2),-r(3)),3)/(r(0)*r(1));
      break;
    case -2:
      return sqrt(105./8.)*(dcomplex(r(2),-r(3))*dcomplex(r(2),-r(3)))*r(4)/(r(0)*r(1));
      break;
    case -1:
      return sqrt(21./16.)*dcomplex(r(2),-r(3))*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1));
      break;
    case 0:
      return sqrt(7./4.)*dcomplex(r(4)*(2*r(4)*r(4)-3*r(2)*r(2)-3*r(3)*r(3)),0.0)/(r(0)*r(1));
      break;
    case 1:
      return -sqrt(21./16.)*dcomplex(r(2),r(3))*(4*r(4)*r(4)-r(2)*r(2)-r(3)*r(3))/(r(0)*r(1));
      break;
    case 2:
      return sqrt(105./8.)*(dcomplex(r(2),r(3))*dcomplex(r(2),r(3)))*r(4)/(r(0)*r(1));
      break;
    case 3:
      return -sqrt(35./16.)*pow(dcomplex(r(2),r(3)),3)/(r(0)*r(1));
      break;
    default:
      return 0.0;
    }
    break;

  case 4:
    switch(m) {
    case -4:
      return 3*sqrt(35./128.)*pow(dcomplex(r(2),-r(3)),4)/(r(1)*r(1));
      break;
    case -3:
      return 3*sqrt(35./16.)*pow(dcomplex(r(2),-r(3)),3)*r(4)/(r(1)*r(1));
      break;
    case -2:
      return 3*sqrt(5./32.)*(dcomplex(r(2),-r(3))*dcomplex(r(2),-r(3)))*(7*r(4)*r(4)-r(1))/(r(1)*r(1));
      break;
    case -1:
      return 3*sqrt(5./16.)*(dcomplex(r(2),-r(3))*r(4)*(7*r(4)*r(4)-3*r(1)))/(r(1)*r(1));
      break;
    case 0:
      return 3*sqrt(1./64.)*dcomplex(35*pow(r(4),4)-30*r(4)*r(4)*r(1) + 3*r(1)*r(1),0)/(r(1)*r(1));
      break;
    case 1:
      return -3*sqrt(5./16.)*(dcomplex(r(2),r(3))*r(4)*(7*r(4)*r(4)-3*r(1)))/(r(1)*r(1));
      break;
    case 2:
      return 3*sqrt(5./32.)*(dcomplex(r(2),r(3))*dcomplex(r(2),r(3)))*(7*r(4)*r(4)-r(1))/(r(1)*r(1));
      break;
    case 3:
      return -3*sqrt(35./16.)*pow(dcomplex(r(2),r(3)),3)*r(4)/(r(1)*r(1));
      break;
    case 4:
      return 3*sqrt(35./128.)*pow(dcomplex(r(2),r(3)),4)/(r(1)*r(1));
      break;
    default:
      return 0.0;
      break;
    }
    break;

  default:
    error("Spherical Harmonics not implemented");
    break;
  }

}

dcomplex cSphericalHarmonic(int l, int m, const Array1<doublevar> & r) {
  return conj(SphericalHarmonic(l,m,r));
}


dcomplex Yjml(doublevar j, doublevar mj, int l, const Array1<doublevar> & r, doublevar s) {
  dcomplex up,dn;
  if (j==l+0.5) {
    up = sqrt((l+mj+0.5)/(2.*l+1.))*SphericalHarmonic(l,int(mj-0.5),r)*exp(I*s);
    dn = sqrt((l-mj+0.5)/(2.*l+1.))*SphericalHarmonic(l,int(mj+0.5),r)*exp(-I*s);
  }
  else if (j==l-0.5) {
    up = -sqrt((l-mj+0.5)/(2.*l+1.))*SphericalHarmonic(l,int(mj-0.5),r)*exp(I*s);
    dn = sqrt((l+mj+0.5)/(2.*l+1.))*SphericalHarmonic(l,int(mj+0.5),r)*exp(-I*s);
  }
  else {
    error("Invalid j for spin spherical harmonic");
  }
  return up + dn;
}

dcomplex cYjmlSeparated(doublevar j, doublevar mj, int l, const Array1<doublevar> & r, int s) {
  dcomplex val;
  if (j == l+0.5) {
    if (s==0) val = sqrt((l+mj+0.5)/(2.*l+1.))*cSphericalHarmonic(l,int(mj-0.5),r);
    else val = sqrt((l-mj+0.5)/(2.*l+1.))*cSphericalHarmonic(l,int(mj+0.5),r);
  }
  else if (j==l-0.5) {
    if (s==0) val = -sqrt((l-mj+0.5)/(2.*l+1.))*cSphericalHarmonic(l,int(mj-0.5),r);
    else val = sqrt((l+mj+0.5)/(2.*l+1.))*cSphericalHarmonic(l,int(mj+0.5),r);
  }
  else {
    error("Invalid j for spin spherical harmonic conjugate");
  }
  return val;
}



RelPseudopotential::~RelPseudopotential()  {
  for(int i=0; i< radial_basis.GetDim(0); i++)
    for(int j=0; j < radial_basis.GetDim(1); j++)  
      if(radial_basis(i,j)) delete radial_basis(i,j);
}
//----------------------------------------------------------------------



void RelPseudopotential::calcNonlocTmove(Wavefunction_data * wfdata, System * sys,
                     Sample_point * sample,
                     Wavefunction * wf,
                     Array1 <doublevar> & totalv,  //total p.e. from the psp
                     vector <Tmove> & tmoves  //variables for T-moves of Casula
                     ) { 
  int tot=nTest();
  Array1 <doublevar> test(tot);
  for(int i=0; i< tot; i++) {
    test(i)=rng.ulec();
  }
  Array1 <doublevar> parm_deriv;
  calcNonlocWithAllvariables(wfdata,sys,sample, wf, test,totalv, true, tmoves,false, parm_deriv);
}
//----------------------------------------------------------------------

void RelPseudopotential::calcNonlocWithTest(Wavefunction_data *wfdata , System * sys, 
                                         Sample_point * sample, Wavefunction *wf ,
                                         const Array1 <doublevar> & accept_var,
                                         Array1 <doublevar> & totalv) { 
  vector<Tmove>  tmoves;
  Array1 <doublevar> parm_deriv;
  calcNonlocWithAllvariables(wfdata,sys, sample, wf, accept_var,totalv, false, tmoves,false, parm_deriv);
}

void RelPseudopotential::calcNonlocParmDeriv(Wavefunction_data * wfdata, System * sys,
                                          Sample_point * sample,
                                          Wavefunction * wf,
                                          const Array1 <doublevar> & accept_var,
                                          Array1 <doublevar> & totalv, Array1 <doublevar> & parm_deriv) { 
  vector<Tmove>  tmoves;
  calcNonlocWithAllvariables(wfdata,sys,sample, wf, accept_var,totalv, false, tmoves,true, parm_deriv);
  
}

void RelPseudopotential::calcNonlocWithAllvariables(Wavefunction_data * wfdata,
                                                 System * sys,
                                                 Sample_point * sample,
                                                 Wavefunction * wf,
                                                 const Array1 <doublevar> & accept_var,
                                                 Array1 <doublevar> & totalv,
                                                 bool do_tmoves,vector <Tmove> & tmoves,
                                                 bool parm_derivatives, Array1 <doublevar> & parm_deriv
                                          )
{
  int natoms=sample->ionSize();
  int nwf=wf->nfunc();
  assert(accept_var.GetDim(0) >= nTest());
  assert(totalv.GetDim(0) >= nwf);
  assert(nelectrons == sample->electronSize());
  Array1 <doublevar> ionpos(3), oldpos(3), newpos(3);
  Array1 <doublevar> newdist(5), olddist(5);

  totalv=0;
  doublevar accum_local=0;
  doublevar accum_nonlocal=0;
  wf->updateVal(wfdata, sample);
  wfStore.initialize(sample, wf);
  Parm_deriv_return base_deriv;
  if(parm_derivatives) {
    error("parm_derivatives not implemented in RelPseudopotential");
    /*
    parm_deriv.Resize(wfdata->nparms());
    parm_deriv=0;
    base_deriv.nparms_start=0;
    base_deriv.nparms_end=wfdata->nparms();
    base_deriv.need_hessian=0;    
    wf->getParmDeriv(wfdata, sample, base_deriv);
    */
  }
  int accept_counter=0;
  Array1 <doublevar>  nonlocal(nwf);
  Array2 <dcomplex> integralpts(nwf, maxaip);
//minyi: save old inverse transpose matrix
  Array2 <Array2 <dcomplex> > inverse_old;
  wf->getInverseTranspose(wfdata, inverse_old); //get transpose inverse matrix
  Array1 <doublevar> ori_s;
  for(int at=0; at< natoms; at++){
    if(numL(at) != 0) {
      Array1 <doublevar> v_l(numL(at));
      
      sample->getIonPos(at, ionpos);
      ori_s.Resize(sample->electronSize());
      for(int e=0; e < sample->electronSize(); e++)  {
        sample->getElectronSpin(e,ori_s(e));
        sample->getElectronPos(e, oldpos);
        sample->updateEIDist();
        sample->getEIDist(e,at, olddist);
        nonlocal=0;
        int spin=0;
        getRadial(at,spin, sample, olddist, v_l);
        //Start integral
        int accept;
        if(deterministic) {
          accept= olddist(0) < cutoff(at);
        }
        else {
          doublevar strength=0;
          const doublevar calculate_threshold=10;
          for(int l=0; l<numL(at)-1; l++)
            strength+=calculate_threshold*(2*l+1)*fabs(v_l(l));
   
          strength=min((doublevar) 1.0, strength);
          doublevar rand=accept_var(accept_counter++);
          for(int l=0; l<numL(at)-1; l++)
            v_l(l)/=strength;
          accept=strength>rand;
        }
//-------------------------------------------------------
//read Molecorb
        if(accept)  {
          wfStore.saveUpdate(sample, wf, e);
          Array4 <dcomplex> Val;
	        Array4 <dcomplex> oldWfVal;
          wf->getSpinorComponents(wfdata,sample,oldWfVal,e);
          int nmo=oldWfVal.GetDim(2);
          int ndet=oldWfVal.GetDim(1); //trick to get nmo and ndet
          Wf_return WfVal_old(nwf,2);
          wf->getVal(wfdata, e,WfVal_old);
//cout<<exp(WfVal_old.amp(0,0))*exp(I*WfVal_old.phase(0,0))*inverse_old(0,0,0)(0,0)<<endl;
//nlj is the number of the projector,nmo #of mo,=n(0)+n(1)=n(0)
          for(int i=0; i< aip(at); i++) {
            sample->setElectronPos(e, oldpos);
            doublevar base_sign=sample->overallSign();
            doublevar base_phase=sample->overallPhase();
            for(int d=0; d < 3; d++) 
              newpos(d)=integralpt(at,i,d)*olddist(0)-olddist(d+2);
//  cout<<olddist(2)<<" "<<olddist(3)<<" "<<olddist(4)<<" "<<oldpos(0)<<" "<<oldpos(1)<<" "<<oldpos(2)<<endl;
//  cout<<integralpt(at,i,0)<<" "<<integralpt(at,i,1)<<" "<<integralpt(at,i,2)<<endl;
            sample->translateElectron(e, newpos);
            sample->setElectronSpin(e,ori_s(e));
            sample->updateEIDist();
            sample->getEIDist(e,at,newdist);
            doublevar new_sign=sample->overallSign();
            doublevar new_phase=sample->overallPhase();
            wf->updateVal(wfdata, sample);
            Wf_return WfVal_new(nwf,2);
            wf->getVal(wfdata, e,WfVal_new);
            wf->getSpinorComponents(wfdata, sample, Val, e);
//------------------------------------------------------- 
            dcomplex A_up;dcomplex B_do;

            assert(numL(at)%2==0);
            for (int n=0; n < numL(at)-1; n++) {
              int l = (n%2+n)/2;
              doublevar j;
              if (n%2==1) j = l+0.5;
              else if (n%2==0 && n>0) j = l-0.5;
              else if (n==0) j = l+0.5;
              else error("Invalid value for l,j");
              dcomplex pup,pdn;
              for (doublevar mj = -j; mj <= j; mj += 1) {
                //The cYjmlSeparated has already carried out spin integral, so only need the
                //term that doesn't vanish
                pup += Yjml(j,mj,l,olddist,ori_s(e))*cYjmlSeparated(j,mj,l,newdist,0);
                pdn += Yjml(j,mj,l,olddist,ori_s(e))*cYjmlSeparated(j,mj,l,newdist,1);
              }
              A_up += pup*v_l(n);
              B_do += pdn*v_l(n);
            }

           for(int w=0; w< nwf; w++) { 
             dcomplex ratio1=0;
             dcomplex ratio2=0;//for spin up and down chanel
             for(int det=0;det<ndet;det++){
              for(int n=0;n<nmo;++n) {
                ratio1=ratio1+Val(w,det,n,0)*(inverse_old(w,det)(e,n));
                ratio2=ratio2+Val(w,det,n,1)*(inverse_old(w,det)(e,n)); 
               }
              }
              integralpts(w,i)=0.0;
              integralpts(w,i)=(A_up*ratio1+B_do*ratio2)*integralweight(at, i)*Val(w,0,0,2)
                               *exp(-I*WfVal_old.phase(w,0))
                               /exp(WfVal_old.amp(w,0));
             }
//------------------------------------------------------- 
            
            for(int w=0; w< nwf; w++)  {
             if(!do_tmoves) { 
                nonlocal(w)+=integralpts(w,i).real();
              }
              else { 
                error("Tmove not supported");
              }
              //-----------parameter derivatives
              if(parm_derivatives) { 
                error("parm_derivatives not supported in this ecp evaluation");
                /*
                Parm_deriv_return deriv;
                deriv.nparms_start=0;
                deriv.nparms_end=wfdata->nparms();
                deriv.need_hessian=0;
                wf->getParmDeriv(wfdata, sample, deriv);
                int np=wfdata->nparms();
                for(int p=0; p < np; p++) { 
                  parm_deriv(p)+=(deriv.gradient(p)-base_deriv.gradient(p))*integralpts(w,i).real();
                }
                */
             }
              //------
            }
            sample->setElectronPos(e, oldpos);
          } 
          //--------------------
          wfStore.restoreUpdate(sample, wf, e);
          sample->setElectronSpin(e,ori_s(e));
        }
        //----------------------------------------------
        //now do the local part
        doublevar vLocal=0;
        int localL=numL(at)-1; //The l-value of the local part is
                               //the last part.
        vLocal=v_l(localL);
        accum_local+=vLocal;
        accum_nonlocal+=nonlocal(0);
  // cout << "atom " << at << " r " << olddist(0) <<   " vLocal  " << vLocal
  // << "    nonlocal   " << nonlocal(0) << endl<<endl;
        for(int w=0; w< nwf; w++) {
          totalv(w)+=vLocal+nonlocal(w);
          //totalv(w)+=nonlocal(w); 
          //totalv(w)+=vLocal;
//totalv(w)=0;
        }
      }  //electron loop
    }  //if atom has any psp's
  }  //atom loop

  //cout << "psp: local part " << accum_local
  //<< "  nonlocal part " << accum_nonlocal << endl;
}


//------------------------------------------------------------------------

#include "System.h"

/*!
*/
void RelPseudopotential::read(vector <vector <string> > & pseudotext,
                           System * sys)
{

  sys->getAtomicLabels(atomnames);
  int natoms=atomnames.size();
  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  //Get the atomic integration points
  aip.Resize(natoms);
  aip=6; //default value for atomic integration points
  integralpt.Resize(natoms, maxaip, 3);
  integralpt_orig.Resize(natoms, maxaip, 3);
  integralweight.Resize(natoms, maxaip);
  addzeff.Resize(natoms);
  addzeff=false;
  Array1 <int> atom_has_psp(natoms);
  atom_has_psp=0;
  int maxL=0;
  
  //Assuming we have spin 1/2 particles here, as often assumed..
  radial_basis.Resize(natoms,2);
  radial_basis=NULL;
  /*
  vector < vector < string > > basistxt;
  for(unsigned int i=0; i< pseudotext.size(); i++) {
    unsigned int pos=0;
    basistxt.push_back(basistmp);
  }
  */

  for(unsigned int i=0; i< pseudotext.size(); i++ ) {
    for(int at=0; at<natoms; at++) {
      unsigned int pos=0;
      if( pseudotext[i][0] == atomnames[at] ) {
        if(atom_has_psp(at))
          error("There are two PSEUDO sections for ", atomnames[at]);
        atom_has_psp(at)=1;
        vector <string> basistmp;
        if(haskeyword(pseudotext[i], pos=0, "SPIN_DEP")) { 
          if(!readsection(pseudotext[i], pos=0, basistmp, "BASIS_UP"))
            error("Need BASIS_DN section in pseudopotential for ", pseudotext[i][0]);
          allocate(basistmp, radial_basis(at,0));
          basistmp.clear();
          if(!readsection(pseudotext[i], pos=0, basistmp, "BASIS_DN"))
            error("Need BASIS_UP section in pseudopotential for ", pseudotext[i][0]);
          allocate(basistmp, radial_basis(at,1));
          
        }
        else { 
          if(!readsection(pseudotext[i], pos=0, basistmp, "BASIS"))
            error("Need Basis section in pseudopotential for ", pseudotext[i][0]);
          allocate(basistmp, radial_basis(at,0));
          allocate(basistmp, radial_basis(at,1));

        }
        
//        string type;
        //vector <string> tmpbasis(basistxt[i]);
        assert(radial_basis(at,0)->nfunc()==radial_basis(at,1)->nfunc());
        int nlval=radial_basis(at,0)->nfunc();
        if(maxL < nlval) maxL=nlval;
        readvalue(pseudotext[i], pos=0, aip(at), "AIP");
        if(haskeyword(pseudotext[i], pos=0, "ADD_ZEFF")) {
          addzeff(at)=true;
        }

      }
    }
  }



  for(int at=0; at<natoms; at++) {
    int aiptemp=aip(at);
    Array1 <doublevar> xpt(maxaip);
    Array1 <doublevar> ypt(maxaip);
    Array1 <doublevar> zpt(maxaip);
    Array1 <doublevar> weight(maxaip);


    gesqua(aiptemp, xpt, ypt, zpt, weight);

    aip(at)=aiptemp;

    for(int i=0; i< aip(at); i++)
    {
      integralpt(at,i,0)=integralpt_orig(at, i, 0)=xpt(i);
      integralpt(at,i,1)=integralpt_orig(at, i, 1)=ypt(i);
      integralpt(at,i,2)=integralpt_orig(at, i, 2)=zpt(i);
      integralweight(at, i)=weight(i);
    }

  }


  //allocate the storage

  numL.Resize(natoms);

  numL=0;
  for(int at=0; at < natoms; at++ )
  {
    if(radial_basis(at,0) != NULL)
      numL(at)=radial_basis(at,0)->nfunc();
  }


  //find the cutoff radius for a static calculation

  Sample_point * tempsample=NULL;
  sys->generateSample(tempsample);
  cutoff.Resize(natoms);
  const doublevar cutoff_threshold=1e-5;
  const doublevar cutoff_max=20.0;
  cutoff=0.0;
  const doublevar cutoff_interval=.05;
  for(int at=0; at< natoms; at++)
  {
    if(numL(at) > 0) {
    Array1 <doublevar> cutoffL(numL(at));
    Array1 <doublevar> v_l(numL(at));
    Array1 <doublevar> v_l2(numL(at));
    Array1 <doublevar> tempr(5);
    Array1 <int> foundcutoff(numL(at));
    foundcutoff=0;
    tempr=0;
    cutoffL=cutoff_max;
    for(doublevar r=cutoff_max; r> 0; r-=cutoff_interval)
    {
      //going along the z axis; watch this if we ever try non-spherical
      //potentials

      tempr(0)=r; tempr(1)=r*r; tempr(4)=r;
      getRadial(at,0, tempsample, tempr, v_l);
      getRadial(at,1, tempsample, tempr, v_l2);

      for(int l=0; l< numL(at)-1; l++)
      {
        if( (fabs(v_l(l)) > cutoff_threshold 
             || fabs(v_l2(l)) > cutoff_threshold) && !foundcutoff(l))
        {
          //It seems to me that it should be r+cutoff_interval,
          //but for compatability with the f90 code, we'll put it
          //at r-cutoff_interval.  Shouldn't be too critical, anyway.
          cutoffL(l)=r-cutoff_interval;
          foundcutoff(l)=1;
        }
      }
    }

    for(int l=0; l< numL(at)-1; l++)
    {
      cutoff(at)=max(cutoff(at), cutoffL(l));
    }
    }
  }
  delete tempsample;

}

void RelPseudopotential::getRadial(int at, int spin,Sample_point * sample,
                                Array1 <doublevar> & r, Array1 <doublevar> & v_l) {
  assert(radial_basis(at,spin) != NULL);
  radial_basis(at,spin)->calcVal(r, v_l);
  if(addzeff(at)) {
    //for(int l=0; l < numL(at); l++) {
    int l=numL(at)-1;
    doublevar cutoff_rad=radial_basis(at,spin)->cutoff(l);
    if(r(0) < cutoff_rad) {
      v_l(l) += sample->getIonCharge(at)/r(0);
    }
    //}
  }
}


void RelPseudopotential::getRadial(int at, int spin, Sample_point * sample,
                                Array1 <doublevar> & r, 
                                Array2 <doublevar> & v_l){
  assert(radial_basis(at,spin) != NULL);
  radial_basis(at,spin)->calcLap(r, v_l);
  if(addzeff(at)) {
    int l=numL(at)-1;
    doublevar cutoff_rad=radial_basis(at,spin)->cutoff(l);
    if(r(0) < cutoff_rad) {
      v_l(l,0) += sample->getIonCharge(at)/r(0);
      for(int d=0; d< 3; d++) {
        v_l(l,d+1) -= sample->getIonCharge(at)*r(d+2)/(r(0)*r(1));
      }
    }
  }
}

//----------------------------------------------------------------------


//------------------------------------------------------------------------
