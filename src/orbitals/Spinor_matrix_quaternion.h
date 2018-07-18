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

#ifndef SPINOR_MATRIX_QUATERNION_H_INCLUDED
#define SPINOR_MATRIX_QUATERNION_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"

class System;
class Sample_point;

struct quaternion {
    pair<dcomplex,dcomplex> val;
    quaternion() { this->val.first = dcomplex(0.0); this->val.second = dcomplex(0.0);}
    quaternion(dcomplex x, dcomplex y) {
	this->val.first = x;
	this->val.second = y;
    }
    quaternion(doublevar a, doublevar b, doublevar c, doublevar d) {
	dcomplex x,y;
	x = dcomplex(a,b);
	y = dcomplex(c,d);
	this->val.first = x;
	this->val.second = y;
    }
    doublevar norm() {
	doublevar tmp = 0.0;
	tmp += abs(this->val.first)*abs(this->val.first);
	tmp += abs(this->val.second)*abs(this->val.second);
	tmp = sqrt(tmp);
	return tmp;
    }
};

doublevar abs(quaternion x) {
    return x.norm();
}

//quaternion times int,real,or complex
template<class T>
quaternion operator*(T x, quaternion q) {
    quaternion tmp;
    tmp.val.first = x*q.val.first;
    tmp.val.second = x*q.val.second;
    return tmp;
}

template<class T>
quaternion operator*(quaternion q, T x) { return x*q; };

//Quaternion times quaternion
quaternion operator*(quaternion q, quaternion r) {
    quaternion tmp;
    tmp.val.first  = q.val.first*r.val.first - r.val.second*q.val.second;
    tmp.val.second = q.val.first*r.val.second + r.val.first*q.val.second;
    return tmp;
}

quaternion operator+(quaternion q, quaternion r) { q.val.first += r.val.first; q.val.second += r.val.second; return q; }
quaternion operator+=(quaternion q, quaternion r) { q = q + r; return q; }

template <class CharT,class Traits>
basic_ostream<CharT,Traits>& 
operator<<(basic_ostream<CharT,Traits>& os, const quaternion & x) {
    basic_ostringstream<CharT, Traits> s;
    s.flags(os.flags());
    s.imbue(os.getloc());
    s.precision(os.precision());
    s << "( " << x.val.first.real() << ", " << x.val.first.imag() << ", " << x.val.second.real() << ", " << x.val.second.imag() << ")";
    return os << s.str();
}

template<class CharT, class Traits>
basic_istream<CharT, Traits>&
operator>>(basic_istream<CharT, Traits>& is, quaternion & x) {
    doublevar a,b,c,d;
    CharT ch;
    is >> ch;
    if (ch == '(') {
	is >> a >> ch;
	if (ch == ',') {
	    is >> b >> ch;
	    if (ch == ',') {
		is >> c >> ch;
		if (ch == ',') {
		    is >> d >> ch;
		    if (ch == ')') {
			x = quaternion(a,b,c,d);
		    }
		    else
			is.setstate(ios_base::failbit);
		}
		else
		    is.setstate(ios_base::failbit);
	    }
	    else
		is.setstate(ios_base::failbit);
	}
	else
	    is.setstate(ios_base::failbit);
    }
    else
        is.setstate(ios_base::failbit);
    return is;
}

#ifdef USE_MPI
inline void overloaded_broadcast(Array1 <quaternion> & v) {
    //CM
    //quaternion is basically an ordered pair of complex numbers.
    //Trick, created Array1 <dcomplex> of double size and broadcast that way.
    //then simply reassemble back into quaternion array
    Array1 <dcomplex> w;
    w.Resize(2*v.GetDim(0));
    for (int i = 0; i < v.GetDim(0); i++) {
	w(2*i)   = v(i).val.first;
	w(2*i+1) = v(i).val.second;
    }
    overloaded_broadcast(w);
    for (int i = 0; i < v.GetDim(0); i++) {
	v(i).val.first  = w(2*i);
	v(i).val.second = w(2*i+1);
    }
}
#endif


//----------------------------------------------------------------------------

class Spinor_matrix_quaternion: public Complex_MO_matrix {
protected:
  void init();

private:
  //Center_set centers;
  //Array1 <Basis_function *> basis;
  //int nmo;
  //int totbasis;
  //int maxbasis;
 // doublevar magnification_factor;
  //string orbfile;
  Array2 <int> mofill;
  Array2 <quaternion> moCoeff2;
  Array1 <int> nbasis;

  Array1 <doublevar> obj_cutoff; //!< cutoff for each basis object
  Array1 <doublevar> cutoff;  //!< Cutoff for individual basis functions
  Array1 <int> nfunctions; //!< number of functions in each basis
  //Array1 <int> basismo;
  //Array2 <doublevar> moCoeff;
  //Array2 <int> basisfill;

  Array2 <int>  basisfill_list;
  Array2 <quaternion>  moCoeff_list;
  Array1 <int>  basismo_list;

 Array1 <doublevar> symmvals_temp1d;
 Array2 <doublevar> symmvals_temp2d;



public:

  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1< Array1 <int> > & occupations) {error("Spinors do not need a spin index"); }
  virtual void buildLists(Array1 <int> & occupations);

  /*!
    get the number of molecular orbitals
   */
  //virtual int getNmo()
  //{
  //  return nmo;
  //}

  virtual int showinfo(ostream & os);

  virtual int writeinput(string &, ostream &);


  //virtual void read(vector <string> & words, unsigned int & startpos, System * sys);

  //! Takes an ORB file and inserts all the coefficients.
  //virtual int readorb(istream &);


  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, Array1 <int> &);

  virtual void getMoCoeff(Array2 <dcomplex> & coeff) {
    error("Cutoff_Mo doesn't support optimization yet");
  }
  
  /*!
    
   */
  virtual void setMoCoeff(Array2 <dcomplex> & coeff) {
    error("Cutoff MO doesn't support optimization yet");
  }

  virtual int nMoCoeff() {
    error("Need to implement MO_matrix_quaternion::nMoCoeff()");
    return 0;
  }


  virtual void updateVal(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    Array2 <dcomplex> & newvals
    //!< The return: in form (MO)
  ) { error("Spinors do not need a spin index"); }

  virtual void updateVal(
    Sample_point * sample,
    int e,
    //!< electron number
    Array2 <dcomplex> & newvals
    //!< The return: in form (MO)
  );
  
  virtual void getBasisVal(
    Sample_point * sample,
    int e,
    Array1 <dcomplex> & newvals
    ){
    error("Need to implement MO_matrix_quaternion::getBasisVal()");
  }

  virtual void updateLap(
    Sample_point * sample,
    int e,
    int listnum,
    Array2 <dcomplex> & newvals
  ) { error("Spinors do not need a spin index"); }

  virtual void updateLap(
    Sample_point * sample,
    int e,
    Array2 <dcomplex> & newvals
  );

  virtual void updateSpinLap(
    Sample_point * sample,
    int e,
    //!< electron number
    Array2 <dcomplex> & newvals
    //!< The return: in form (MO)
  );
  virtual void updateHessian(Sample_point * sample,
			     int e,
			     int listnum,
			     Array2 <dcomplex>& newvals
			     //!< in form ([value gradient, dxx,dyy,dzz,dxy,dxz,dyz], MO)
			     ) { error("Spinors do not need a spin index"); }

  virtual void updateHessian(Sample_point * sample,
			     int e,
			     Array2 <dcomplex>& newvals
			     //!< in form ([value gradient, dxx,dyy,dzz,dxy,dxz,dyz], MO)
			     );
  Spinor_matrix_quaternion()
  {}

};

//######################################################################

#include "Qmc_std.h"
#include "Spinor_matrix_quaternion.h"
#include "Sample_point.h"
#include "qmc_io.h"



void Spinor_matrix_quaternion::init() {

  
  //Determine where to cut off the basis functions
  
  cutoff.Resize(totbasis);
  int basiscounter=0;
  for(int i=0; i< centers.size(); i++)
  {
    for(int j=0; j< centers.nbasis(i); j++)
    {
      Basis_function* tempbasis=basis(centers.basis(i,j));
      for(int n=0; n< tempbasis->nfunc(); n++)
      {
        //cout << "cutoff " << endl;
        cutoff(basiscounter)=tempbasis->cutoff(n);
        //cout << "rcut(" << basiscounter << ") "
        //     << cutoff(basiscounter) << endl;
        basiscounter++;
      }
    }
  }

  obj_cutoff.Resize(basis.GetDim(0));
  for(int b=0; b< basis.GetDim(0); b++) {
    int nf=basis(b)->nfunc();
    doublevar maxcut=basis(b)->cutoff(0);
    for(int n=1; n< nf; n++) {
      doublevar cut=basis(b)->cutoff(n);
      if(cut > maxcut) maxcut=cut;
    }
    obj_cutoff(b)=maxcut;
  }
  
  
  nfunctions.Resize(basis.GetDim(0));
  for(int b=0; b< basis.GetDim(0); b++) {
    nfunctions(b)=basis(b)->nfunc();
  }
  

  nbasis.Resize(nmo);
  mofill.Resize(nmo, totbasis);
  moCoeff2.Resize(nmo, totbasis);

  ifstream ORB(orbfile.c_str());

  if(!ORB) {
    error("couldn't find orb file ", orbfile);
  }

  Array3 <int> coeffmat;
  Array1 <quaternion> coeff;
  
  readorb(ORB,centers, nmo, maxbasis,kpoint, coeffmat, coeff);
  string in;
  ORB.close();

  
  //Find the cutoffs

  int totfunc=0;
  nbasis=0;
  //basismo=0;
  const doublevar threshold=1e-12;
  for(int ion=0; ion<centers.size(); ion++)
  {
    int f=0;

    doublevar dot=0;
    for(int d=0; d<3; d++) dot+=centers.centers_displacement(ion,d)*kpoint(d);
    
    //T kptfac=T(exp(dcomplex(0,1.0)*dot*pi));
    dcomplex kptfac=eval_kpoint_fac<dcomplex>(dot);
    
    for(int n=0; n< centers.nbasis(ion); n++) {
      
      int fnum=centers.basis(ion,n);
      int imax=basis(fnum)->nfunc();

      for(int i=0; i<imax; i++) { //sum over the symmetries
        for(int mo=0; mo<nmo; mo++) {      //and the MO's
          quaternion temp;
          if(coeffmat(mo,ion, f) == -1) {
            temp = quaternion();
            //cout << "missing MO pointer: mo# " << mo << " ion # " << ion
            //<< " function on ion: " << f << endl;
            //error("In the orb file, there is a missing pointer. It might "
            //      "be a badly structured file.");
          }
          else temp=coeff(coeffmat(mo,ion,f));
          if(abs(temp) > threshold) {
            mofill(mo, nbasis(mo))=totfunc;
            moCoeff2(mo, nbasis(mo))=kptfac*magnification_factor*temp;
            nbasis(mo)++;
          }

        }//mo
        f++;  //keep a total of functions on center
        totfunc++;
      } //i
    } //n
  }  //ion
  symmvals_temp1d.Resize(maxbasis);
  symmvals_temp2d.Resize(maxbasis,10);

}

//---------------------------------------------------------------------------------------------

void Spinor_matrix_quaternion::writeorb(ostream & os, 
    Array2 <doublevar> & rotation, Array1 <int>  &moList) {


  int nmo_write=moList.GetDim(0);
  assert(rotation.GetDim(0)==nmo_write);
  assert(rotation.GetDim(1)==nmo_write);

  os.precision(15);
  int counter=0;
  int totfuncs=0;
  for(int ion=0; ion<centers.equiv_centers.GetDim(0); ion++) {
    for(int n=0; n< centers.nbasis(centers.equiv_centers(ion,0)); n++) {
      int fnum=centers.basis(centers.equiv_centers(ion,0),n);
      totfuncs+=basis(fnum)->nfunc();
    }
  }
  
  for(int m=0; m < nmo_write; m++) {
    for(int ion=0; ion<centers.equiv_centers.GetDim(0); ion++) {
      int f=0;
      for(int n=0; n< centers.nbasis(centers.equiv_centers(ion,0)); n++) {
        int fnum=centers.basis(centers.equiv_centers(ion,0),n);
        int imax=basis(fnum)->nfunc();

        for(int i=0; i<imax; i++) {
          os << m+1 << "  "   << f+1 << "   " << ion+1 << "   " << counter+1 << endl;
          f++;  //keep a total of functions on center
          counter++;
        } //i
      } //n
    }  //ion
  }
  os << "COEFFICIENTS\n";
  ifstream orbin(orbfile.c_str());

  //cout << "orbfile " << orbfile << endl;
  Array1 <quaternion> coeff;
  Array3 <int> coeffmat;
  readorb(orbin,centers, nmo, maxbasis,kpoint, coeffmat, coeff);
  orbin.close();

  Array2 <quaternion> moCoeff(nmo,totbasis);
  for (int i =0; i < nmo; i++) {
      for (int j = 0; j < totbasis; j++) 
	  moCoeff(i,j) = quaternion();
  }
  int nmo_read=coeffmat.GetDim(0);
  int ncenter=coeffmat.GetDim(1);
  int maxfunc=coeffmat.GetDim(2);
  int currfunc;
  for(int mo=0; mo < nmo_read; mo++) { 
    currfunc=0;
    for(int ion=0; ion<centers.equiv_centers.GetDim(0); ion++) {
      int f=0;
      int equiv_center=centers.equiv_centers(ion,0);
      for(int n=0; n< centers.nbasis(centers.equiv_centers(ion,0)); n++) {
        int fnum=centers.basis(centers.equiv_centers(ion,0),n);
        int imax=basis(fnum)->nfunc();
        for(int i=0; i<imax; i++) {
          if(coeffmat(mo,equiv_center,f)!=-1)
            moCoeff(mo,currfunc)=coeff(coeffmat(mo,equiv_center,f));
          f++;
          currfunc++;
        }
      }
    }
  }
  // Now we rotate and write out the coefficients
  Array2 <quaternion> rotatedMO(nmo_write,currfunc);
  for (int i =0; i < nmo; i++) {
      for (int j = 0; j < totbasis; j++) 
	  rotatedMO(i,j) = quaternion();
  }
  for(int mo=0; mo < nmo_write; mo++) { 
    for(int f=0; f< currfunc; f++) { 
      for(int mo2=0; mo2 < nmo_write; mo2++) { 
        int realmo=moList(mo2);
        rotatedMO(mo,f)+=rotation(mo,mo2)*moCoeff(realmo,f);
        //cout << "mo1 " << mo << " mo2 " << mo2 << " f " << f 
        //  << " realmo " << realmo << " rotation " << rotation(mo,mo2) 
        //     << " rotatedmo " << rotatedMO(mo,f) << " coeff " << moCoeff(realmo,f) << endl;
      }
    }
  }
  int counter2=1;
  for(int m=0; m < nmo_write; m++) {
    for(int f=0; f< currfunc; f++) {
      os << rotatedMO(m, f) << "   ";
      if(counter2 % 5 ==0) os << endl;
      counter2++;
    }
  }
}
//---------------------------------------------------------------------

void Spinor_matrix_quaternion::buildLists(Array1 <int> & occupations){
  int nmo_list=occupations.GetDim(0);
  basisfill_list.Resize(totbasis, nmo_list);
  moCoeff_list.Resize(totbasis, nmo_list);
  basismo_list.Resize(totbasis);
  basismo_list=0;
  for(int i=0; i < nmo_list; i++)
  {
    int mo=occupations(i);
    for(int bas=0; bas < nbasis(mo); bas++)
    {
      int func=mofill(mo, bas);

      //basisfill_list(lis)(func, basismo_list(lis)(func))=mo;
      basisfill_list(func, basismo_list(func))=i;
      moCoeff_list(func, basismo_list(func))=moCoeff2(mo, bas);
      //cout << "basisfill_list " << 2 << "  f  "
      //     << func << "  mo " <<  mo;
      //cout << "  basis coeff " << moCoeff2(mo, bas) << endl;
      //cout << "real basisfill " << basisfill(func, basismo_list(lis)(func))
      //     << " coeff  " << moCoeff(func, basismo_list(lis)(func)) << endl;
      basismo_list(func)++;
    }
  }
}


//----------------------------------------------------------------------

int Spinor_matrix_quaternion::showinfo(ostream & os)
{
  os << "Cutoff MO " << endl;
  os << "Number of molecular orbitals: " << nmo << endl;
  string indent="  ";
  os << "Basis functions: \n";
  for(int i=0; i< basis.GetDim(0); i++)
  {
    basis(i)->showinfo(indent, os);
  }
  return 1;
}

int Spinor_matrix_quaternion::writeinput(string & indent, ostream & os)
{
  os << indent << "QUATERNION_MO" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "ORBFILE " << orbfile << endl;
  //if(oldsofile!="") 
  //  os << indent << "OLDSOFILE " << oldsofile << endl;
  os << indent << "MAGNIFY " << magnification_factor << endl;
  string indent2=indent+"  ";
  for(int i=0; i< basis.GetDim(0); i++)
  {
    os << indent << "BASIS { " << endl;
    basis(i)->writeinput(indent2, os);
    os << indent << "}" << endl;
  }

  os << indent << "CENTERS { " << endl;
  centers.writeinput(indent2, os);
  os << indent << "}" << endl;
  return 1;
}
//------------------------------------------------------------------------

void Spinor_matrix_quaternion::updateVal(
  Sample_point * sample,  int e,  Array2 <dcomplex> & newvals) {
  //cout << "start updateval " << endl;
  int centermax=centers.size();

  static Array1 <doublevar> R(5);
  doublevar spin;
  sample->getElectronSpin(e,spin);

  //Array1 <doublevar> symmvals_temp(maxbasis);

  //Make references for easier access to the list variables.
  Array1 <int> & basismotmp(basismo_list);
  Array2 <int> & basisfilltmp(basisfill_list);
  Array2 <quaternion> & moCoefftmp(moCoeff_list);
  assert(newvals.GetDim(1) >= 1);

  newvals=0;
  Basis_function * tempbasis;

  //int fn;
  quaternion c;
  int mo=0;
  int scalebasis=basisfill_list.GetDim(1);
  int totfunc=0;
  int b; //basis
  //cout << "here " << endl;
  centers.updateDistance(e, sample);
  //int retscale=newvals.GetDim(1);
  for(int ion=0; ion < centermax; ion++) {
    //sample->getECDist(e, ion, R);
    centers.getDistance(e,ion,R);
    for(int n=0; n< centers.nbasis(ion); n++) {
      b=centers.basis(ion,n);
      tempbasis=basis(b);
      if(obj_cutoff(b) > R(0)) {
        tempbasis->calcVal(R, symmvals_temp1d);
        //cout << "ion " << ion << "b " << b << " mo "<< mo << endl;
        int imax=nfunctions(b);
        for(int i=0; i< imax; i++) {
          int reducedbasis=scalebasis*totfunc;
          if(R(0) < cutoff(totfunc)) {
            for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++) {
              //mo=basisfill(totfunc, basmo);
              //c=moCoeff(totfunc, basmo);
              //cout << "basisfill reducedbasis "<< reducedbasis
              // << "  basmo " << basmo << endl;
              mo=basisfilltmp.v[reducedbasis+basmo];
              //cout << "mocoeff (mo=" << mo <<  endl;
              //mo_counter(mo)++;
              c=moCoefftmp.v[reducedbasis+basmo];

              newvals(mo, 0)+=(c.val.first*exp(I*spin)+c.val.second*exp(-I*spin))*symmvals_temp1d(i);
              //newvals.v[retscale*mo]+=c*symmvals_temp.v[i];

            }
          }
          totfunc++;
        }
      }
      else {
        totfunc+=nfunctions(b);
      }
    }
  }
  //n_calls++;
    //cout << "done updateVal " << endl;
}

//------------------------------------------------------------------------

void Spinor_matrix_quaternion::updateLap( Sample_point * sample,
  int e, Array2 <dcomplex> & newvals) {

  //cout << "updateLap" << endl;
  int centermax=centers.size();
  //int momax=occupation.GetDim(0);
  newvals=0;
  //assert(momax <= nmo);
  assert(e < sample->electronSize());

  assert(newvals.GetDim(1) >=5);

  // cout << "array " << endl;
  Array1 <doublevar> R(5);
  doublevar spin;
  sample->getElectronSpin(e,spin);
  // cout << "symvals " << endl;
   //cout << "arrayref " << endl;

  //References to make the code easier to read and slightly faster.
  Array1 <int> & basismotmp(basismo_list);
  Array2 <int> & basisfilltmp(basisfill_list);
  Array2 <quaternion> & moCoefftmp(moCoeff_list);

  Basis_function * tempbasis;

  quaternion c;
  int scaleval=0, scalesymm=0;
  int mo=0;
  int scalebasis=basisfilltmp.GetDim(1);
  centers.updateDistance(e, sample);
  int totfunc=0;
  int b;
  int symmvals_stride=symmvals_temp2d.GetDim(1);
  for(int ion=0; ion < centermax; ion++) {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++) {
      b=centers.basis(ion, n);
      tempbasis=basis(b);
      if(R(0) < obj_cutoff(b)) {
       // cout << "basis " << endl;
      tempbasis->calcLap(R, symmvals_temp2d);

      int imax=nfunctions(b);
      for(int i=0; i< imax; i++) {
        //cout << "i " << i << endl;

        int reducedbasis=scalebasis*totfunc;
        scalesymm=i*symmvals_stride;
        if(R(0) < cutoff(totfunc)) {
          for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++) {
            //mo=basisfill(basis, basmo);
            //c=moCoeff(basis, basmo);

            mo=basisfilltmp.v[reducedbasis+basmo];
            c=moCoefftmp.v[reducedbasis+basmo];
            //cout << c << "   ";

            //mo_counter(mo)++;

            scaleval=mo*5;
            //cout << "coeff " << c << endl;
            //cout << "reducedbasis " << reducedbasis << endl;
            // cout << mo << "   " << basmo << "   " << c << endl;


            //for(int j=0; j< 5; j++) {
            for(int j=0; j< 5; j++) {
              //newvals(mo,j)+=c*symmvals(fn,j);
              newvals.v[scaleval+j]+=(c.val.first*exp(I*spin)+c.val.second*exp(-I*spin))*symmvals_temp2d.v[scalesymm+j];
              //cout << "newvals(" << mo << "," << j << ")  " << newvals(mo,j) << endl;
              //cout << "symmvals(" << j << ")  " << symmvals(fn,j) << endl;
            }
	    //CM
	    //use scalesymm in symmvals_temp2d in order to just get spatial value. 
	    //newvals.v[scaleval+5]+=(I*c.val.first*exp(I*R(5))-I*c.val.second*exp(-I*R(5)))*symmvals_temp2d.v[scalesymm]; //ds
          }
          //cout << endl;

        }

        totfunc++;
      }
    }
    else {
      totfunc+=nfunctions(b);
    }
    }
  }


  //n_calls++;
  
  //cout << "newvals " << endl;
  //output_array(newvals);

}

void Spinor_matrix_quaternion::updateSpinLap(
  Sample_point * sample,  int e, Array2 <dcomplex> & newvals) {
  //cout << "start updateval " << endl;
  int centermax=centers.size();

  static Array1 <doublevar> R(5);
  doublevar spin;
  sample->getElectronSpin(e,spin);

  //Array1 <doublevar> symmvals_temp(maxbasis);

  //Make references for easier access to the list variables.
  Array1 <int> & basismotmp(basismo_list);
  Array2 <int> & basisfilltmp(basisfill_list);
  Array2 <quaternion> & moCoefftmp(moCoeff_list);
  assert(newvals.GetDim(1) >= 3);

  newvals=0;
  Basis_function * tempbasis;

  //int fn;
  quaternion c;
  int mo=0;
  int scalebasis=basisfill_list.GetDim(1);
  int totfunc=0;
  int b; //basis
  //cout << "here " << endl;
  centers.updateDistance(e, sample);
  //int retscale=newvals.GetDim(1);
  for(int ion=0; ion < centermax; ion++) {
    //sample->getECDist(e, ion, R);
    centers.getDistance(e,ion,R);
    for(int n=0; n< centers.nbasis(ion); n++) {
      b=centers.basis(ion,n);
      tempbasis=basis(b);
      if(obj_cutoff(b) > R(0)) {
        tempbasis->calcVal(R, symmvals_temp1d);
        //cout << "ion " << ion << "b " << b << " mo "<< mo << endl;
        int imax=nfunctions(b);
        for(int i=0; i< imax; i++) {
          int reducedbasis=scalebasis*totfunc;
          if(R(0) < cutoff(totfunc)) {
            for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++) {
              //mo=basisfill(totfunc, basmo);
              //c=moCoeff(totfunc, basmo);
              //cout << "basisfill reducedbasis "<< reducedbasis
              // << "  basmo " << basmo << endl;
              mo=basisfilltmp.v[reducedbasis+basmo];
              //cout << "mocoeff (mo=" << mo <<  endl;
              //mo_counter(mo)++;
              c=moCoefftmp.v[reducedbasis+basmo];

              newvals(mo, 0)+=(c.val.first*exp(I*spin)+c.val.second*exp(-I*spin))*symmvals_temp1d(i);
              newvals(mo, 1)+=(I*c.val.first*exp(I*spin)-I*c.val.second*exp(-I*spin))*symmvals_temp1d(i);
              newvals(mo, 2)+=(-c.val.first*exp(I*spin)-c.val.second*exp(-I*spin))*symmvals_temp1d(i);
              //newvals.v[retscale*mo]+=c*symmvals_temp.v[i];

            }
          }
          totfunc++;
        }
      }
      else {
        totfunc+=nfunctions(b);
      }
    }
  }
  //n_calls++;
    //cout << "done updateVal " << endl;
}
//--------------------------------------------------------------------------

void Spinor_matrix_quaternion::updateHessian(
  Sample_point * sample,
  int e,
  //const Array1 <int> & occupation,
  //!<A list of the MO's to evaluate
  Array2 <dcomplex> & newvals
  //!< The return: in form (MO, [val, grad, dxx,dyy,...])
)
{

  int centermax=centers.size();
  newvals=0;
  assert(e < sample->electronSize());
  assert(newvals.GetDim(1)==10);
  

  Array1 <doublevar> R(5);
  doublevar spin;
  sample->getElectronSpin(e,spin);
 // static Array2 <doublevar> symmvals_temp(maxbasis,10);

  //References to make the code easier to read and slightly faster.
  Array1 <int> & basismotmp(basismo_list);
  Array2 <int> & basisfilltmp(basisfill_list);
  Array2 <quaternion> & moCoefftmp(moCoeff_list);

  Basis_function * tempbasis;

  quaternion c;
  int scaleval=0, scalesymm=0;
  int mo=0;
  int scalebasis=basisfilltmp.GetDim(1);
  centers.updateDistance(e, sample);
  int totfunc=0;
  int b;
  int symmvals_stride=symmvals_temp2d.GetDim(1);
  for(int ion=0; ion < centermax; ion++)  {
    centers.getDistance(e, ion, R);
    for(int n=0; n< centers.nbasis(ion); n++)  {
      b=centers.basis(ion, n);
      tempbasis=basis(b);
      if(R(0) < obj_cutoff(b)) {
        tempbasis->calcHessian(R, symmvals_temp2d);

        int imax=nfunctions(b);
        for(int i=0; i< imax; i++)  {
          int reducedbasis=scalebasis*totfunc;
	  //This is spatial derivitive counter of basis, leave at 10
          scalesymm=i*10;
          if(R(0) < cutoff(totfunc))  {
            for(int basmo=0; basmo < basismotmp.v[totfunc]; basmo++)   {
              mo=basisfilltmp.v[reducedbasis+basmo];
              c=moCoefftmp.v[reducedbasis+basmo];
              scaleval=mo*symmvals_stride;
              for(int j=0; j< 10; j++) {
                newvals.v[scaleval+j]+=(c.val.first*exp(I*spin)+c.val.second*exp(-I*spin))*symmvals_temp2d.v[scalesymm+j];
              }
	      //CM
	      //spin derivitives. Don't add anything in symmvals in order to get just spatial value
	      //newvals.v[scaleval+10]+=(I*c.val.first*exp(I*R(5))-I*c.val.second*exp(-I*R(5)))*symmvals_temp2d.v[scalesymm]; //ds
	      //newvals.v[scaleval+11]-=(c.val.first*exp(I*R(5))+c.val.second*exp(-I*R(5)))*symmvals_temp2d.v[scalesymm]; //dss
            }
          }

          totfunc++;
        }
      }
      else {
        totfunc+=nfunctions(b);
      }
    }
  }
}

//--------------------------------------------------------------------------



#endif // SPINOR_MATRIX_QUATERNION_H_INCLUDED

