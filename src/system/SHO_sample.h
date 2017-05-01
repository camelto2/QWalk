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


#ifndef SHO_SAMPLE_H_INCLUDED
#define SHO_SAMPLE_H_INCLUDED


#include "Sample_point.h"
class Wavefunction;
#include "SHO_system.h"
class Sample_storage;


class SHO_sample : public Sample_point
{
public:

  ~SHO_sample()
  {}


  void init(System * sys);

  void randomGuess();

  int electronSize()
  {
    return nelectrons;
  }
  int ionSize()
  {
    return 0;
  }
  int centerSize()
  {
    return 0;
  }
  doublevar getIonCharge(const int ion)
  {
    return 0;
  }


  int ndim() { return nd; } 

  /*!
    Moves electron e
    position should be an array of length at least 3

  */
  void setElectronPos(const int e,const Array1 <doublevar> & position);

  void getElectronPos(const int e, Array1 <doublevar> & R)
  {
    assert( R.GetDim(0) >= 3 );
    for(int i=0; i< 3; i++)
    {
      //R(i) = points.r(i,e);
      R(i)=elecpos(e,i);
    }
  }

  void moveIon(const int ion, const Array1 <doublevar> & r);


  void getIonPos(const int ion, Array1 <doublevar> & r) {
    error("getIonPos not supported yet");
  }



  void updateEIDist();
  void updateEEDist();
  void updateECDist() {
    updateEIDist();
  }




  void getECDist(const int e, const int cent,
                 Array1 <doublevar> & distance)
  {
    getEIDist(e,cent, distance);
  }

  void getEIDist(const int e,const int ion, Array1 <doublevar> & distance)
  {
    error("getEIDist not supported yet");
  }
  void getEEDist(const int e1,const int e2, Array1 <doublevar> & distance)
  {
    //cout << "getEEDist" << endl;
    assert(distance.GetDim(0) >= 5);
    assert( ! elecDistStale(e1));
    assert( e1 < e2 );
    for(int i=0; i< 5; i++)
    {
      distance(i)=pointdist(i,e1,e2);
    }
  }

  void rawOutput(ostream &);
  void rawInput(istream &);

  void generateStorage(Sample_storage *& );
  void saveUpdate(int, Sample_storage * );
  void restoreUpdate(int, Sample_storage *);

  //CM:
  //Dynamic Spin functions below
  void getElectronSpin(const int e, doublevar & s) {
    error("getElectronSpin not implemented for SHO sample");
  }

  void setElectronSpin(const int e, const doublevar & s) {
    error("setElectronSpin not implemented for SHO sample");
  }

  void translateSpin(const int e, const doublevar & trans) {
    error("translateSpin not implemented for SHO sample");
  }

private:

  int nelectrons;

  Array2 <doublevar> elecpos; //electron positions
  Array1 <int> elecDistStale;
  Array3 <doublevar> pointdist;  //this is an upper triangular matrix;
                                 //the lower is currently wasted
  SHO_system * parent;
  int nd;
  
};



#endif //RING_SAMPLE_H_INCLUDED
//-------------------------------------------------------------------------
