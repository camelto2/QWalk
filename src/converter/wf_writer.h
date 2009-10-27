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

#ifndef WF_WRITER_H_INCLUDED
#define WF_WRITER_H_INCLUDED
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include "converter.h"
#include "basis_writer.h"

class Wf_writer {
  public:
    virtual void print_wavefunction(std::ostream & os)=0;
};


class Slat_wf_writer:public Wf_writer {
public:
  Slat_wf_writer() { magnification=-1; orbtype="ORBITALS"; } 
  int nup, ndown; //number of up and down electrons
  int spin_dwn_start;  //the MO at which the down orbitals start(for a UHF calctype only)
  double magnification;
  std::string calctype; //RHF, ROHF, UHF, or GVB
  std::string mo_matrix_type; //CUTOFF_MO, BLAS_MO, STANDARD_MO, etc
  std::string orbname;  //name of the orbitals file
  std::string orbtype;  //ORBITALS or CORBITALS
  std::string centername; //name of the centers file(if used, probably not)
  std::string basisname; //name of the basis file
  bool write_centers; //whether or not to use the centers file
  bool use_global_centers; //whether to use the global centers generated by the System object
                          //(for a molecular system, it doesn't matter, but for a periodic 
                          // system, this should be set to true, assuming we're using an atom-centered basis)
  
  std::vector < double > kpoint; //the k-point at which we're doing the calculation
                              //This really should be kept in the system, so don't use it.
                              //It's just kept around because it pervades crystal2qmc..


  std::vector < double > detwt; //the determinantal weights.  If this is not set, a single determinant will
                                //be assumed, and occ_up and occ_down will automatically be filled in.
                                //If detwt is set, then it is expected that occ-up and occ_down are filled in
                                //appropriately for each determinant. 
  std::vector < std:: vector < int > > occ_up;  //the occupation of the up molecular orbitals for each determinant,
                                                //starting from 1.
                                                //For single determinant calculations (RHF, UHF, ROHF), this will
                                                //automatically be filled as 1-nup.  
  std::vector < std:: vector < int > > occ_down;//same as occ_up, except it will be filled as 1-ndown for RHF and ROHF
                                                // and as spin_dwn_start-spin_down_start+ndown for UHF
  

  void print_wavefunction(std::ostream & inputfile );
};

class Atom;
class Basis_writer;

class Jastrow_wf_writer:public Wf_writer {
public:
  std::vector <std::string > atomnames;
  std::vector <Basis_writer * > basis;
  void print_wavefunction(std::ostream & inputfile);

  /*!
    Gets the unique names from the list of atoms.
  */
  void set_atoms(std::vector <Atom > & atoms);

  /*!
    Add a basis to the list of ones to print in the Jastrow section
  */
  void add_basis(Basis_writer & b);
};

class Jastrow2_wf_writer:public Wf_writer {
public:
  int ngroups;
  std::vector < int > spindiff; //!< Whether we have different like/unlike channels for a group


  std::vector <std::string > atomnames;
  std::vector < std::vector <Basis_writer * >  > ee_basis;
  std::vector < std::vector <Basis_writer * >  > ei_basis;
  void print_wavefunction(std::ostream & inputfile);

  /*!
    Gets the unique names from the list of atoms.
  */
  void set_atoms(std::vector <Atom > & atoms);

  void set_groups(int ngroups_) {
    ee_basis.resize(ngroups_);
    ei_basis.resize(ngroups_);
    ngroups=ngroups_;
    spindiff.resize(ngroups);
    for(int g=0; g< ngroups; g++) spindiff[g]=0;
  }

  void set_spin_diff(int group){
    assert(group < ngroups);
    spindiff[group]=1;
  }



  /*!
    Add a basis to the list of ones to print in the Jastrow section
  */
  void add_ee_basis(Basis_writer & b, int group);
  void add_ei_basis(Basis_writer & b, int group);
};



void print_std_jastrow2(Jastrow2_wf_writer & jast2writer, std::ostream & os,
                        double basis_cutoff);

void print_3b_jastrow2(std::ostream & os,std::vector<std::string> & unique_atoms, double cutoff=7.5);
/*
 Doubles a periodic cell in the specified direction.  latvec, moCoeff, and the atoms are all changed
 to reflect the new cell.  Assumes atom-centered basis functions..
*/ 

using namespace std;

template <typename T>
void fold_kpoint(Slat_wf_writer & slwriter, 
                 std::vector <std::vector <double> > & latvec,
                 int dir,
                 std::vector <std::vector <T> > & moCoeff,
                 std::vector <Atom> & atoms) { 
  //Note: this won't work for UHF wavefunctions..

  vector <Atom> natoms;
  for(vector<Atom>::iterator i=atoms.begin(); i!=atoms.end(); i++) { 
    Atom tmp=*i;
    for(int d=0;d < 3; d++) { 
      tmp.pos[d]+=latvec[dir][d];
    }
    natoms.push_back(tmp);
  }
  atoms.insert(atoms.end(), natoms.begin(), natoms.end());

  vector <vector <T> > nmocoeff;
  vector<T> motmp;
  int count=1;
  for(typename vector<vector<T> >::iterator i=moCoeff.begin(); i!=moCoeff.end();
      i++) 
    { 
      motmp.clear();
      //gamma point
      for(typename vector<T>::iterator j=i->begin(); j!= i->end(); j++) { 
	motmp.push_back(*j);
      }
      for(typename vector<T>::iterator j=i->begin(); j!= i->end(); j++) { 
	motmp.push_back(*j);
      }
      nmocoeff.push_back(motmp);
      motmp.clear();
      cout << "count " << count << " "  << nmocoeff.size() << " ";
      //X point
      for(typename vector<T>::iterator j=i->begin(); j!= i->end(); j++) { 
	motmp.push_back(*j);
      }
      for(typename vector<T>::iterator j=i->begin(); j!= i->end(); j++) { 
	motmp.push_back(-*j);
      }
      nmocoeff.push_back(motmp);
      cout << " " << nmocoeff.size() << endl;
      count++;
    }
  moCoeff=nmocoeff;
  
  for(int d=0; d< 3; d++) {
    latvec[dir][d]*=2;
  }
  
}

//----------------------------------------------------------------------
#endif