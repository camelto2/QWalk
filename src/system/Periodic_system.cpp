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

#include "Periodic_system.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include "Periodic_sample.h"
#include "qmc_io.h"

void Periodic_system::notify(change_type change, int n)
{
  switch(change)
  {
  case sample_static:
    //pseudo.setStatic(1);
    break;
  default:
    cout << "WARNING: Periodic system got a signal that it doesn't know: "
    << change << endl;
  }
}


int Periodic_system::generateSample(Sample_point *& samptr)
{
  assert(samptr==NULL);
  samptr=new Periodic_sample;
  samptr->init(this);
  return 1;
}

int Periodic_system::showinfo(ostream & os)
{
  const int ndim=3;
  os << "Periodic system " << endl;
  os << "Lattice vectors:" << endl;
  for(int i=0; i< ndim; i++) {
    os << i << " : ";
    for(int j=0; j < ndim; j++) {
      os << latVec(i,j) << "   " ;
    }
    os << endl;
  }

  os << "reciprocal vectors(multiplied by 2*pi) : " << endl;
  for(int i=0; i< ndim; i++) {
    os << i << " : ";
    for(int j=0; j< ndim; j++) {
      os << 2*pi*recipLatVec(i,j) << "    ";
    }
    os << endl;
  }
  os << "total number of points in reciprocal ewald sum: "
  << ngpoints << endl;

  os << "Self e-i " << self_ei << endl;
  os << "Self e-e " << self_ee << endl;
  os << "Self i-i " << self_ii << endl;
  os << "xc correction " << xc_correction <<" Warning: not added to E_total!"<< endl;
  os << "ion sum of ewald " << ion_ewald << endl;
  os << endl;
  ions.showinfo(os);

  os << "ionic polarization ";
  for(int d=0; d< 3; d++) os << ion_polarization(d) << "  ";
  os << endl;
  return 1;
}


//----------------------------------------------------------------------
#include "MatrixAlgebra.h"

void getCross(Array2 <doublevar> & latVec, Array2<doublevar> & crossProduct) { 
  crossProduct.Resize(3,3);
  assert(latVec.GetDim(0)>=3 && latVec.GetDim(1) >=3) ;
  crossProduct(0,0)=(latVec(1,1)*latVec(2,2)-latVec(1,2)*latVec(2,1));
  crossProduct(0,1)=(latVec(1,2)*latVec(2,0)-latVec(1,0)*latVec(2,2));
  crossProduct(0,2)=(latVec(1,0)*latVec(2,1)-latVec(1,1)*latVec(2,0));
  
  crossProduct(1,0)=(latVec(2,1)*latVec(0,2)-latVec(2,2)*latVec(0,1));
  crossProduct(1,1)=(latVec(2,2)*latVec(0,0)-latVec(2,0)*latVec(0,2));
  crossProduct(1,2)=(latVec(2,0)*latVec(0,1)-latVec(2,1)*latVec(0,0));
  
  crossProduct(2,0)=(latVec(0,1)*latVec(1,2)-latVec(0,2)*latVec(1,1));
  crossProduct(2,1)=(latVec(0,2)*latVec(1,0)-latVec(0,0)*latVec(1,2));
  crossProduct(2,2)=(latVec(0,0)*latVec(1,1)-latVec(0,1)*latVec(1,0));
  
}
/*!
*/
int Periodic_system::read(vector <string> & words,
                          unsigned int & pos)
{
  const int ndim=3;
  int startpos=pos;

  vector <string> latvectxt;

  vector <string> spintxt;
  if(!readsection(words, pos, spintxt, "NSPIN")) {
    error("Need NSPIN in system");
  }
  nspin.Resize(2);
  nspin(0)=atoi(spintxt[0].c_str());
  nspin(1)=atoi(spintxt[1].c_str());


  totnelectrons=nspin(0)+nspin(1);
  pos=startpos;
  ions.read(words, pos);

  int natoms=ions.size();
  for(int i=0; i< natoms; i++) {
    atomLabels.push_back(ions.getLabel(i));
  }


  vector <string> ktxt;
   
  if(readsection(words, pos=0, ktxt, "KPOINT")) {
    if(ktxt.size()!=3) error("KPOINT must be a section of size 3");
    kpt.Resize(3);
    for(int i=0; i< 3; i++) 
      kpt(i)=atof(ktxt[i].c_str());
  }
  else {
    kpt.Resize(3);
    kpt=0;
  }


  pos=startpos;
  if(!readsection(words, pos, latvectxt, "LATTICEVEC"))
    error("LATTICEVEC is required in PERIODIC");

  if(latvectxt.size() != ndim*ndim)
    error("LATTICEVEC must have exactly ",ndim*ndim, " values");
  latVec.Resize(ndim, ndim);

  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      latVec(i,j)=atof(latvectxt[i*ndim+j].c_str());
    }
  }
  //2D
  dim2 = 0;
  if (latVec(2,0) == 0 && latVec(2,1) == 0 && latVec(2,2) == 100) {
      single_write(cout,"Detected the 3rd lattice vector as (0,0,100)\n");
      single_write(cout,"Running as a 2D system with 2D ewald sums\n");
      dim2 = 1;
  }

  inverseLatVec.Resize(ndim,ndim);
  InvertMatrix(latVec,inverseLatVec,ndim);
  latvectxt.clear();

  origin.Resize(3);
  vector <string> origintxt;
  if(readsection(words, pos=0, origintxt, "ORIGIN")) {
    if(origintxt.size() < 3) error("ORIGIN section must have at least 3 elements.");
    for(int i=0; i< 3; i++) origin(i)=atof(origintxt[i].c_str());
  }
  else {
    origin=0;   //defaulting the origin to zero
  }


  //--Find the ghost centers for MO evaluation

  doublevar cutoff_divider;
  if(!readvalue(words, pos=0, cutoff_divider, "CUTOFF_DIVIDER"))
    cutoff_divider=2.000001;

  int ewald_gmax=200;
  readvalue(words, pos=0, ewald_gmax,"EWALD_GMAX");

  Array2 <doublevar> atompos(natoms, 3);
  for(int i=0; i< natoms; i++) {
    for(int d=0; d< 3; d++) {
      atompos(i,d)=ions.r(d,i);
    }
  }


  find_centers(origin, latVec, atompos, centerpos, equiv_atom, 
               center_displacement, cutoff_divider);
  assert(centerpos.GetDim(0)==equiv_atom.GetDim(0));


  
  int ncenters=centerpos.GetDim(0);

  ncenters_atom.Resize(natoms);
  ncenters_atom=0;
  for(int cen=0; cen < ncenters; cen++) {
    assert(equiv_atom(cen) < natoms);

    ncenters_atom(equiv_atom(cen))++;
  }


  int maxcent=0;
  for(int at=0; at< natoms; at++) {
    if(ncenters_atom(at) > maxcent) maxcent=ncenters_atom(at);
  }

  ncenters_atom=0;
  equiv_centers.Resize(natoms,maxcent);
  for(int cen=0; cen < ncenters; cen++) {
    int at=equiv_atom(cen);
    equiv_centers(at, ncenters_atom(at) )=cen;
    ncenters_atom(at)++;
  }
  
  //Test the center calculation..
  //
  //ofstream cenout("tmp.centers");
  //for(int cen=0; cen < ncenters; cen++) {
  //  int at=equiv_atom(cen);
  //  cenout << atomLabels[at] << "   ";
  //  for(int d=0; d< 3; d++)
  //    cenout << centerpos(cen,d) << "   ";
  //  cenout << endl;
  //}
  //cenout.close();


  //-------------cross products

  //cross product:  0->1x2, 1->2x0, 2->0x1
  Array2 <doublevar> crossProduct(ndim, ndim);
  getCross(latVec,crossProduct);


  //------------reciprocal cell
  recipLatVec.Resize(ndim, ndim);
  doublevar det=Determinant(latVec, ndim);

  if(det <0) det*=-1;
  debug_write(cout, "cell volume ", det,"\n");
  cellVolume=det;
  if (dim2 == 1) {
      cellVolume = 0;
      cellVolume += crossProduct(2,0)*crossProduct(2,0);
      cellVolume += crossProduct(2,1)*crossProduct(2,1);
      cellVolume += crossProduct(2,2)*crossProduct(2,2);
      cellVolume = sqrt(cellVolume);
      single_write(cout,"cell Area: ", cellVolume, "\n");
  }

  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      recipLatVec(i,j)=crossProduct(i,j)/det;
    }
  }

  //------------primitive lattice vectors (for polarization calculation)
  if(readsection(words, pos=0, latvectxt, "PRIMLATTICEVEC")) { 
    primlat.Resize(ndim,ndim);
    
    if(latvectxt.size() != ndim*ndim)
      error("LATTICEVEC must have exactly ",ndim*ndim, " values");
    
    for(int i=0; i< ndim; i++) {
      for(int j=0; j< ndim; j++) {
        primlat(i,j)=atof(latvectxt[i*ndim+j].c_str());
      }
    }
    Array2 <doublevar> cross(ndim, ndim);
    getCross(primlat,cross);
    
    
    
   //------------reciprocal cell
    prim_recip_vec.Resize(ndim, ndim);
    doublevar det=Determinant(primlat, ndim);
    
    debug_write(cout, "primitive cell volume ", det,"\n");
    
    for(int i=0; i< ndim; i++) {
      for(int j=0; j< ndim; j++) {
        prim_recip_vec(i,j)=cross(i,j)/det;
      }
    }
  }
  else {
    prim_recip_vec=recipLatVec;
    primlat=latVec;
  }

  //-------------------normal vectors

  normVec.Resize(ndim, ndim);

  for(int i=0; i< ndim; i++) {
    for(int j=0; j < ndim; j++) {
      normVec(i,j)=crossProduct(i,j);
    }

    //Check to make sure the direction is facing out
    doublevar dotprod=0;
    for(int j=0; j < ndim; j++) {
      dotprod+=normVec(i,j)*latVec(i,j);
    }
    if(dotprod < 0) {
      for(int j=0; j< ndim; j++) {
        normVec(i,j)= -normVec(i,j);
      }
    }
  }

  corners.Resize(ndim, ndim);
  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      corners(i,j)=origin(j)+latVec(i,j);
    }
  }

  
  
  //Make sure the ions are inside the simulation cell
  Array1 <int> nshift;
  for(int at=0; at< natoms; at++) {
    Array1 <doublevar> ionpos(3);
    for(int d=0; d< 3; d++) ionpos(d)=ions.r(d,at);
    if(enforcePbc(ionpos, nshift))  {
      single_write(cout,"**Warning** atom ", at ,
                   " is out of the simulation cell.  Putting it inside.\n");
      for(int d=0; d< 3; d++) ions.r(d,at)=ionpos(d);
    }
  }



  //----Set up reciprocal ewald sum
  smallestheight=1e99;
  for(int i=0; i< ndim; i++) {
    doublevar tempheight=0;
    doublevar length=0;
    for(int j=0; j< ndim; j++) {
      tempheight+=crossProduct(i,j)*latVec(i,j);
      length+=crossProduct(i,j)*crossProduct(i,j);
    }
    tempheight=fabs(tempheight)/sqrt(length);
    if(tempheight < smallestheight ) smallestheight=tempheight;
  }
  debug_write(cout, "elsu ", smallestheight,"\n");

  alpha=5.0/smallestheight; //Heuristic?  Stolen from Lubos's code.
  single_write(cout, "alpha ", alpha, "\n");

  debug_write(cout, "alpha ", alpha, "\n");

  ///--------------------------------------------
  // Finding the g-points for the Ewald sum
  //
 //  ngpoints=27*ewald_gmax*ewald_gmax*ewald_gmax ;
  ngpoints=0;
  for(int ig=0; ig <= ewald_gmax; ig++) {
    int jgmin=-ewald_gmax;
    if(ig==0) jgmin=0;
    for(int jg=jgmin; jg <= ewald_gmax; jg++) {
      int kgmin=-ewald_gmax;
      if(ig==0 && jg==0) kgmin=0;
      for(int kg=kgmin; kg <= ewald_gmax && dim2 == 0; kg++) {
        doublevar gsqrd=0;
        for(int i=0; i< ndim; i++) {
	  if (dim2 == 1) recipLatVec(2,i) = 0;
          doublevar tmp=2*pi*(ig*recipLatVec(0,i)
                                       +jg*recipLatVec(1,i)
                                       +kg*recipLatVec(2,i));
          gsqrd+=tmp*tmp;
        }
        if(4.0 * pi*exp(-gsqrd/(4*alpha*alpha))
                               /(cellVolume*gsqrd) > 1e-10) ngpoints++;
      }
      if (dim2 == 1) {
	  doublevar gmag = 0;
	  for (int i = 0; i < 2; i++) {
              doublevar tmp=2*pi*(ig*recipLatVec(0,i)
                                           +jg*recipLatVec(1,i));
              gmag+=tmp*tmp;
	  }
	  gmag = sqrt(gmag);
	  if (pi/cellVolume/gmag > 5.e-3) ngpoints++;
      }
    }
  }
  gpoint.Resize(ngpoints, 3);
  gweight.Resize(ngpoints);
  int currgpt=0;
  for(int ig=0; ig <= ewald_gmax; ig++) {
    int jgmin=-ewald_gmax;
    if(ig==0) jgmin=0;
    for(int jg=jgmin; jg <= ewald_gmax; jg++) {
      int kgmin=-ewald_gmax;
      if(ig==0 && jg==0) kgmin=0;
      for(int kg=kgmin; kg <= ewald_gmax && dim2 == 0; kg++) {
        for(int i=0; i< ndim; i++) {
          gpoint(currgpt, i)=2*pi*(ig*recipLatVec(0,i)
                                       +jg*recipLatVec(1,i)
                                       +kg*recipLatVec(2,i));
        }
        //cout << "there" << endl;
        doublevar gsqrd=0; //|g|^2
        for(int i=0; i< ndim; i++) {
          gsqrd+=gpoint(currgpt,i)*gpoint(currgpt,i);
        }

        if(gsqrd > 1e-8) {
          gweight(currgpt)=4.0 * pi*exp(-gsqrd/(4*alpha*alpha))
                               /(cellVolume*gsqrd);
          if(gweight(currgpt) > 1e-10) { //
            currgpt++;
          }
        }
      }
      if (dim2 == 1) {
        for(int i=0; i< 2; i++) {
          gpoint(currgpt, i)=2*pi*(ig*recipLatVec(0,i)
                                       +jg*recipLatVec(1,i));
        }
	gpoint(currgpt,2) = 0;
	doublevar gmag = 0;
	for (int i = 0; i < 2; i++) {
	    gmag += gpoint(currgpt,i)*gpoint(currgpt,i);
	}
	gmag = sqrt(gmag);
        if(gmag*gmag > 1e-8) {
	    gweight(currgpt) = pi/cellVolume/gmag;
	    if (gweight(currgpt) > 5.e-3) currgpt++;
        }
      }
    }
  }
  
  //Array2 <doublevar> gpointtemp(ngpoints,3);
  //Array1 <doublevar> gweighttemp(ngpoints);
  

  single_write(cout,"Ewald sum using ",ngpoints," reciprocal points\n");
  //Done finding the g-points.
  //---------------------------------------

  //Adjust to the correct number of kpoints..
  /*
  ngpoints=currgpt;
  gpoint.Resize(ngpoints, 3);
  gweight.Resize(ngpoints);



  for(int i=0; i< ngpoints; i++) {
    for(int d=0; d< ndim; d++) {
      gpoint(i,d)=gpointtemp(i,d);
    }
    gweight(i)=gweighttemp(i);
  }
*/
  //Resize the stored ion variables
  ion_sin.Resize(ngpoints);
  ion_cos.Resize(ngpoints);
  constEwald();
  ion_ewald=ewaldIon();

  vector < vector <string> > pseudotext;
  vector <string> pseudotexttmp;
  if(readsection(words, pos, pseudotexttmp, "PSEUDO") != 0) {
    error("pseudo section is now in the global space");
  }
  
  ion_polarization.Resize(3);
  ion_polarization=0.0;
  Array1 <doublevar> ion_pos(3);
  for(int at=0; at < natoms; at++) {
    getIonPos(at, ion_pos);
    for(int d=0; d< 3; d++) {
      ion_polarization(d)+=ion_pos(d)*ions.charge(at) ;
    }
  }


  if (dim2 == 1) {
      self_ei = ewaldSelf();
      ion_ewald = ewaldII();
  }

  return 1;
}


//----------------------------------------------------------------------

void Periodic_system::setIonPos(int ion, Array1 <doublevar> & r) {
  assert(r.GetDim(0) >= 3);


  Array1 <doublevar> temp(r.GetDim(0));
  temp=r;
  Array1<int> nshift;
  enforcePbc(temp, nshift);

  Array1 <doublevar> oldr(3);
  getIonPos(ion, oldr);
  Array1 <doublevar> deltar(3);
  for(int i=0; i< 3; i++) {
    deltar(i)=temp(i)-oldr(i);
  }

  for(int d=0; d< 3; d++) {
    ions.r(d,ion)=temp(d);
    for(int n=0; n< ncenters_atom(ion); n++) {
      centerpos(equiv_centers(ion, n), d)+=deltar(d);
    }
  }

  ion_ewald=ewaldIon();
}

//------------------------------------------------------------------------

doublevar Periodic_system::minDistance(Array1 <doublevar> pos1, Array1 <doublevar> pos2, Array1 <doublevar> &rmin ) {
  doublevar height2=smallestheight*smallestheight*.25;
  int nlat=primlat.GetDim(0);
  //static Array1 <doublevar> tmpvec(3),rmin(3);
  doublevar tmpvec[3]; 
  doublevar tmpdis=0,dismin=1e99;
  Array1<doublevar> r(3); 
  for (int d=0; d<3; d++) {
    r[d] = pos2[d] - pos1[d];
    rmin(d) = r[d]; 
    tmpdis+=r[d]*r[d];
  }
  if(tmpdis < height2) return tmpdis;
  
  for(int a=0; a < nlat; a++) { 
    tmpdis=0;
    for(int d=0; d < 3; d++) { 
      tmpvec[d]=r(d)+primlat(a,d);
      tmpdis+=tmpvec[d]*tmpvec[d];
    }
    if(tmpdis < dismin) { 
      for(int d=0; d< 3; d++) rmin[d]=tmpvec[d];
      dismin=tmpdis;
      //cout << "dismin " << dismin << endl;
      if(tmpdis< height2)
        break;
    }
    if(tmpdis < height2) break;
  }
  
  return dismin;
}

int Periodic_system::enforcePbc(Array1 <doublevar> & pos, Array1 <int> & nshifted) {
  assert(pos.GetDim(0) >=3);
  int shifted=0;
  nshifted.Resize(3);
  nshifted=0;
  //cout << "pos " << pos(0) << "  " << pos(1) << "  " << pos(2) << endl;
  // Array1 <doublevar> oldpos=pos;

  int nshift=0;
  for(int i=0; i< 3; i++) {
    int shouldcontinue=1;
    while(shouldcontinue) {

      //Whether we're past the origin side
      doublevar tooshort=0;
      for(int j=0; j<3;j++) tooshort+=normVec(i,j)*(pos(j)-origin(j));

      //Whether we're past the lattice vector
      doublevar toofar=0;
      for(int j=0; j< 3; j++) toofar+=normVec(i,j)*(pos(j)-corners(i,j));

      //the 1e-12 seems to help avoid numerical problems, esp 
      //when integrating over a grid(which tends to hit the edges)
      if(tooshort < -1e-12) {
	//cout <<"tooshort " << tooshort << endl;
        for(int j=0; j< 3; j++) pos(j)+=latVec(i,j);
        shifted=1;
        nshifted(i)+=1;
        nshift++;
      }
      else if(toofar > 1e-12) {
	//cout << "toofar " << toofar << endl;
        for(int j=0; j< 3; j++) pos(j)-=latVec(i,j);
        shifted=1;
	// JK: here was +1, which works fine for real k-points that
	// correspond to standing waves. For general (complex) k-points the
	// wavefunction is "directional", however.
        nshifted(i)-=1;
        nshift++;
      }
      else {
        shouldcontinue=0;
      }
      if(nshift > 1000)
	
        error("Did over 1000 shifts and we're still out of the simulation cell."
	            "  There's probably something wrong.  Position : ", pos(i));
    }
  }

  //if(nshift>0) { 
  //  cout << "shifted " << nshifted(0) << "  " << nshifted(1) << "  " << nshifted(2) <<endl;
  //  cout << "oldpos " << oldpos(0) << "  " << oldpos(1) << "  " << oldpos(2) << endl;
  //  cout << "newpos " << pos(0) << "  " << pos(1) << "  " << pos(2) << endl;
  //}

  return shifted;

}



//----------------------------------------------------------------------


doublevar Periodic_system::calcLoc(Sample_point * sample)
{
    if (dim2 == 0) { //3D
        doublevar ewalde=ewaldElectron(sample);
        //cout << "ion_ewald " << ion_ewald << " self_ii " << self_ii
        //     << " self_ee " << self_ee << " self_ei " << self_ei << endl;
        //cout << " ewalde " << ewalde << " xc_correction " << xc_correction << endl;
        //we do not want the xc_correction in the total energy in order to compare 
        //to all other qmc codes, it is still printed out so can be added by hand 
        
	/*
        // Fraser
        doublevar electronelectron = ewaldEE(sample);
        doublevar ionion = ewaldII(sample);
        doublevar electronion = ewaldEI(sample);
        doublevar self = ewaldSelf(sample);
    
    
        cout << endl;
        cout << "xc:           " << xc_correction << endl;
        cout << "QWalk + xc    " << ion_ewald+self_ii+self_ee+self_ei+ewalde+xc_correction << endl;
        cout << "QWalk Ewald:  " << ion_ewald+self_ii+self_ee+self_ei+ewalde << endl;
        cout << "Fraser Ewald: " << electronelectron+ionion+electronion+self << endl;
        cout << endl;
	*/
    
        return ion_ewald+self_ii+self_ee+self_ei+ewalde; //+xc_correction;
    }
    else {//2D
        doublevar electronelectron = ewaldEE(sample);
        doublevar electronion = ewaldEI(sample);
	return electronelectron+electronion+self_ei+ion_ewald; //2D EWald
    }

}

doublevar Periodic_system::psi(Array1 <doublevar> & pos1, Array1 <doublevar> & pos2) {

    Array1 <doublevar> rij(3);
    for (int d = 0; d < 3; d++) rij(d) = pos1(d) - pos2(d);

    doublevar recip = 0.0;
    for (int gpt = 0; gpt < ngpoints; gpt++) {
	doublevar gdotr = 0.0;
	for (int d = 0; d < 3; d++) gdotr += gpoint(gpt,d)*rij(d);
	recip += gweight(gpt)*cos(gdotr);
    }

    doublevar constpsi = -pi/alpha/alpha/cellVolume;

    doublevar real = 0.0;
    const int nlatvec=2;
    for(int kk=-nlatvec; kk <=nlatvec; kk++) {
      for(int jj=-nlatvec; jj <=nlatvec; jj++) {
        for(int ii=-nlatvec; ii <=nlatvec; ii++) {
	  Array1 <doublevar> pos(3);
	  pos = 0;
          for(int d=0; d< 3; d++) {
            pos(d) = rij(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
          }
          doublevar r = pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2);
	  r= sqrt(r);
          real+=erfcm(alpha*r)/r;
        }
      }
    }

    return 2.0*recip + constpsi + real;
    // multiply by two on recip for symmetry +G / -G
}

doublevar Periodic_system::psi2d(Array1 <doublevar> & pos1, Array1 <doublevar> & pos2) {

    Array1<doublevar> rij(3);
    for (int d = 0; d < 3; d++) rij(d) = pos1(d) - pos2(d);

    doublevar real = 0;
    const int nlatvec = 4;
    for (int ii = -nlatvec; ii <= nlatvec; ii++) {
	for (int jj = -nlatvec; jj <= nlatvec; jj++) {
	    Array1 <doublevar> pos(3);
	    pos = rij;
	    doublevar r = 0;
	    for (int d = 0; d < 2; d++) {
		pos(d) -= ii*latVec(0,d) + jj*latVec(1,d);
	    }
	    for (int d = 0; d < 3; d++) 
		r += pos(d)*pos(d);
	    r = sqrt(r);
	    real += erfcm(r*alpha)/r;
	}
    }

    doublevar constant = 0;
    constant = rij(2)*erf(rij(2)*alpha);
    constant += exp(-rij(2)*rij(2)*alpha*alpha)/alpha/sqrt(pi);
    constant *= -2*pi/cellVolume;

    doublevar recip = 0;
    for (int gpt = 0; gpt < ngpoints; gpt++) {
	doublevar gdotr = 0;
	doublevar gmag = 0;
	for (int d = 0; d < 2; d++) {
	    gdotr += gpoint(gpt,d)*rij(d);
	    gmag += gpoint(gpt,d)*gpoint(gpt,d);
	}
	gmag = sqrt(gmag);
	doublevar tmp;
	doublevar argp,argm;
	doublevar z = abs(rij(2));
	argp = gmag/2/alpha + z*alpha;
	argm = gmag/2/alpha - z*alpha;
	if (gmag*z < 500) {
	    tmp = exp(-gmag*z)*erfc(argm) + exp(gmag*z)*erfc(argp);
	}
	else {
	    tmp = exp(-gmag*z-argm*argm)/sqrt(pi)*(2*alpha/gmag+4*z*alpha*alpha*alpha/gmag/gmag);
	    tmp += exp(gmag*z-argp*argp)/sqrt(pi)*(2*alpha/gmag-4*z*alpha*alpha*alpha/gmag/gmag);
	}
	recip += gweight(gpt)*tmp*cos(gdotr);
    }
    return 2.0*recip + constant + real;
}

doublevar Periodic_system::zeta() {

    doublevar recip = 0.0;
    for (int gpt = 0; gpt < ngpoints; gpt++) recip += gweight(gpt);

    doublevar constzeta = 0.0;
    constzeta = -pi/alpha/alpha/cellVolume;
    constzeta -= 2.0*alpha/sqrt(pi);

    doublevar real = 0.0;
    const int nlatvec=2;
    for(int kk=-nlatvec; kk <=nlatvec; kk++) {
      for(int jj=-nlatvec; jj <=nlatvec; jj++) {
        for(int ii=-nlatvec; ii <=nlatvec; ii++) {
	  if ( ii != 0 && jj != 0 && kk != 0 ) {
	    Array1 <doublevar> Rvec(3);
	    Rvec = 0;
            for(int d=0; d< 3; d++) {
              Rvec(d) = kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
            }
            doublevar r = Rvec(0)*Rvec(0) + Rvec(1)*Rvec(1) + Rvec(2)*Rvec(2);
	    r = sqrt(r);
            real += erfcm(alpha*r)/r;
	  }
        }
      }
    }
    
    return 2.0*recip + constzeta + real;
    //G symmetry
}

doublevar Periodic_system::zeta2d() {

    doublevar real = 0;
    const int nlatvec = 4;
    for (int ii = -nlatvec; ii <= nlatvec; ii++) {
	for (int jj = -nlatvec; jj <= nlatvec; jj++) {
	    if ( ii != 0 && jj != 0) {
		Array1 <doublevar> Rvec(2);
		doublevar R = 0;
		Rvec = 0;
		for (int d = 0; d < 2; d++) {
		    Rvec(d) = ii*latVec(0,d) + jj*latVec(1,d);
		    R += Rvec(d)*Rvec(d);
		}
		R = sqrt(R);
		real += erfcm(alpha*R)/R;
	    }
	}
    }

    doublevar constant = 0;
    constant -= 2*alpha/sqrt(pi) + 2*sqrt(pi)/alpha/cellVolume;

    doublevar recip = 0;
    for (int gpt = 0; gpt < ngpoints; gpt ++) {
	doublevar gmag = 0;
	for (int d = 0; d < 2; d++) gmag += gpoint(gpt,d)*gpoint(gpt,d);
	gmag = sqrt(gmag);
	recip += 2.0*gweight(gpt)*erfc(gmag/alpha/2);
    }

    return 2.0*recip + constant + real;
}

doublevar Periodic_system::ewaldEE(Sample_point * sample) {
  sample->updateEEDist();
  sample->updateEIDist();

  Array1 <doublevar> e1pos(3);
  Array1 <doublevar> e2pos(3);
  doublevar sum = 0.0;
  for (int e1 = 0; e1 < totnelectrons; e1++) {
      for (int e2 = 0; e2 < totnelectrons; e2++) {
	  if ( e1 != e2 ) {
	     sample->getElectronPos(e1,e1pos);
	     sample->getElectronPos(e2,e2pos);
             sum += psi2d(e1pos,e2pos);
	  }
      }
  }
  sum /= 2.0;

  return sum;

}

doublevar Periodic_system::ewaldII() {

   Array1 <doublevar> i1pos(3);
   Array1 <doublevar> i2pos(3);
   doublevar sum = 0.0;
   for (int i1 = 0; i1 < ions.size(); i1++) {
       for (int i2 = 0; i2 < ions.size(); i2++) {
	   if ( i1 != i2 ) {
	       for (int d = 0; d < 3; d++) {
		   i1pos(d) = ions.r(d,i1);
		   i2pos(d) = ions.r(d,i2);
	       }
	       sum += ions.charge(i1)*ions.charge(i2)*psi2d(i1pos,i2pos);
	   }
       }
   }
   sum /= 2.0;
   return sum;
}

doublevar Periodic_system::ewaldEI(Sample_point * sample) {
  sample->updateEEDist();
  sample->updateEIDist();

    Array1 <doublevar> ipos(3);
    Array1 <doublevar> epos(3);
    doublevar sum = 0.0;
    for (int e = 0; e < totnelectrons; e++) {
	for (int i = 0; i < ions.size(); i++) {
            sample->getElectronPos(e,epos);
	    sample->getIonPos(i,ipos);
	    sum -= ions.charge(i)*psi2d(epos,ipos);
	}
    }

   return sum;
}

doublevar Periodic_system::ewaldSelf() {

    doublevar eself, iself;
    doublevar z = zeta2d();
    eself = 0.5*totnelectrons*z;

    doublevar chargesum = 0.0;
    for (int i = 0; i < ions.size(); i++) chargesum += ions.charge(i)*ions.charge(i);
    
    iself = 0.5*chargesum*z;

    return eself + iself;
}

//----------------------------------------------------------------------

void Periodic_system::calcLocSeparated(Sample_point * sample, Array1 <doublevar> & totalv)
{
  //  cout << "Not implemented" << endl; 
  //  doublevar ewalde=ewaldElectron(sample);
  
  //cout << "ion_ewald " << ion_ewald << " self_ii " << self_ii
  //     << " self_ee " << self_ee << " self_ei " << self_ei << endl;
  //cout << " ewalde " << ewalde << " xc_correction " << xc_correction << endl;
  //we do not want the xc_correction in the total energy in order to compare 
  //to all other qmc codes, it is still printed out so can be added by hand 
  for (int e=0; e<totnelectrons; e++) {
    totalv(e) = self_e_single + ewalde_separated(e); 
  }
    //  return ion_ewald+self_ii+self_ee+self_ei+ewalde; //+xc_correction;
}

//----------------------------------------------------------------------
/*!
\todo
Add in the correction terms for summing over non-zero l
 */

void Periodic_system::constEwald() {


  self_ee=0;
  self_ei=0;

  int nions=ions.size();

  doublevar ionIonSum=0; // The sums over i < j
  doublevar ionElecSum=0;
  doublevar elecElecSum=0;

  elecElecSum=totnelectrons*(totnelectrons-1);
  for(int ion=0; ion < nions; ion++) {
    ionElecSum+=ions.charge(ion)*(-totnelectrons);
    for(int ion2=ion+1; ion2 < nions; ion2++) {
      ionIonSum+=ions.charge(ion)*ions.charge(ion2);
    }
  }

  doublevar ionSum2=0; //The sums over the squares of q_i
  doublevar elecSum2=0;
  elecSum2=totnelectrons;
  for(int ion=0; ion < nions; ion++) {
    ionSum2+=ions.charge(ion)*ions.charge(ion);
  }

  doublevar squareconst=-.5*(2*alpha/sqrt(pi)+pi/(cellVolume*alpha*alpha));
  doublevar ijconst=-pi/(cellVolume*alpha*alpha);
  self_ii=ionIonSum*ijconst+ionSum2*squareconst;
  self_ei=ionElecSum*ijconst; // I guess this term is the charge+ and charge- background interactions
  self_ee=.5*elecElecSum*ijconst+elecSum2*squareconst;
  //                  e-e_i_bg                       e-Ion_bg                       e_e_self
  self_e_single = (totnelectrons - 1)*ijconst + ionElecSum/totnelectrons*ijconst + squareconst; 
  self_e_single_test= ionElecSum/totnelectrons*ijconst + squareconst; 
  ijbg = ijconst; 
  //Correct for the exchange-correlation false interaction with
  //the uniform electron gas number.  cexc was fitted to several
  //systems..
  doublevar rs=(3./(4*pi*totnelectrons/cellVolume));
  rs=pow(rs, 1.0/3.0);
  doublevar cexc=0.36;
  xc_correction=cexc/rs;

  //cout << "self ii " << self_ii << endl;
  //cout << "self ei " << self_ei << endl;
  //cout << "self ee " << self_ee << endl;

}


//----------------------------------------------------------------------

/*!

 */
doublevar Periodic_system::ewaldIon() {
  assert(ion_sin.GetDim(0) >= ngpoints);
  assert(ion_cos.GetDim(0) >= ngpoints);



  int nions=ions.size();
  const int nlatvec=2;  //Number of lattice vectors to iterate over

  Array1 <doublevar> r1(3), r2(3);
  doublevar IonIon=0;
  for(int i=0; i< nions; i++)
  {
    for(int j=0; j<i; j++)
    {
      for(int d=0; d< 3; d++) {
        r1(d)=ions.r(d,i)-ions.r(d,j);
      }


      //----over 2 lattice vectors
      for(int kk=-nlatvec; kk <=nlatvec; kk++) {
        for(int jj=-nlatvec; jj <=nlatvec; jj++) {
          for(int ii=-nlatvec; ii <=nlatvec; ii++) {
            for(int d=0; d< 3; d++) {
              r2(d)=r1(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
            }
            doublevar r=sqrt(r2(0)*r2(0)+r2(1)*r2(1)+r2(2)*r2(2));

            IonIon+=ions.charge(i)*ions.charge(j)*erfcm(alpha*r)/r;
//            cout << "r " << r << "  ionion " << IonIon << " i " << i << " j " << j << endl;
          }
        }
      }
      //----done lattice vectors
    }
  }


  debug_write(cout, "real space ion-ion ", IonIon,"\n");


  ion_sin=0;
  ion_cos=0;

  for(int gpt=0; gpt < ngpoints; gpt++) {

    for(int ion=0; ion <nions; ion++) {
      double dotprod=0;
      for(int d=0; d< 3; d++) dotprod+=ions.r(d,ion)*gpoint(gpt, d);
      //if(gpt==0) cout << "dot product " << dotprod << endl;
      ion_sin(gpt)+=ions.charge(ion)*sin(dotprod);
      ion_cos(gpt)+=ions.charge(ion)*cos(dotprod);
      //cout << "sine " << sin(dotprod) << "  cos " << cos(dotprod) << " charge " << ions.charge(ion)
      //     << "  total sin " << ion_sin(gpt) << "  ion_cos " << ion_cos(gpt)
      //     << endl;
    }
  }


  doublevar eion=0;
  for(int gpt=0; gpt < ngpoints; gpt++) {
    eion+=gweight(gpt)*(ion_sin(gpt)*ion_sin(gpt)+ion_cos(gpt)*ion_cos(gpt))/2;
    //cout << "gweight " << gweight(gpt) << " ion_sin " << ion_sin(gpt) << " ion_cos " << ion_cos(gpt)

      //   << endl;
  }

  eion*=2;

  debug_write(cout,"reciprocal space ion ", eion,"\n");
  //cout << "eion " << eion << " IonIon " << IonIon << endl;
  return eion+IonIon;
}

//----------------------------------------------------------------------


void Periodic_system::calcLocWithTestPos(Sample_point * sample, Array1 <doublevar> &tpos, Array1<doublevar> & Vtest) {
  sample->updateEEDist();
  sample->updateEIDist();
  Array1 <int> nshift; 
  if(enforcePbc(tpos, nshift))  { 
    cout <<"Warning, tpos is not in the cell" << endl; 
  }
  int nions=ions.size();
  int nlatvec=1;
  Array1 <doublevar> r1(3), r2(3), eidist(5), eedist(5);
  Vtest.Resize(totnelectrons+1);
  for (int e=0; e<totnelectrons; e++) {
    Vtest(e) = ijbg; 
  }
  Vtest(totnelectrons) = self_e_single_test; 
  
  Array1 <doublevar> rion(3); 
  Array1 <doublevar> oldpos(3); 
  //  sample->getElectronPos(0, oldpos);
  //sample->setElectronPosNoNotify(0, tpos);
  //sample->updateEIDist(); 
  for(int ion=0; ion < nions; ion++) {
    sample->getIonPos(ion, rion); 
    //sample->getEIDist(0, ion, eidist); 
    sample->minDist(tpos, rion, r1); //the minimum distance between ion and the test electron
    //----over  lattice vectors
    //for (int d=0; d<3; d++) {
    //  r1(d) = eidist(d+2); 
    //}
    for(int kk=-nlatvec; kk <=nlatvec; kk++) {
      for(int jj=-nlatvec; jj <=nlatvec; jj++) {
	for(int ii=-nlatvec; ii <=nlatvec; ii++) {
	  for(int d=0; d< 3; d++) {
	    r2(d)=r1(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
	  }
	  doublevar r=sqrt(r2(0)*r2(0)+r2(1)*r2(1)+r2(2)*r2(2));
	  Vtest(totnelectrons) -= ions.charge(ion)*erfcm(alpha*r)/r;
	  //            elecIon_real-=ions.charge(ion)*erfcm(alpha*r)/r;
	}
      }
    }
      //----done lattice vectors
  }

  //cout << "electron-electron " << endl;
  //-------------Electron-electron real part
  //sample->updateEEDist(); 
  
  Array2 <doublevar> elecpos(totnelectrons, 3);
  sample->getAllElectronPos(elecpos);
  for(int e =0; e < totnelectrons; e++) {
    //sample->getEEDist(0, e, eedist); 
    //for (int d=0; d<3; d++) r1(d)=eedist(d+2); 
    for (int d=0; d<3; d++) rion(d) = elecpos(e, d); 
    sample->minDist(tpos, rion, r1); 
    //----over  lattice vectors
    for(int kk=-nlatvec; kk <=nlatvec; kk++) {
      for(int jj=-nlatvec; jj <=nlatvec; jj++) {
	for(int ii=-nlatvec; ii <=nlatvec; ii++) {
	  for(int d=0; d< 3; d++) {
	    r2(d)=r1(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
	  }
	  doublevar r=sqrt(r2(0)*r2(0)+r2(1)*r2(1)+r2(2)*r2(2));
	  Vtest(e) += erfcm(alpha*r)/r;
	}
      }
    }
    //----done lattice vectors
  }

  //cout << "electron recip " << endl;
  
  //---------electron reciprocal part
  
  
  doublevar rdotg;

  for(int gpt=0; gpt < ngpoints; gpt++) {
    rdotg = gpoint(gpt, 0)*tpos(0)
            +gpoint(gpt, 1)*tpos(1)
            +gpoint(gpt, 2)*tpos(2);
    doublevar test_sin = sin(rdotg);
    doublevar test_cos = cos(rdotg);
    Vtest(totnelectrons) += 2.0*(-ion_cos(gpt)*test_cos - ion_sin(gpt)*test_sin + 0.5)*gweight(gpt);
    for(int e=0; e< totnelectrons; e++) {
      rdotg=gpoint(gpt, 0)*elecpos(e,0)
            +gpoint(gpt, 1)*elecpos(e,1)
            +gpoint(gpt, 2)*elecpos(e,2);
      Vtest(e) += 2.0*(test_cos*cos(rdotg) + test_sin*sin(rdotg))*gweight(gpt);
    }
  }
}


doublevar Periodic_system::ewaldElectron(Sample_point * sample) {
  //cout << sample << endl;
  sample->updateEEDist();
  //cout << "updateEIDIST" << endl;
  sample->updateEIDist();
  //cout << "done updatedist" << endl;

  Array1 <doublevar> eidist(5);
  int nions=ions.size();
  int nlatvec=1;
  Array1 <doublevar> r1(3), r2(3);

  //cout << "electron-ion " << endl;
  //-------------Electron-ion real part
  doublevar elecIon_real=0;
  Array1<doublevar> elecIon_real_separated(totnelectrons); 
  elecIon_real_separated = 0.0; 
  for(int e=0; e< totnelectrons; e++) {
    for(int ion=0; ion < nions; ion++) {

      sample->getEIDist(e,ion, eidist);
      for(int d=0; d< 3; d++) r1(d)=eidist(d+2);

      //----over  lattice vectors
      for(int kk=-nlatvec; kk <=nlatvec; kk++) {
        for(int jj=-nlatvec; jj <=nlatvec; jj++) {
          for(int ii=-nlatvec; ii <=nlatvec; ii++) {
            for(int d=0; d< 3; d++) {
              r2(d)=r1(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
            }
            doublevar r=sqrt(r2(0)*r2(0)+r2(1)*r2(1)+r2(2)*r2(2));
	    elecIon_real_separated(e) -= ions.charge(ion)*erfcm(alpha*r)/r;
            elecIon_real-=ions.charge(ion)*erfcm(alpha*r)/r;
          }
        }
      }
      //----done lattice vectors
    }
  }
  //cout << "electron-electron " << endl;
  //-------------Electron-electron real part

  doublevar elecElec_real=0;
  Array1 <doublevar> elecElec_real_separated(totnelectrons); 
  elecElec_real_separated = 0.0; 
  for(int e1=0; e1< totnelectrons; e1++) {
    for(int e2 =e1+1; e2 < totnelectrons; e2++) {
      sample->getEEDist(e1,e2, eidist);
      for(int d=0; d< 3; d++) r1(d)=eidist(d+2);

      //----over  lattice vectors
      for(int kk=-nlatvec; kk <=nlatvec; kk++) {
        for(int jj=-nlatvec; jj <=nlatvec; jj++) {
          for(int ii=-nlatvec; ii <=nlatvec; ii++) {
            for(int d=0; d< 3; d++) {
              r2(d)=r1(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
            }
            doublevar r=sqrt(r2(0)*r2(0)+r2(1)*r2(1)+r2(2)*r2(2));

            elecElec_real+=erfcm(alpha*r)/r;
	    elecElec_real_separated(e1) += erfcm(alpha*r)/r;
	    elecElec_real_separated(e2) += erfcm(alpha*r)/r;
          }
        }
      }
      //----done lattice vectors
    }
  }

  //cout << "electron recip " << endl;

  //---------electron reciprocal part


  doublevar rdotg;
  Array2 <doublevar> elecpos(totnelectrons, 3);
  sample->getAllElectronPos(elecpos);
  doublevar elecIon_recip=0, elecElec_recip=0;
  Array1 <doublevar> elecIon_recip_separated(totnelectrons), 
    elecElec_recip_separated(totnelectrons);
  elecIon_recip_separated = 0.0; 
  elecElec_recip_separated = 0.0; 
  for(int gpt=0; gpt < ngpoints; gpt++) {
    doublevar sum_sin=0, sum_cos=0;
    Array1<doublevar> sn(totnelectrons); 
    Array1<doublevar> cs(totnelectrons); 
    
    for(int e=0; e< totnelectrons; e++) {
      rdotg=gpoint(gpt, 0)*elecpos(e,0)
            +gpoint(gpt, 1)*elecpos(e,1)
            +gpoint(gpt, 2)*elecpos(e,2);
      sn(e) = sin(rdotg); 
      cs(e) = cos(rdotg);
      sum_sin += sn(e); 
      sum_cos += cs(e); 
    }
    elecIon_recip-=(ion_cos(gpt)*sum_cos + ion_sin(gpt)*sum_sin)*gweight(gpt);
    elecElec_recip+=(sum_cos*sum_cos + sum_sin*sum_sin)*gweight(gpt)/2;
    for (int e=0; e<totnelectrons; e++) {
      elecIon_recip_separated(e) -= (ion_cos(gpt)*cs(e) 
				     + ion_sin(gpt)*sn(e))*gweight(gpt);
      elecElec_recip_separated(e) += (sum_cos*cs(e) + sum_sin*sn(e)-0.5)*gweight(gpt);
      // this term is essentially the contribution from 
    }
  }

  elecIon_recip*=2;
  elecElec_recip*=2;
  //cout << "done " << endl;

  //cout << "elecElec_real " << elecElec_real << endl;
  //cout << "elecIon_real " << elecIon_real << endl;
  //cout << "elecElec_recip " << elecElec_recip << endl;
  //cout << "elecIon_recip " << elecIon_recip << endl;
  ewalde_separated.Resize(totnelectrons); 
  doublevar cc; 
  for (int e=0; e<totnelectrons; e++) {
    ewalde_separated(e) = elecElec_real_separated(e) + elecIon_real_separated(e)
         + 2.0*(elecIon_recip_separated(e) + elecElec_recip_separated(e)); 
    //
    //ewalde_separated(e) = elecElec_real_separated(e) 
    //  + 2.0*(elecElec_recip_separated(e)); 
        //noted that the elecElec_recip_separated term already removed the self interaction term; 
  }
  return elecElec_real + elecIon_real + elecElec_recip+elecIon_recip;
}

//------------------------------------------------------------------------


