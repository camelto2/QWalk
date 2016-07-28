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


#include "Maximize_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include <iomanip>
#include "Program_options.h"
#include "System.h"
#include "Split_sample.h"
#include "Properties.h"
#include "Generate_sample.h"

//----------------------------------------------------------------------

void Maximize_method::read(vector <string> words,
            unsigned int & pos,
            Program_options & options) { 
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);

  debug_write(cout, "wfdata allocate\n");
  wfdata=NULL;
  if(options.twftext.size() < 1) error("Need TRIALFUNC section for OPTIMIZE");
  allocate(options.twftext[0], sys, wfdata);

  if(!readvalue(words,pos=0,nconfig,"NCONFIG"))
    nconfig=100;
  

#ifdef USE_MPI
  int extra_configs = nconfig % mpi_info.nprocs;
  if(mpi_info.node < extra_configs) {
    nconfigs_per_node = (nconfig/mpi_info.nprocs) + 1;
  }
  else {
    nconfigs_per_node = (nconfig/mpi_info.nprocs);
  }
#else
  nconfigs_per_node = nconfig;
#endif
}

//----------------------------------------------------------------------

void Maximize_method::run(Program_options & options, ostream & output) { 
  Wavefunction * wf=NULL;
  Sample_point * sample=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  
  Wf_return wfret(1,2);

  Array1 <Config_save_point> config_pos(nconfigs_per_node);
  Array1 <Maximize_config> maximize_config(nconfigs_per_node);
  Array1 <doublevar> epos(3);

  Primary guidewf;
  generate_sample(sample,wf,wfdata,&guidewf,nconfig,config_pos);
  int nelectrons=sample->electronSize();
  Properties_gather mygather;
  
  int count = 1;
  Array2 <doublevar> tempconfig(nelectrons,3);
  Array2 <doublevar> temphessian(3*nelectrons,3*nelectrons);
  
  for(int i=0; i < nconfigs_per_node; i++) { 
    stringstream tableout;
    maximize(sample,wf,config_pos(i),temphessian);
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,epos);
      for(int d=0; d< 3; d++) {
        tempconfig(e,d) = epos(d);
        //tableout << epos(d) << " ";
      }
    }

    Properties_point pt;
    mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                            sample, &guidewf);
    
    wf->getVal(wfdata,0,wfret);
    
    //tableout << pt.energy(0) << " ";
    //tableout << wfret.amp(0,0) << endl;
    
    maximize_config(i).nelectrons = nelectrons;
    maximize_config(i).logpsi = wfret.amp(0,0);
    maximize_config(i).energy = pt.energy(0);
    maximize_config(i).config = tempconfig;
    maximize_config(i).hessian= temphessian;
     
    config_pos(i).write(cout);
  }
  
  string outfilename=options.runid+".json";
  write_configurations_maximize_json(outfilename, maximize_config);
  
  delete wf;
  delete sample;
}
//----------------------------------------------------------------------

#include "macopt.h"
class Maximize_fn:public Macopt {
public:
  Maximize_fn(int _n,int _verbose, double _tol,
              int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) { 
              } 
  ~Maximize_fn(){ };


  //-----------------------------------------
  void update_positions(double * _p) { 
    int count=1; 
    Array1 <doublevar> epos(3);
    int nelec=sample->electronSize();
    for(int e=0; e< nelec; e++) { 
      for(int d=0; d< 3; d++) { 
        epos(d)=_p[count++];
      }
      sample->setElectronPos(e,epos);
    }
    
  }
  //-----------------------------------------

  double func(double * _p) { 
    update_positions(_p);
    wf->updateVal(wfdata,sample);
    Wf_return wfret(1,2);
    wf->getVal(wfdata,0,wfret);
    return -wfret.amp(0,0);
  }
  //-----------------------------------------
  
  double dfunc(double * _p, double * _g) { 
    update_positions(_p);
    wf->updateLap(wfdata,sample);
    int count=0;
    int nelec=sample->electronSize();
    Wf_return lap(1,5);
    wf->getVal(wfdata,0,lap);
    //doublevar cutoff=0.2;
    _g[count++]=-lap.amp(0,0);
    for(int e=0; e< nelec; e++) { 
      wf->getLap(wfdata,e,lap);
      for(int d=0; d< 3; d++) {
        doublevar grad=-lap.amp(0,d+1);
        //if(abs(grad) > cutoff) {
        //  cout << "big grad " << grad << endl;
        //  grad=sign(grad)*cutoff;
        //}
        _g[count++]=grad;
      }
    }
    return _g[0];
  }
  //-----------------------------------------
  //move along the gradient for time tstep and evaluate the function
  doublevar eval_tstep(double * x, doublevar tstep,Array1 <doublevar> & grad,
                       int n) {
    Array1 <doublevar> xnew(n+1);
    for(int i=1; i<=n; i++)
      xnew(i)=x[i]-tstep*grad(i);
    return func(xnew.v);
  }

  //-----------------------------------------
  // check if the magnitude of the gradient is within tol
  bool grad_under_tol(Array1 <doublevar> & grad, int n, doublevar tol) {
    doublevar total = 0;
    for(int i=1; i<=n; i++)
      total += grad(i)*grad(i);
    return total < tol*tol;
  }

  doublevar grad_abs(Array1 <doublevar> & grad, int n) {
    doublevar total = 0;
    for(int i=1; i<=n; i++)
      total += grad(i)*grad(i);
    return sqrt(total);
  }
  
  //-----------------------------------------
  void macoptII(double * x,int n) {
    Array1 <doublevar> grad(n+1);
    //Array1 <doublevar> xnew(n+1);
    int max_it = 100;
    int max_big_it = 10000;
    for(int big_it=0; big_it < max_big_it; big_it++) {
      //find the direction of the gradient at x
      dfunc(x,grad.v);
      doublevar gradlen = grad_abs(grad,n);
      doublevar unit;
      doublevar fbase=grad(0);
      doublevar tol = 1e-12;
      doublevar outer_tol = tol*10; // make outer loop tolerance looser than inner
      doublevar bracket_tstep=0.0,last_func=fbase;
      //bracket the minimum in this direction (in 1D units of tstep) for bisection method
      for(doublevar tstep=1e-8; tstep < 20.0; tstep*=2.0) {
        
        doublevar f=eval_tstep(x,tstep,grad,n);
        cout << "tstep " << tstep << " func " << f << " fbase " << fbase << endl;
        if(f > fbase or f > last_func) {
          bracket_tstep=tstep;
          break;
        }
        else last_func=f;
      }
      
      cout << "bracket_tstep " << bracket_tstep << endl;
      //bisection method works best using the golden ratio
      doublevar resphi=2.-(1.+sqrt(5.))/2.;
      doublevar a=0, b=resphi*bracket_tstep, c=bracket_tstep;
      doublevar af=fbase,
      bf=eval_tstep(x,b,grad,n),
      cf=eval_tstep(x,c,grad,n);
      //if(af-bf<tolbf and cf-bf<tolbf and bf-af<tolbf and bf-cf<tolbf){
      //  break; //exit big loop if there is small change in the direction of gradient 
      //}
      //(grad psi)/psi = grad (log psi), small change in grad direction ( var grad is grad(log psi))
      if(gradlen < outer_tol) {break;} //Exit big loop if (grad psi)/psi is small 
      cout << "first step  a,b,c " << a << " " << b << "  " << c
      << " funcs " << af << " " << bf << " " << cf << endl;
      //bisection method iteration, (a, b, c)
      for(int it=0; it < max_it; it++) {
        doublevar d,df;
        //choose point d on the bigger interval (a,b) or (b,c)
        if( (c-b) > (b-a))
          d=b+resphi*(c-b);
        else
          d=b-resphi*(b-a);
        df=eval_tstep(x,d,grad,n);
        //check if the function is smaller at d than at b to choose next bracket
        if(df < bf) {
          if( (c-b) > (b-a) ) {
            // How to check tolerance: log(1+x) ~ x. Want to make sure to check relative change in psi
            //(af-bf) = log psi_a - log psi_b = log(psi_a/psi_b) ~ (psi_a-psi_b)/psi_b < tol
            unit = gradlen*(b-a);
            if(af/unit-bf/unit<tol and bf/unit-af/unit<tol) { cout << it << " (af-bf)/(b-a)=" << af/unit-bf/unit << endl; break; }
            a=b;
            af=bf;
            b=d;
            bf=df;
          }
          else {
            unit = gradlen*(c-b);
            if(cf/unit-bf/unit<tol and bf/unit-cf/unit<tol) { cout << it << " (cf-bf)/(c-b)=" << cf/unit-bf/unit << endl; break; }
            c=b;
            cf=bf;
            b=d;
            bf=df;
          }
        }
        //if function is bigger at d than at b, make d the new bracket boundary
        else {
          if( (c-b) > (b-a) ) {
            unit = gradlen*(c-d);
            if(cf/unit-df/unit<tol and df/unit-cf/unit<tol) { cout << it << " (df-cf)/(d-c)=" << df/unit-cf/unit << endl; break; }
            c=d;
            cf=df;
          }
          else {
            unit = gradlen*(d-a);
            if(af/unit-df/unit<tol and df/unit-af/unit<tol) { cout << it << " (af-df)/(d-a)=" << af/unit-df/unit << endl; break; }
            a=d;
            af=df;
          }
        }
        cout << "step " << it << " a,b,c " << a << " " << b << "  " << c
        << " funcs " << af << " " << bf << " " << cf << endl;
        if(it==max_it-1) {
          cout << "Warning: inner loop did not reach tolerance" << endl;
        }
      }
      cout << "last step b-a,c-b " << b-a << " " << c-b
      << " func diffs " << af-bf << " " << cf-bf << " " << endl;
      //finished bisection search, minimum at b; compute x for tstep b
      doublevar best_tstep=b;
      for(int i=1; i<=n; i++)
        x[i]=x[i]-best_tstep*grad[i];

      if(big_it==max_big_it-1) {
        doublevar grad_squared = 0;
        for (int i=1;i<n;i++){
          grad_squared += grad(i)*grad(i);
        }
        cout << "Warning: outer loop did not reach tolerance. |grad| = " << sqrt(grad_squared) << endl;
      }
    }
  }
  //-----------------------------------------
  
  void iteration_print(double f, double gg, double tol,  int itn) {
    cout << "It " << endl;
  }
  
  //-----------------------------------------
  void calc_hessian(double * x, Array2 <doublevar> & hessian, int n, double h=1e-3) {
    Array1 <doublevar> grad_plus(n+1);
    Array1 <doublevar> grad_minus(n+1);
    Array1 <doublevar> xnew(n+1);
    if(h<=0) { h = 1e-3; } // default step size 
    
    for(int i=1; i<=n; i++) {
      xnew(i)=x[i];
    }
    
    for(int i=1; i<=n; i++) {
      // forward step
      xnew(i) = x[i] + h;
      dfunc(xnew.v,grad_plus.v);
      // backward step
      xnew(i) = x[i] - h;
      dfunc(xnew.v,grad_minus.v);
      // reset position
      xnew(i) = x[i];
      // compute derivative
      for(int j=1; j<=n; j++) {
        hessian(i-1,j-1) = (grad_plus(j) - grad_minus(j)) / (2*h);
      }   
    }
  }
  //-----------------------------------------

  Wavefunction * wf;
  Wavefunction_data * wfdata;
  Sample_point * sample;
  System * system;
};

//----------------------------------------------------------------------

void Maximize_method::maximize(Sample_point * sample,Wavefunction * wf,Config_save_point & pt,Array2 <doublevar> & hessian) { 
  int nelectrons=sample->electronSize();
  pt.restorePos(sample);
  Maximize_fn maximizer(nelectrons*3,1,1e-12,1000,1);
  maximizer.wf=wf;
  maximizer.wfdata=wfdata;
  maximizer.sample=sample;
  maximizer.system=sys;

  Array1 <doublevar> allpos(nelectrons*3+1);
  Array1 <doublevar> epos(3);
  int count=1;
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) { 
      allpos(count++)=epos(d);
    }
  }

  maximizer.maccheckgrad(allpos.v,nelectrons*3,0.001,nelectrons*3);
  maximizer.macoptII(allpos.v,nelectrons*3);  
  maximizer.calc_hessian(allpos.v,hessian,nelectrons*3);

  pt.savePos(sample);
}

//----------------------------------------------------------------------
  
int Maximize_method::showinfo(ostream & os) { 
  return 1;
}

//----------------------------------------------------------------------
//void Maximize_method::calc_hessian(double * x, Array2 <doublevar> & hessian, int n) {
//  Array1 <doublevar> grad_plus(n+1);
//  Array1 <doublevar> grad_minus(n+1);
//  Array1 <doublevar> xnew(n+1);
//  double h = 1e-3; // step size
//  
//  for(int i=1; i<=n; i++) {
//    xnew(i)=x[i];
//  }
//  
//  for(int i=1; i<=n; i++) {
//    // forward step
//    xnew(i) = x[i] + h;
//    dfunc(xnew,grad_plus.v);
//    // backward step
//    xnew(i) = x[i] - h;
//    dfunc(xnew,grad_minus.v);
//    // reset position
//    xnew(i) = x[i];
//    // compute derivative
//    for(j=1; j<=n; j++) {
//      hessian(i-1,j-1) = (grad_plus(j) - grad_minus(j)) / (2*h)
//    }   
//  }
//}

//----------------------------------------------------------------------
void write_configurations_maximize_json(string & filename, Array1 <Maximize_config> configs) { 
  int nconfigs=configs.GetDim(0);
  time_t starttime;
  time(&starttime);
  string tmpfilename=filename; //+".qw_tomove";
  string backfilename=filename+".backup";
  if(mpi_info.node==0) { rename(tmpfilename.c_str(),backfilename.c_str()); }

  #ifdef USE_MPI
  stringstream os;
  os.precision(15);
  for(int i=0; i< nconfigs; i++) {
    if(i!=0) { os << ", "; }
    configs(i).write(os,false);
  }

  string walkstr=os.str();
  int nthis_string=walkstr.size();
  if(mpi_info.node==0) {
    ofstream os(tmpfilename.c_str());
    os << "[";
    os << walkstr;
    MPI_Status status;
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Recv(nthis_string,i);
      char * buf=new char[nthis_string+1];
      MPI_Recv(buf,nthis_string,MPI_CHAR, i, 0, MPI_Comm_grp, & status);
      buf[nthis_string]='\0';
      os << ", " << buf;
      delete [] buf;
    }
    os << "]";
  }
  else {
    MPI_Send(nthis_string,0);
    //we know that MPI_Send obeys const-ness, but the interfaces are not clean 
    // and so...casting!
    MPI_Send((char *) walkstr.c_str(),nthis_string, MPI_CHAR, 0,0,MPI_Comm_grp);
  }
#else
  ofstream os(tmpfilename.c_str());
  os.precision(15);
  os << "[";
  for(int i=0; i< nconfigs; i++) {
    if(i!=0) { os << ", "; }
    configs(i).write(os,false);
  }
  os << "]";
  os.close();
#endif
  time_t endtime;
  time(&endtime);
  if(mpi_info.node==1) {
    debug_write(cout, "Write took ", difftime(endtime, starttime), " seconds\n");
  }
}


void write_configurations_maximize_yaml(string & filename, Array1 <Maximize_config> configs) { 
  int nconfigs=configs.GetDim(0);
  time_t starttime;
  time(&starttime);
  string tmpfilename=filename; //+".qw_tomove";
  string backfilename=filename+".backup";
  if(mpi_info.node==0) { rename(tmpfilename.c_str(),backfilename.c_str()); }

  #ifdef USE_MPI
  stringstream os;
  os.precision(15);
  for(int i=0; i< nconfigs; i++) {
     configs(i).write(os,false);
  }

  string walkstr=os.str();
  int nthis_string=walkstr.size();
  if(mpi_info.node==0) {
    ofstream os(tmpfilename.c_str());
    os << "results:" << endl;
    os << walkstr;
    MPI_Status status;
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Recv(nthis_string,i);
      char * buf=new char[nthis_string+1];
      MPI_Recv(buf,nthis_string,MPI_CHAR, i, 0, MPI_Comm_grp, & status);
      buf[nthis_string]='\0';
      os << buf;
      delete [] buf;
    }
  }
  else {
    MPI_Send(nthis_string,0);
    //we know that MPI_Send obeys const-ness, but the interfaces are not clean 
    // and so...casting!
    MPI_Send((char *) walkstr.c_str(),nthis_string, MPI_CHAR, 0,0,MPI_Comm_grp);
  }
#else
  ofstream os(tmpfilename.c_str());
  os.precision(15);
  os << "results:" << endl;
  for(int i=0; i< nconfigs; i++) {
    configs(i).write(os,false);
  }
  os.close();
#endif
  time_t endtime;
  time(&endtime);
  if(mpi_info.node==1) {
    debug_write(cout, "Write took ", difftime(endtime, starttime), " seconds\n");
  }
}

//----------------------------------------------------------------------
void write_configurations_maximize(string & filename, 
                              Array1 <Maximize_config> configs) { 
  int nconfigs=configs.GetDim(0);
  time_t starttime;
  time(&starttime);
  string tmpfilename=filename; //+".qw_tomove";
  string backfilename=filename+".backup";
  if(mpi_info.node==0) { rename(tmpfilename.c_str(),backfilename.c_str()); }

  #ifdef USE_MPI
  stringstream os;
  os.precision(15);
  for(int i=0; i< nconfigs; i++) {
     //os << " walker { \n";
     configs(i).write(os,false);
     //os << " } \n";
   }

  string walkstr=os.str();
  int nthis_string=walkstr.size();
  if(mpi_info.node==0) {
    ofstream os(tmpfilename.c_str());
    for(int e=0; e< configs(0).nelectrons; e++) {
      for(int d=0; d< 3; d++) {
        os << "e" << e << "_"<< d << " ";
      }
    }
    os << "elocal" << " ";
    os << "logpsi" << endl;
    
    os << walkstr;
    MPI_Status status;
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Recv(nthis_string,i);
      char * buf=new char[nthis_string+1];
      MPI_Recv(buf,nthis_string,MPI_CHAR, i, 0, MPI_Comm_grp, & status);
      buf[nthis_string]='\0';
      os << buf;
      delete [] buf;
    }
  }
  else {
    MPI_Send(nthis_string,0);
    //we know that MPI_Send obeys const-ness, but the interfaces are not clean 
    // and so...casting!
    MPI_Send((char *) walkstr.c_str(),nthis_string, MPI_CHAR, 0,0,MPI_Comm_grp);
  }
  #else
    ofstream os(tmpfilename.c_str());
    os.precision(15);
    for(int i=0; i< nconfigs; i++) {
      //os << " walker { \n";
      configs(i).write(os,false);
      //os << " } \n";
    }
    os.close();
  #endif
    time_t endtime;
    time(&endtime);
    if(mpi_info.node==1) {
      debug_write(cout, "Write took ", difftime(endtime, starttime), " seconds\n");
    }
}

//---------------------------------------------------------------------
//YAML format
//void Maximize_config::write(ostream & os) {
//  os << "- " << "psi: " << logpsi << endl;
//  os << "  " << "energy: " << energy << endl;
//  os << "  " << "config: " << endl;
//  for(int e=0; e<nelectrons; e++) {
//    os << "  - - " << config(e,0) << endl;
//    for(int d=1; d<3; d++) {
//      os << "    - " << config(e,d) << endl;
//    }
//  }
//  os << "  " << "hessian: " << endl;
//  for(int c=0; c<3*nelectrons; c++) {
//    os << "  - - " << hessian(c,0) << endl;
//    for(int d=1; d<3*nelectrons; d++) {
//      os << "    - " << hessian(c,d) << endl;
//    }
//  }
//}

// JSON format
void Maximize_config::write(ostream & os,bool write_hessian) {
  os << "{";
  os << "\"psi\": " << logpsi;
  os << "," << "\"energy\": " << energy;
  os << "," << "\"config\": " << "[";
  for(int e=0; e<nelectrons; e++) {
    if(e!=0) { os << ", "; }
    os << "[";
    for(int d=0; d<3; d++) {
      if(d!=0) { os << ", "; }
      os << config(e,d); 
    }
    os << "]";
  }
  os << "]";
  if(write_hessian){
    os << "," << "\"hessian\": " << "[";
    for(int c=0; c<3*nelectrons; c++) {
      if(c!=0) { os << ", "; }
      os << "[";
      for(int d=0; d<3*nelectrons; d++) {
        if(d!=0) { os << ", "; }
        os << hessian(c,d); 
      }
      os << "]";
    }
    os << "]";
  }
  os << "}";
}

void Maximize_config::read(istream & is) {}
void Maximize_config::mpiSend(int node) {}
void Maximize_config::mpiReceive(int node) {}