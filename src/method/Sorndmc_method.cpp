/*
 
Copyright (C) 2007 Lucas K. Wagner
addition of the PURE DMC: Michal Bajdich and Fernando A. Reboredo

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



#include "Sorndmc_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include "Program_options.h"
#include "average.h"
#include "Generate_sample.h"


void Sorndmc_method::read(vector <string> words,
                      unsigned int & pos,
                      Program_options & options)
{
  ndim=3;

  have_read_options=1;

  vector <string> offsetwords;

  //required options

  if(!readvalue(words, pos=0, timestep, "TIMESTEP"))
    error("Need TIMESTEP in METHOD section");

  if(!readvalue(words, pos=0, nblock, "NBLOCK"))
    nblock=100;
  
  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
    nconfig=max(2048/mpi_info.nprocs,1);

  if(!readvalue(words, pos=0, nstep, "NSTEP")) 
    nstep=1;
  else
    error("Number of steps is set to 1 in RNDMC");
  
  if(!readvalue(words, pos=0, readconfig, "READCONFIG"))
      //CM
    //readconfig=options.runid+".config";
    //error("Must give READCONFIG for DMC");
    error("Must provide a config file for SORNDMC");

  //optional options

  if(!readvalue(words, pos=0, eref, "EREF"))
    eref=0.0;

  if(haskeyword(words, pos=0, "CDMC")) do_cdmc=1;
  else do_cdmc=0;
  if(haskeyword(words, pos=0, "TMOVES")) tmoves=1; 
  else tmoves=0;

  if(!readvalue(words, pos=0, nhist, "CORR_HIST")) 
    nhist=-1;

  if(!readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    storeconfig=options.runid+".config";

  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="sorndmc";

  if(!readvalue(words, pos=0, start_feedback, "START_FEEDBACK"))
    start_feedback=1;

  if(readvalue(words, pos=0, feedback_interval, "FEEDBACK_INTERVAL")) {
    if(feedback_interval < 1) 
      error("FEEDBACK_INTERVAL must be greater than or equal to 1");
  }
  else feedback_interval=5;

  if(!readvalue(words, pos=0, feedback, "FEEDBACK"))
    feedback=1.0;

  if(!readvalue(words, pos=0, branch_start_cutoff, "BRANCH_START_CUTOFF")) 
    branch_start_cutoff=10;
  

  branch_stop_cutoff=branch_start_cutoff*1.5;
  
  
  vector <string> proptxt;
  if(readsection(words, pos=0, proptxt, "PROPERTIES")) 
    mygather.read(proptxt);

  vector<string> tmp_dens;
  pos=0;
  while(readsection(words, pos, tmp_dens, "DENSITY")) {
    dens_words.push_back(tmp_dens);
  }

  pos=0;
  while(readsection(words, pos, tmp_dens, "NONLOCAL_DENSITY")) {
    nldens_words.push_back(tmp_dens);
  }

  pos=0;
  while(readsection(words, pos, tmp_dens, "AVERAGE")) {
    avg_words.push_back(tmp_dens);
  }
  
  vector <string> dynamics_words;
  if(!readsection(words, pos=0, dynamics_words, "DYNAMICS") ) 
    dynamics_words.push_back("SPLIT");

  low_io=0;
  if(haskeyword(words, pos=0,"LOW_IO")) low_io=1;

  if (!readsection(words, pos=0,guiding_words, "GUIDING_WF"))
      error("Must include GUIDING_WF");

  if(!readvalue(words, pos=0, nblock, "MAX_NODAL_CROSS_AGE"))
      max_nodal_cross_age = 0;

  allocate(dynamics_words, dyngen);

}

//----------------------------------------------------------------------

int Sorndmc_method::generateVariables(Program_options & options) {

  if(!have_read_options) 
    error("need to call Sorndmc_method::read before generateVariables");
  if(have_allocated_variables) 
    error("already allocated variables in Sorndmc_method");

  have_allocated_variables=1;
  allocate(options.systemtext[0], mysys);
  mysys->generatePseudo(options.pseudotext, mypseudo);
  allocate(options.twftext[0], mysys, mywfdata); 
  allocate(guiding_words, mysys, mygwfdata);
  
  densplt.Resize(dens_words.size());
  for(int i=0; i< densplt.GetDim(0); i++) {
    allocate(dens_words[i], mysys, options.runid,densplt(i));
  }
  nldensplt.Resize(nldens_words.size());
  for(int i=0; i< nldensplt.GetDim(0); i++) {
    allocate(nldens_words[i], mysys, options.runid,nldensplt(i));
  }
  
  return 1;
}

//----------------------------------------------------------------------



int Sorndmc_method::allocateIntermediateVariables(System * sys,
                                              Wavefunction_data * wfdata,
					      Wavefunction_data * gwfdata) {
  if(wf) delete wf;
  wf=NULL;
  if(gwf) delete gwf;
  gwf=NULL;
  if(sample) delete sample;
  sample=NULL;
  if(gsample) delete gsample;
  gsample=NULL;
  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  wfdata->generateWavefunction(wf);
  gwfdata->generateWavefunction(gwf);
  sys->generateSample(sample);
  sys->generateSample(gsample);
  sample->attachObserver(wf);
  gsample->attachObserver(gwf);
  nwf=wf->nfunc();

  if(wf->nfunc() >1) 
    error("DMC doesn't support more than one guiding wave function");

  guidingwf=new Primary; //Keep primary, the node release GF is sqrt(1+alpha^2)

  pts.Resize(nconfig);
  for(int i=0; i < nconfig; i++) {
    pts(i).age.Resize(nelectrons);
    pts(i).age=0;
  }
  
  average_var.Resize(avg_words.size());
  average_var=NULL;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    allocate(avg_words[i], sys, wfdata, average_var(i));
  }
  

  return 1;
}

//----------------------------------------------------------------------

int Sorndmc_method::showinfo(ostream & os)
{

  if(have_allocated_variables) {
    mysys->showinfo(os);
    mypseudo->showinfo(os);
    mywfdata->showinfo(os);
    mygwfdata->showinfo(os);
  }
  os << "###########################################################\n";
  os << "Spin-orbital Release Node Diffusion Monte Carlo:\n";
  os << endl;
  os << "This method departs from the traditional release node method" << endl;
  os << "by replacing the guiding wave function from normal RN with " << endl;
  os << "a spinor guiding function (so it is bosonic)" << endl;
  os << "The configuration space is expanded by adding the spin degree " << endl;
  os << "of freedom and we control the exploration of that by the spin time step" << endl;
  os << "Number of processors " <<           mpi_info.nprocs << endl;
  os << "Blocks: " <<                        nblock    << endl;
  os << "Steps per block: " <<               nstep     << endl;
  os << "Timestep: " <<                      timestep  << endl;
  if(tmoves) 
    os << "T-moves turned on" << endl;
  string indent="  ";

  dyngen->showinfo(indent, os);

  os << "###########################################################" << endl;
  return 1;
}

//----------------------------------------------------------------------

void Sorndmc_method::find_cutoffs() {
  doublevar eaverage=0;
  doublevar rneaverage = 0;
  
  doublevar totweight=0;
  doublevar rntotweight = 0;
  for(int i=0; i< nconfig; i++) {
      //eaverage+=(prop.trace(step,i).energy(w)+offset(w))
      //  *guidingwf->getOperatorWeight(prop.trace(step,i).wf_val,w);
    eaverage+=pts(i).prop.energy(0)*pts(i).weight;
    rneaverage+=exp(pts(i).prop.wf_val.amp(0,0)-pts(i).gprop.wf_val.amp(0,0))*pts(i).prop.energy(0)*pts(i).weight;
    totweight+=pts(i).weight;
    rntotweight+=exp(pts(i).prop.wf_val.amp(0,0)-pts(i).gprop.wf_val.amp(0,0))*pts(i).weight;
    //cout << "en " << pts(i).prop.energy(0) << endl;
  }
  //cout << mpi_info.node << "par sum " << nconfig << endl; 
  totweight=parallel_sum(totweight);
  rntotweight=parallel_sum(rntotweight);
  //cout << mpi_info.node << "eaverage " << endl;
  eaverage=parallel_sum(eaverage)/totweight;
  rneaverage=parallel_sum(rneaverage)/rntotweight;

  eref=eaverage;
  single_write(cout, " setting eref= ", eref, "\n");
  single_write(cout, " rneref= ", rneaverage, "\n");
  int totconf=parallel_sum(nconfig);
  doublevar eaverage2=0;
  for(int i=0; i < nconfig; i++){
    doublevar effenergy=pts[i].prop.energy(0);
    eaverage2+=(effenergy-eaverage)
               *(effenergy-eaverage);
  }
  
  eaverage2=parallel_sum(eaverage2)/totconf;

  //The variance of the starting distribution.
  doublevar sigmac=sqrt(eaverage2);

  branchcut_start=branch_start_cutoff*sigmac; //start of cutoff region for branching
  branchcut_stop=branch_stop_cutoff*sigmac; //end of cutoff region  

  single_write(cout, " start branch cut at ", branchcut_start, "\n");
  single_write(cout, " stop branch cut at ", branchcut_stop, "\n");


}

//----------------------------------------------------------------------

void Sorndmc_method::run(Program_options & options, ostream & output) {
  if(!have_allocated_variables) 
    error("Must generate variables to use Sorndmc_method::run");
  string logfile=options.runid+".fermionic.log";
  string logfile2=options.runid+".bosonic.log";
  string logfile3=options.runid+".absolute.log";
  
  if(mpi_info.node==0 ) {
    ofstream logout(logfile.c_str(), ios::app);
    logout << "#-------------------------------------------------\n";
    logout << "#SORNDMC fermionic: timestep " << timestep 
           << endl;
    logout << "#-------------------------------------------------\n\n\n";
    logout.close();
    ofstream logout2(logfile2.c_str(), ios::app);
    logout2 << "#-------------------------------------------------\n";
    logout2 << "#SORNDMC bosonic: timestep " << timestep 
           << endl;
    logout2 << "#-------------------------------------------------\n\n\n";
    logout2.close();
    ofstream logout3(logfile3.c_str(), ios::app);
    logout3 << "#-------------------------------------------------\n";
    logout3 << "#SORNDMC absolute: timestep " << timestep 
	   << endl;
    logout3 << "#-------------------------------------------------\n\n\n";
    logout3.close();
  }

  myprop_f.setLog(logfile, log_label);
  myprop_b.setLog(logfile2, log_label);
  myprop_a.setLog(logfile3, log_label);
  runWithVariables(myprop_f, mysys, mywfdata, mygwfdata, mypseudo,output);
}

//----------------------------------------------------------------------

/*!

 */
void Sorndmc_method::runWithVariables(Properties_manager & prop,
                                  System * sys, 
                                  Wavefunction_data * wfdata,
				  Wavefunction_data * gwfdata,
                                  Pseudopotential * pseudo,
                                  ostream & output)
{

  allocateIntermediateVariables(sys, wfdata, gwfdata);
  if(!wfdata->supports(laplacian_update) || !gwfdata->supports(laplacian_update))
    error("DMC doesn't support all-electron moves..please"
          " change your wave function to use the new Jastrow");

  cout.precision(15);
  output.precision(10);

  string maxEntname = log_label+"_maxEnt.dat";
  ofstream maxEnt(maxEntname.c_str());

  string stats_name = log_label+"_accept_efficiency.dat";
  ofstream stats(stats_name.c_str());

  
  prop.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata); 
  myprop_b.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata); 
  myprop_a.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, 
	       wfdata); 

  restorecheckpoint(readconfig, sys, wfdata, gwfdata, pseudo);
  prop.initializeLog(average_var);
  myprop_b.initializeLog(average_var);
  myprop_a.initializeLog(average_var);
  
  // Get signs from the trial wave function
  Array1<doublevar> walk_en_psiG(nconfig);
  Array1<doublevar> ini_en_psiT(nconfig);
  Array1<doublevar> ini_ratio(nconfig);
  for(int walker=0; walker < nconfig; walker++) {
    pts(walker).sign=pts(walker).prop.wf_val.sign(0);
    walk_en_psiG(walker) = pts(walker).gprop.energy(0);
    ini_en_psiT(walker) = pts(walker).prop.energy(0);
    ini_ratio(walker) = exp(pts(walker).prop.wf_val.amp(0,0) - pts(walker).gprop.wf_val.amp(0,0));
  }

  nhist=1;
  
  doublevar teff=timestep;
  doublevar time = 0.0;
  int total_config = parallel_sum(nconfig);

  for(int block=0; block < nblock; block++) {

    int totkilled=0;  
    int totbranch=0;
    int totpoints=0;
    for(int step=0; step < nstep; ) {
      int npsteps=min(feedback_interval, nstep-step);

      Dynamics_info dinfo;
      doublevar acsum=0;
      doublevar deltar2=0;
      Array1 <doublevar> epos(3);
      
      doublevar avg_acceptance=0;
      
      for(int walker=0; walker < nconfig; walker++) {
        //pts(walker).config_pos.restorePos(sample);
        //wf->updateLap(wfdata, sample);
	pts(walker).config_pos.restorePos(gsample);
	gwf->updateLap(gwfdata,gsample);
	//------Do several steps without branching
        for(int p=0; p < npsteps; p++) {
          pseudo->randomize();
          
          for(int e=0; e< nelectrons; e++) {
            int acc;
            acc=dyngen->sample(e, gsample, gwf, gwfdata, guidingwf,
                               dinfo, timestep); //MOVE USING THE GWF
            
            if(dinfo.accepted) 
              deltar2+=dinfo.diffusion_rate/(nconfig*nelectrons*npsteps);
            
            
            if(dinfo.accepted) { 
              pts(walker).age(e)=0;
            }
            else { 
              pts(walker).age(e)++;
            }
            avg_acceptance+=dinfo.acceptance/(nconfig*nelectrons*npsteps);
            
            if(acc>0) acsum++;
          }
          totpoints++;
          Properties_point pt;
	  Properties_point gpt;
          vector <Tmove> tmov;
          doublevar subtract_out_enwt=0;
          if(tmoves) {  //------------------T-moves
	      //Not supported for this
	    error("T-moves not supported currently for this");
	      /*
            pt.setSize(nwf);
            wf->getVal(wfdata,0,pt.wf_val);
            sys->calcKinetic(wfdata,sample,wf,pt.kinetic);
            pt.potential=sys->calcLoc(sample);
            pt.weight=1.0; //this gets set later anyway
            pt.count=1;
            pseudo->calcNonlocTmove(wfdata,sys,sample,wf,pt.nonlocal,tmov);
            //cout << "choosing among " <<  tmov.size() << " tmoves " << endl;
            //Now we do the t-move
            doublevar sum=1; 
            for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) { 
              assert(mov->vxx < 0);
              sum-=timestep*mov->vxx;  
            }
            pt.nonlocal(0)-=(sum-1)/timestep;
            subtract_out_enwt=-(sum-1)/timestep;
            //cout << "sum " << sum <<  " nonlocal " << pt.nonlocal(0) << " ratio " << sum/pt.nonlocal(0) << endl;
            assert(sum >= 0);
            doublevar rand=rng.ulec()*sum;
            sum=1; //reset to choose the move
            if(rand > sum) { 
              for(vector<Tmove>::iterator mov=tmov.begin(); mov!=tmov.end(); mov++) { 
                sum-=timestep*mov->vxx;
                if(rand < sum) { 
                  Array1 <doublevar> epos(3);
                  sample->getElectronPos(mov->e, epos);
                  //cout << "moving electron " << mov->e << " from " << epos(0) << " " << epos(1)
                  //  << " " << epos(2) << " to " << mov->pos(0) << " " << mov->pos(1) 
                  //  << " " << mov->pos(2) << endl;
                  sample->setElectronPos(mov->e,mov->pos);
                  break;
                }
              }
            }
            //wf->updateLap(wfdata, sample);
	    */
          } ///---------------------------------done with the T-moves
          else { 
	    Config_save_point tmp_config;
	    tmp_config.savePos(gsample); //transfer information from gsample to sample
	    tmp_config.restorePos(sample);
            mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                                sample, guidingwf); //Gathering data (local energy) from PsiT
            mygather.gatherData(gpt, pseudo, sys, gwfdata, gwf, 
                                gsample, guidingwf); //Gathering data (local energy) from PsiG
	    //CM
	    //gatherData sets myprop.weight to be |Psi_T|/|Psi_T| because of the way they structured it. 
	    //myprop.weight = 1
          }
          Dmc_history new_hist;
          new_hist.main_en=pts(walker).prop.energy(0); //previous local energy with PsiT
          pts(walker).past_energies.push_front(new_hist);
          deque<Dmc_history> & past(pts(walker).past_energies);
          if(past.size() > nhist) 
            past.erase(past.begin()+nhist, past.end());

	  new_hist.main_en=pts(walker).gprop.energy(0); //previous local energy with PsiG
	  pts(walker).gpast_energies.push_front(new_hist);
	  deque<Dmc_history> & gpast(pts(walker).gpast_energies);
	  if(gpast.size() > nhist)
	      gpast.erase(gpast.begin()+nhist,gpast.end());
          

          pts(walker).prop=pt;   //Prop weight 1
	  pts(walker).gprop=gpt; //Prop weight 1
	  pts(walker).prop.weight(0)*=exp(pts(walker).prop.wf_val.amp(0,0)-pts(walker).gprop.wf_val.amp(0,0)); // |psiT/psiG|
	  pts(walker).prop.weight(0)*=pts(walker).sign*pts(walker).prop.wf_val.sign(0); // s_i * sign(psiT)


	  //CM We are trying to do normal DMC, but now guiding function is used for weights. Need to update with PSI_G,
	  pts(walker).weight*=getWeight(pts(walker),teff,etrial); //cummulative walker weight for branch process

	  walk_en_psiG(walker) += pts(walker).gprop.energy(0);


          if(pts(walker).ignore_walker) {
            pts(walker).ignore_walker=0;
            pts(walker).weight=1;
            pts(walker).prop.count=0;
          }


	  //pts(walker).weight*=pts(walker).prop.weight(0);


          //pts(walker).prop.weight=pts(walker).weight;
	  //multiply the property weight by the accumulated walker weight.
          pts(walker).prop.weight(0)*=pts(walker).weight;
	  pts(walker).gprop.weight(0)*=pts(walker).weight;
          //This is somewhat inaccurate..will need to change it later
          //For the moment, the autocorrelation will be slightly
          //underestimated
          pts(walker).prop.parent=walker;
          pts(walker).prop.nchildren=1;
          pts(walker).prop.children(0)=walker;
          pts(walker).prop.avgrets.Resize(1,average_var.GetDim(0));

	  pts(walker).gprop.parent=walker;
	  pts(walker).gprop.nchildren=1;
	  pts(walker).gprop.children(0)=walker;

          for(int i=0; i< average_var.GetDim(0); i++) { 
            average_var(i)->evaluate(wfdata, wf, sys, sample, pts(walker).prop.avgrets(0,i));
          }

	  Properties_point tmp_prop_absolute;
	  tmp_prop_absolute = pts(walker).prop;
	  tmp_prop_absolute.weight(0) = fabs(pts(walker).prop.weight(0));

	  //Insert current properties for averaging
          prop.insertPoint(step+p, walker, pts(walker).prop);
          myprop_b.insertPoint(step+p, walker, pts(walker).gprop);
	  myprop_a.insertPoint(step+p, walker, tmp_prop_absolute);
          for(int i=0; i< densplt.GetDim(0); i++)
            densplt(i)->accumulate(sample,pts(walker).prop.weight(0));
          for(int i=0; i< nldensplt.GetDim(0); i++)
            nldensplt(i)->accumulate(sample,pts(walker).prop.weight(0),
                                     wfdata,wf);

        }

        pts(walker).config_pos.savePos(sample);

	if (max_nodal_cross_age != 0) {
	  if (pts(walker).sign*pts(walker).prop.wf_val.sign(0) > 0)
	      pts(walker).nodal_cross_age = 0;
	  else
	      pts(walker).nodal_cross_age += 1;

	  if (pts(walker).nodal_cross_age >= max_nodal_cross_age) 
	     pts(walker).sign *= -1;
	}

      }
      //---Finished moving all walkers

      doublevar accept_ratio=acsum/(nconfig*nelectrons*npsteps);
      teff=timestep*accept_ratio; //deltar2/rf_diffusion; 

      updateEtrial(feedback);

      step+=npsteps;

      int nkilled=0;
      //nkilled=calcBranch();  //NO BRANCHING
      
      totkilled+=nkilled;
      totbranch+=nkilled;

      time += timestep;
      doublevar h0 = 0.0;
      doublevar h1 = 0.0;
      for (int walker = 0; walker < nconfig; walker++) {
	  doublevar ratio = exp(pts(walker).prop.wf_val.amp(0,0)-pts(walker).gprop.wf_val.amp(0,0));
	  h0 += ini_ratio(walker)*ratio*exp(-timestep*walk_en_psiG(walker));
	  h1 += ini_ratio(walker)*ratio*(0.5*(ini_en_psiT(walker)+pts(walker).prop.energy(0)))*exp(-timestep*walk_en_psiG(walker));
      }
      doublevar total_h0 = parallel_sum(h0)/total_config;
      doublevar total_h1 = parallel_sum(h1)/total_config;
      single_write(maxEnt,time,"  ",total_h0);
      single_write(maxEnt,"  ",total_h1,"\n");

    }

    ///----Finished block
    
    if(!low_io || block==nblock-1) {
      savecheckpoint(storeconfig,sample);
      for(int i=0; i< densplt.GetDim(0); i++)
        densplt(i)->write();
      for(int i=0; i< nldensplt.GetDim(0); i++)
        nldensplt(i)->write(log_label);
    }

    prop.endBlock();
    myprop_b.endBlock();
    myprop_a.endBlock();

    Properties_block lastblock;
    Properties_block lastblock2;
    Properties_block lastblock3;
    prop.getLastBlock(lastblock);
    myprop_b.getLastBlock(lastblock2);
    myprop_a.getLastBlock(lastblock3);
    doublevar efficiency;
    doublevar weight = lastblock.avg(Properties_types::weight,0);
    doublevar weight_abs = lastblock3.avg(Properties_types::weight,0);
    efficiency = weight/weight_abs;

    totbranch=parallel_sum(totbranch);
    totkilled=parallel_sum(totkilled);
    totpoints=parallel_sum(totpoints);

    Properties_final_average finavg;
    Properties_final_average finavg2;
    Properties_final_average finavg3;
    prop.getFinal(finavg);
    myprop_b.getFinal(finavg2);
    myprop_a.getFinal(finavg3);
    eref=finavg.avg(Properties_types::total_energy,0);
    updateEtrial(feedback);
    
    doublevar maxage=0;
    doublevar avgage=0;
    for(int w=0;w < nconfig; w++) {
      for(int e=0; e< nelectrons; e++) { 
        if(maxage<pts(w).age(e)) maxage=pts(w).age(e);
        avgage+=pts(w).age(e);
      }
    }
    avgage/=(nconfig*nelectrons);

    if(output) {
      //cout << "Block " << block 
      //       << " nconfig " << totconfig
      //       << " etrial " << etrial << endl;
      if(global_options::rappture ) { 
	    ofstream rapout("rappture.log");
        rapout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(block+1)/doublevar(nblock))
               << "  Diffusion Monte Carlo" << endl;
        cout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(block+1)/doublevar(nblock))
             << "  Diffusion Monte Carlo" << endl;
        rapout.close();
      }
      output << "***" << endl;
      output << "Block " << block 
             << " etrial " << etrial << endl;
      output << "maximum age " << maxage 
	     << " average age " << avgage << endl;
      dyngen->showStats(output);

      prop.printBlockSummary(output);

      output << "Branched "
	     << totbranch << " times.  So a branch every " 
	     << doublevar(totpoints)/doublevar(totbranch)
	     << " steps " << endl;

      output << "RN efficiency: " <<  efficiency << endl;
    }

    doublevar acc = dyngen->get_accept();
    single_write(stats,time,"  ",acc);
    single_write(stats,"  ",efficiency,"\n");

    dyngen->resetStats();

  }
  
  if(output) {
    output << "\n ----------Finished DMC------------\n\n";
    prop.printSummary(output,average_var);  
  }
  wfdata->clearObserver();
  gwfdata->clearObserver();
  deallocateIntermediateVariables();

  maxEnt.close();
  stats.close();
}


//----------------------------------------------------------------------


void Sorndmc_method::savecheckpoint(string & filename,                     
                                 Sample_point * config) {
  if(filename=="") return;
  
  write_configurations(filename, pts);
  return;
  /*
  ofstream checkfile(filename.c_str());
  if(!checkfile) error("Couldn't open", filename );
  checkfile.precision(15);
  
  long int is1, is2;
  rng.getseed(is1, is2);
  checkfile << "RANDNUM " << is1 << "  " << is2 << endl;

  checkfile.precision(15);
  for(int i=0; i< nconfig; i++) { 
    Dmc_point & mypt(pts(i));
    checkfile << "SAMPLE_POINT { \n";
    mypt.config_pos.restorePos(config);
    write_config(checkfile, config);
    checkfile << "   DMC { \n";
    checkfile << "DMCWEIGHT " << mypt.weight << endl;
    checkfile << "VALEN " << nwf << endl; 
    for(int w=0; w< nwf; w++) {
      checkfile << mypt.prop.wf_val.phase(w,0) << "  "
		<< mypt.prop.wf_val.amp(w,0) << "  "
		<< mypt.prop.energy(w)
		<< endl;
    }

    checkfile << "   } \n";
    checkfile << "}\n\n";
  }

  checkfile.close();
   */
}


//----------------------------------------------------------------------



void Sorndmc_method::restorecheckpoint(string & filename, System * sys,
                                    Wavefunction_data * wfdata,
				    Wavefunction_data * gwfdata,
                                    Pseudopotential * pseudo) {

  ifstream is(filename.c_str());
  if(is) { 
    is.close();
    read_configurations(filename, pts);
  }
  else { 
    Array1 <Config_save_point>  configs;
    generate_sample(sample,wf,wfdata,guidingwf,nconfig,configs);
    pts.Resize(nconfig);
    for(int i=0; i< nconfig; i++) 
      pts(i).config_pos=configs(i);
  }
  int ncread=pts.GetDim(0);
  
  //cout << "ncread " << ncread << "  nwread " << nwread << endl;
  if(nconfig < ncread) { 
    Array1 <Rndmc_point> tmp_pts(nconfig);
    for(int i=0; i< nconfig; i++) tmp_pts(i)=pts(i);
    pts=tmp_pts;
  }
  else if(nconfig > ncread) { 
    error("Not enough configurations in ", filename);
  }

  for(int walker=0; walker < nconfig; walker++) {
    pts(walker).config_pos.restorePos(sample);
    pts(walker).config_pos.restorePos(gsample);
    mygather.gatherData(pts(walker).prop, pseudo, sys,
                        wfdata, wf, sample,
                        guidingwf);
    mygather.gatherData(pts(walker).gprop, pseudo, sys,
                        gwfdata, gwf, gsample,
                        guidingwf);
    pts(walker).age.Resize(sys->nelectrons(0)+sys->nelectrons(1));
    pts(walker).age=0;
  }
  find_cutoffs();

  updateEtrial(start_feedback);
    
}


//----------------------------------------------------------------------

void Sorndmc_method::cdmcReWeight(Array2 <doublevar> & energy_temp, 
                              Array1 < Wf_return > & value_temp
                              ) {
  //The following is for C-DMC(from Jeff)  It shouldn't
  //do anything unless the atomic coordinates have moved
  //It's commented out for the moment from the rewrite.
  //It's viewed as experimental code, so watch out!

  //Get the energies for the old ionic configurations and 
  //this configuration
  /*
  Array1 <doublevar> effenergy_temp(nconfig, 0.0);
  Array1 <doublevar> effoldenergy(nconfig,0.0);

  for(int i=0; i< nconfig; i++) {
    for(int w=0; w< nwf; w++) {
      effenergy_temp(i)+=(energy_temp(i,w)+offset(w))
        *guidingwf->getOperatorWeight(value_temp(i),w);

      doublevar olden=trace[0][i].energy(w);
      effoldenergy(i)+=(olden+offset(w))
        *guidingwf->getOperatorWeight(trace[0][i].wf_val, w);
    }
  }

  

  doublevar average_temp, variance_temp;

  average(0, nconfig, effenergy_temp,
          dmcweight, average_temp, variance_temp, 1);

  doublevar sigma_temp;
  if(variance_temp >0) {
    sigma_temp=sqrt(variance_temp);
  }
  else {
    sigma_temp=0;
    error("negative variance when reading in weights");
  }

  debug_write(cout, "average_temp ", average_temp, "\n");
  debug_write(cout, "sigma_temp ", sigma_temp, "\n");

  //Reweight and check whether we had a sign flip
  //(overlap will be either positive or negative

  if(do_cdmc) {
    doublevar norm1=0, norm2=0;
    for(int i=0; i< nconfig; i++) {
      //norm1+=exp(prop.trace(0,i).wf_val(0,1)*2.0);
      norm1+=exp(trace[0][i].wf_val.amp(0,0)*2.0);
      
      norm2+=exp(2.0*value_temp(i).amp(0,0));
    }
    
    int totconfig=parallel_sum(nconfig);
    norm1=parallel_sum(norm1)/totconfig;
    norm2=parallel_sum(norm2)/totconfig;
    
    //Removed the normalization(norm2/norm1), because it causes problems
    //for bigger systems
    
    doublevar overlap=0;
    for(int i=0; i < nconfig; i++) {
      doublevar ratio;
      ratio=guidingwf->getTrialRatio(trace[0][i].wf_val, value_temp(i));
      overlap+=ratio/nconfig; //check if there's an overall sign change
      
      //cout << "walker " << i << " new val " << trace[0][i].wf_val(0,1) 
      //     << " old one " << value_temp(i)(0,1) << endl;

      dmcweight(i)*= ratio*ratio//  *norm2/norm1
        *exp(-(effoldenergy(i)-effenergy_temp(i))*timestep/2);
    }
    
    single_write(cout, "overlap between old and new ", overlap);
    doublevar average_old=0.0, diff_sigma=0.0;
     
    for(int w=0; w< nconfig; w++)  {
       average_old+=effoldenergy(w)/nconfig;
       diff_sigma+=(effoldenergy(w)-effenergy_temp(w))*(effoldenergy(w)-effenergy_temp(w))/nconfig;
    }
    //cout << "difference " << average_temp-average_old <<"  +/-  " << diff_sigma/nconfig << endl;

    //We ignore walkers that either a)cross a node, or 
    //b) are outside of 10 sigmas(for the first move)
    //ofstream diffout("diff_en.dat");
    int ncross=0, nsigma=0;
    for(int i=0; i< nconfig; i++) {
      doublevar ratio;
      ratio=guidingwf->getTrialRatio(trace[0][i].wf_val, value_temp(i));
      //cout << "ratio " << i << "   " << ratio << endl;
      if(ratio*overlap < 0) {
        debug_write(cout, "crossed node ", i , "\n");
        ncross++;
        trace[0][i].count=0;
        dmcweight(i)=1;
        ignore_walker(i)=1;
      }
      else if(fabs(effenergy_temp(i)-average_temp) > 10*sigma_temp) {
        debug_write(cout, "walker out of 10 sigmas ", i, "\n");
	
        nsigma++;
        trace[0][i].count=0;
        dmcweight(i)=1;
        ignore_walker(i)=1;
      }
      // diffout << i << "   " << effenergy_temp(i)-effoldenergy(i) << endl;
    }
    //diffout.close();

    single_write(cout, "nsigma ", parallel_sum(nsigma));
    single_write(cout, "  ncross  ", parallel_sum(ncross), "\n");

  }

  */
  

}


//----------------------------------------------------------------------

void Sorndmc_method::updateEtrial(doublevar feed) {
  
  doublevar totweight=0;
  for(int walker=0; walker < nconfig; walker++)
    totweight+=pts(walker).weight;

  etrial=eref-feed*log(totweight/double(nconfig));

#ifdef USE_MPI
  doublevar mpitot=0;
  MPI_Allreduce(&totweight,&mpitot, 1,
                MPI_DOUBLE, MPI_SUM, MPI_Comm_grp);
  etrial=eref-feed*log(mpitot/double(mpi_info.nprocs*nconfig));
#endif
  
  //cout << "etrial " << etrial << " total weight " << totweight << endl;

}

//----------------------------------------------------------------------

doublevar Sorndmc_method::getWeight(Rndmc_point & pt,
                                doublevar teff, doublevar etr) {

    //USING PSIG ENERGIES FOR WEIGHT
  doublevar teffac=teff/2.0;

  doublevar effenergy=0, effoldenergy=0;

  effenergy=pt.gprop.energy(0);
  effoldenergy=pt.gpast_energies[0].main_en;

  doublevar fbet=max(etr-effenergy, etr-effoldenergy);

  if(fbet > branchcut_stop) {
    teffac=0;
  }
  else if(fbet > branchcut_start) {
    teffac=teffac*(1.-(fbet-branchcut_start)
                   /(branchcut_stop-branchcut_start));
  }

  doublevar return_weight;
  return_weight=exp(teffac*(etr*2-effenergy-effoldenergy));

  return return_weight;
}

//----------------------------------------------------------------------
                                
struct Queue_element { 
  int from_node;
  int to_node;
  Queue_element() { } 
  Queue_element(int from, int to) { from_node=from; to_node=to; } 
};

struct Walker_sort { 
  int node;
  int index; //on the node
  int branch; //how many copies to make
  doublevar weight; 
};



int Sorndmc_method::calcBranch() { 
  int totwalkers=mpi_info.nprocs*nconfig;
  Array1 <doublevar> weights(totwalkers);
  Array1 <doublevar> my_weights(nconfig);
  
  for(int walker=0; walker < nconfig; walker++)
    my_weights(walker)=pts(walker).weight;
#ifdef USE_MPI
  MPI_Allgather(my_weights.v,nconfig, MPI_DOUBLE, weights.v,nconfig,MPI_DOUBLE, MPI_Comm_grp);
#else
  weights=my_weights;
#endif
  Array1 <int> my_branch(nconfig);
  Array1 <int> nwalkers(mpi_info.nprocs);
  nwalkers=0;
  int root=0;
  if(mpi_info.node==root) {  //this if/else clause may be refactored out

    Array1 <int> branch(totwalkers);
    //----Find which walkers branch/die
    //we do it on one node since otherwise we'll have different random numbers!
    //we'll assign the weight for each copy that will be produced
    //this is the core of the branching algorithm..
    //my homegrown algo, based on Umrigar, Nightingale, and Runge
    branch=-1;
    const doublevar split_threshold=1.8;
    for(int w=0; w< totwalkers; w++) { 
      if(weights(w) > split_threshold && branch(w)==-1) { 
        //find branching partner
        doublevar smallestwt=100;
        int smallest=-1;
        for(int w2=0; w2 < totwalkers; w2++) { 
          if(branch(w2)==-1 && w2!= w && weights(w2) < smallestwt) { 
            smallest=w2;
            smallestwt=weights(w2);
          }
        }
        if(smallest != -1) { 
          doublevar weight1=weights(w)/(weights(w)+weights(smallest));
          if(weight1+rng.ulec() >= 1.0) { 
            branch(w)=2;
            branch(smallest)=0;
            weights(w)+=weights(smallest);
            weights(w)/=2.0;
          }
          else { 
            branch(w)=0;
            branch(smallest)=2;
            weights(smallest)+=weights(w);
            weights(smallest)/=2.0;
          }
        }
      }
    }
    for(int w=0; w< totwalkers; w++) {
      if(branch(w)==-1) branch(w)=1;
    }
    //----end homegrown algo
    //count how many walkers each node will have
    //without balancing
    int walk=0;
    for(int n=0; n< mpi_info.nprocs; n++) { 
      for(int i=0; i< nconfig; i++) {
        nwalkers(n)+=branch(walk);
        walk++;
      }
      //cout << "nwalkers " << n << " " << nwalkers(n) << endl;
    }
    //now send nwalkers and which to branch to all processors
    for(int i=0; i< nconfig; i++) { 
      my_branch(i)=branch(i);
      my_weights(i)=weights(i);
    }
#ifdef USE_MPI
    MPI_Bcast(nwalkers.v, mpi_info.nprocs, MPI_INT, mpi_info.node, MPI_Comm_grp);
    for(int i=1; i< mpi_info.nprocs; i++) {
      MPI_Send(branch.v+i*nconfig,nconfig,MPI_INT,i,0,MPI_Comm_grp);
      MPI_Send(weights.v+i*nconfig, nconfig, MPI_DOUBLE, i,0,MPI_Comm_grp);
    }
#endif
               
  }
  else { 
#ifdef USE_MPI
    MPI_Bcast(nwalkers.v, mpi_info.nprocs, MPI_INT, root, MPI_Comm_grp);
    MPI_Status status;
    MPI_Recv(my_branch.v,nconfig, MPI_INT,root,0,MPI_Comm_grp, &status);
    MPI_Recv(my_weights.v,nconfig, MPI_DOUBLE,root,0,MPI_Comm_grp, &status);
#endif	
  }
  //--end if/else clause

  for(int i=0; i< nconfig; i++) { 
    pts(i).weight=my_weights(i);
  }
  
  //Now we all have my_branch and nwalkers..we need to figure out who 
  //needs to send walkers to whom--after this, nwalkers should be a flat array equal to 
  //nconfig(so don't try to use it for anything useful afterwards)
  vector <Queue_element> send_queue;
  int curr_needs_walker=0;
  int nnwalkers=nwalkers(mpi_info.node); //remember how many total we should have
  for(int i=0; i< mpi_info.nprocs; i++) { 
    while(nwalkers(i) > nconfig) {
      if(nwalkers(curr_needs_walker) < nconfig) { 
        nwalkers(curr_needs_walker)++;
        nwalkers(i)--;
        send_queue.push_back(Queue_element(i,curr_needs_walker));
        //cout << mpi_info.node << ":nwalkers " << nwalkers(i) << "  " << nwalkers(curr_needs_walker) << endl;
        //cout << mpi_info.node << ":send " << i << "  " << curr_needs_walker << endl;
      }
      else { 
        curr_needs_walker++;
      }
    }
  }
  
  for(int i=0; i< mpi_info.nprocs; i++) assert(nwalkers(i)==nconfig);
  int killsize=0;
  for(int i=0; i< nconfig; i++) {
    //cout << mpi_info.node << ":branch " << i << "  " << my_branch(i) << " weight " << pts(i).weight <<  endl;
    if(my_branch(i)==0) killsize++;
  }
  //cout << mpi_info.node << ": send queue= " << send_queue.size() << endl;
  //now do branching for the walkers that we get to keep
  Array1 <Rndmc_point> savepts=pts;      // Use Rndmc point, keeps all energies and propagates sign. Only weight gets modified during branching
  int curr=0; //what walker we're currently copying from
  int curr_copy=0; //what walker we're currently copying to
  while(curr_copy < min(nnwalkers,nconfig)) { 
    if(my_branch(curr)>0) { 
      //cout << mpi_info.node << ": copying " << curr << " to " << curr_copy << " branch " << my_branch(curr) << endl;
      my_branch(curr)--;
      pts(curr_copy)=savepts(curr);
      //pts(curr_copy).weight=1;
      curr_copy++;
    }
    else curr++;
  }
  
  //Finally, send or receive spillover walkers 
  if(nnwalkers > nconfig) { 
    vector<Queue_element>::iterator queue_pos=send_queue.begin();
    while(curr < nconfig) { 
      if(my_branch(curr) > 0) { 
        my_branch(curr)--;
        while(queue_pos->from_node != mpi_info.node) { 
          queue_pos++;
        }
        //cout << mpi_info.node << ":curr " << curr << " my_branch " << my_branch(curr) << endl;
        //cout << mpi_info.node << ":sending " << queue_pos->from_node << " to " << queue_pos->to_node << endl;
        savepts(curr).mpiSend(queue_pos->to_node);
        queue_pos++;
      }
      else curr++;
    }
    
  }
  else { //if nnwalkers == nconfig, then this will just get skipped immediately
    vector <Queue_element>::iterator queue_pos=send_queue.begin();
    while(curr_copy < nconfig) { 
      while(queue_pos->to_node != mpi_info.node) queue_pos++;
      //cout << mpi_info.node <<":receiving from " << queue_pos->from_node << " to " << curr_copy << endl;
      pts(curr_copy).mpiReceive(queue_pos->from_node);
      //pts(curr_copy).weight=1;
      curr_copy++;
      queue_pos++;
    }
  }
  return killsize;
  //exit(0);
}
//----------------------------------------------------------------------

void Rndmc_point::read(istream & is) { 
  config_pos.read(is);
  int filepos=is.tellg();
  string dum;
  is >> dum;
  if(!caseless_eq(dum, "SORNDMC") or !caseless_eq(dum, "DMC")) {
    is.seekg(filepos);
    return;
  }
  
  if (caseless_eq(dum, "SORNDMC")) {
    is >> dum; //the {
    is >> dum >> weight;
    is >> dum >> sign;
    is >> dum >> nodal_cross_age;
    
  }
  else if (caseless_eq(dum, "DMC")) {
    is >> dum; //the {
    is >> dum >> weight;
    is >> dum >> sign;
  }
}

void Rndmc_point::write(ostream & os) { 
  string indent="";
  config_pos.write(os);
  os << "SORNDMC { \n";
  //prop.write(indent,os);
  os << "weight " << weight<< endl;
  os << "sign " << sign << endl;
  os << "nodal_cross_age" << nodal_cross_age << endl;
  /*
  for(deque<Dmc_history>::iterator i=past_energies.begin(); 
      i!=past_energies.end(); i++) { 
    os << "past_energies { ";
    i->write(os);
    os << "}\n";    
  }
  for(deque<Dmc_history_avgrets>::iterator i=past_properties.begin(); 
      i!=past_properties.end(); i++) {
    os << "past_properties { ";
    i->write(os);
    os << "}\n";
  }
   */
  os << "}\n";
}
