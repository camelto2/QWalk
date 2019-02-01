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

#include "NodeSVM_method.h"
#include "qmc_io.h"
#include "Wf_return.h"
#include <dlib/svm.h>

void NodeSVM_method::read(vector<string> words,
                          unsigned int & pos,
                          Program_options & options) {

    ndim = 3;

    if(!readvalue(words, pos=0, readconfig, "READCONFIG"))
        error("Must provide config file for NodeSVM_method");

}

int NodeSVM_method::generateVariables(Program_options & options) {
    allocate(options.systemtext[0],sys);
    if (options.twftext.size() < 1)
        error("Need trial wavefunction");
    allocate(options.twftext[0],sys,wfdata);
    sys->generatePseudo(options.pseudotext,pseudo);
    sys->generateSample(sample);
    wfdata->generateWavefunction(wf);
    sample->attachObserver(wf);
    return 1;
}

void NodeSVM_method::run(Program_options & options, ostream & output) {
 
    generateVariables(options);

    ifstream is(readconfig.c_str());
    if (is) {
        is.close();
        read_configurations(readconfig,pts);
    }
    else {
        error("Config file not found");
    }

    totconfig = pts.GetDim(0);
    Wf_return wfval(wf->nfunc(),2);
    for (int walker = 0; walker < totconfig; walker++  ) {
        pts(walker).restorePos(sample);
        wf->updateVal(wfdata,sample);
        wf->getVal(wfdata,0,wfval);
    }
}

