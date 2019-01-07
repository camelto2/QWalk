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

#ifndef NODESVM_METHOD_H_INCLUDED
#define NODESVM_METHOD_H_INCLUDED

#include "Qmc_method.h"
#include "Program_options.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "System.h"
#include "Sample_point.h"
#include "Pseudopotential.h"

class NodeSVM_method: public Qmc_method {
    public:
        NodeSVM_method() {
            wf = NULL;
            wfdata = NULL;
            sys = NULL;
            sample = NULL;
            pseudo = NULL;
        }
        ~NodeSVM_method() {
            if (wf) {
                delete wf;
                wf = NULL;
            }
            if (wfdata) {
                delete wfdata;
                wfdata = NULL;
            }
            if (sys) {
                delete sys;
                sys = NULL;
            }
            if (sample) {
                delete sample;
                sample = NULL;
            }
            if (pseudo) {
                delete pseudo;
                pseudo = NULL;
            }
        }
        void read(vector <string> words,
                          unsigned int & pos,
                          Program_options & options);
        void run(Program_options & options, ostream & output);
        int generateVariables(Program_options & options);
    private:
        int ndim;
        string readconfig;
        int totconfig;
        Array1 <Config_save_point> pts;
        Wavefunction * wf;
        Wavefunction_data * wfdata;
        System * sys;
        Sample_point * sample;
        Pseudopotential * pseudo;
};

#endif
