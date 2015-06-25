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
 
Author: Cody Melton
Date: 06/25/2015

*/
#include "converter.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <cmath>
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include "wf_writer.h"
using namespace std;

template <typename T>
T StringToNumber ( const string &Text ) 
{                              
    stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

template <typename T>
string NumberToString ( T Number )
{
    stringstream ss;
    ss << Number;
    return ss.str();
}

void parse(string &s,vector <string> &parsed_string) {
    stringstream ss(s);
    string buf;
    while (ss >> buf)
	parsed_string.push_back(buf);
}

void read_dirac_mol(string & molfilename,
	            vector <Atom> & atoms,
		    vector <Gaussian_pseudo_writer> & pseudo,
		    vector <Gaussian_basis_set> & basis);

void read_dirac_output(string & outputfilename,
	               vector <Atom> & atoms);

void usage(const char * name) {
    cout << "usage: " << name << " <xx.inp> <xx.mol> " << endl;
    exit(1);
}

void test_basis(vector <Gaussian_basis_set> & basis) {
    for (int i = 0; i < basis.size(); i++) {
	cout << "Basis Number: " << i << endl;
	for (int j = 0; j < basis[i].types.size(); j++) {
	    cout << "    Type: " << basis[i].types[j] << endl;
	    for (int k = 0; k < basis[i].exponents[j].size(); k++) 
		cout << "        " << basis[i].exponents[j][k] << endl;
	}
    }
}


//######################################################################

int main(int argc, char ** argv) {

    if (argc != 3) usage(argv[0]);
    
    ifstream test_out; ifstream test_mol; ifstream test_inp;
    string dirac_out=string(argv[1])+"_"+string(argv[2])+".out";
    string dirac_inp=string(argv[1])+".inp";
    string dirac_mol=string(argv[2])+".mol";
    test_inp.open(dirac_inp.c_str());
    test_mol.open(dirac_mol.c_str());
    test_out.open(dirac_out.c_str());
    if (!test_inp) {
        cerr << "Couldn't find " << dirac_inp << endl;
	exit(1);
    }
    else if (!test_mol) {
	cerr << "Couldn't find " << dirac_mol << endl;
	exit(1);
    }
    else if (!test_out) {
	cerr << "Couldn't find " << dirac_out << endl;
	exit(1);
    }
    test_inp.close(); test_inp.clear();
    test_mol.close(); test_mol.clear();
    test_out.close(); test_out.close();

    vector <Atom> atoms;
    vector <Gaussian_pseudo_writer> pseudo;
    vector <Gaussian_basis_set> basis;
    read_dirac_mol(dirac_mol,atoms,pseudo,basis);

    return 0;    
}

//######################################################################

void read_dirac_mol(string & molfilename,
	               vector <Atom> & atoms,
		       vector <Gaussian_pseudo_writer> & pseudo,
		       vector <Gaussian_basis_set> & basis) {

    ifstream is(molfilename.c_str());
    if (!is) {
	cout << "Couldn't open " << molfilename << endl;
	exit(1);
    }

    string line;
    vector <string> words;
    int linecount = 1;

    while(getline(is,line)) {
	parse(line,words);
	if (linecount == 4) { // Line that tells the number of atoms
	    for (int i = 0; i < StringToNumber<int>(words[1]); i++) {
		atoms.push_back(Atom());
	        atoms.back().basis = i;
	    }
	}
	// Read basis
	if (line.find("LARGE EXPLICIT") != line.npos) {
	    basis.push_back(Gaussian_basis_set());
	    int num_types=StringToNumber<int>(words[2]); //Number of types
	    //Allocate space for number of types and exponents
	    for (int i = 0; i < num_types; i++) { 
		vector <double> tmp;
		basis.back().exponents.push_back(tmp);
		switch (i) {
		    case 0: basis.back().types.push_back("S");
			break;
		    case 1: basis.back().types.push_back("P");
			break;
	   	    case 2: basis.back().types.push_back("6D");
			break;
		    case 3: basis.back().types.push_back("10F");
			break;
		    case 4: basis.back().types.push_back("15G");
			break;
		    default: cout << "Unsupported basis function" << endl;
			 exit(1);
		}
	    }
	    for (int i = 0; i < basis.back().types.size(); i++) {
		words.clear();
		getline(is,line);
		parse(line,words);
		if (line.find("#") != line.npos) {
		    i -= 1;
		    continue;
		}
		else if (line.find("f") != line.npos) {
    		    for (int j = 0; j < StringToNumber<int>(words[1]); j++) {
			getline(is,line);
			double expt = StringToNumber<double>(line);
			basis.back().exponents[i].push_back(expt);
		    }
		}
	    }
	}
	++linecount;
	words.clear();
    }
}
