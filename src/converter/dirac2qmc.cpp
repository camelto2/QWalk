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
#include <complex>
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

template <typename T>
void insertion_sort (vector <T> &list) {
    for(int j=1;j<list.size();j++) {
	T key=list[j];
	T i = j-1;
	while(i>-1 and list[i]>key) {
	    list[i+1]=list[i];
            i=i-1;
        }
        list[i+1]=key;
    }
}

void parse(string &s,vector <string> &parsed_string) {
    stringstream ss(s);
    string buf;
    while (ss >> buf)
	parsed_string.push_back(buf);
}

void skiplines(ifstream &is, int n) {
    string line;
    for (int i = 0; i < n; i++)
	getline(is,line);
}

class Orbital {
    public:
	int natoms;
	vector <string> label;
	vector < vector < vector < complex<double> > > > coeff;
	Orbital(vector <Gaussian_basis_set> &basis) {
	    natoms = basis.size();
	    string blank="";
	    for (int i = 0; i < natoms; i++) label.push_back(blank);
	    vector < vector < complex<double> > > tmp;
	    for (int at = 0; at < natoms; at++) 
		coeff.push_back(tmp);
	    for (int at = 0; at < natoms; at++) {
		vector < complex<double> > x;
		for (int i = 0; i < basis[at].types.size(); i++) {
		    coeff[at].push_back(x);
		}
	    }
	    for (int at = 0; at < natoms; at++) {
		complex<double> x(0.0,0.0);
		for (int l = 0; l < coeff[at].size(); l++) {
		    int n;
		    switch (l) {
			case 0: 
			    n=2;
			    break;
			case 1: 
			    n=6;
			    break;
			case 2:
			    n = 12;
			    break;
			case 3:
			    n = 20;
			    break;
			case 4:
			    n = 30;
			    break;
			default:
			    cout << "Unsupported basis function. Only s,p,d,f,g. Exiting" << endl;
			    exit(1);
		    }
		    for (int i = 0; i < n*basis[at].exponents[l].size(); i++) 
			coeff[at][l].push_back(x);
		}
	    }
	}
        void qwalk_ordering() {
	    for (int at = 0; at < natoms; at++) {
		for (int l = 2; l < coeff[at].size(); l++) {
		    for (int i = 0; i < coeff[at][l].size(); i++) {
			complex<double> tmp;
		        if (l == 2) { 
			    // DIRAC ordering: dxx dxy dxz dyy dyz dzz
			    // QWALK ordering: dxx dyy dzz dxy dxz dyz
			    if (i%12 == 2 || i%12 == 3) {
				tmp = coeff[at][l][i+4];
				coeff[at][l][i+4] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%12 == 4 || i%12 == 5) {
				tmp = coeff[at][l][i+6];
				coeff[at][l][i+6]=coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%12 == 8 || i%12 == 9) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2]=coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    
			}
			else if (l == 3) {
			    // DIRAC ordering: xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
			    // QWALK ordering: xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
			    if (i%20 == 2 || i%20 == 3) {
				tmp = coeff[at][l][i+10];
				coeff[at][l][i+10] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%20 == 4 || i%20 == 5) {
				tmp = coeff[at][l][i+14];
				coeff[at][l][i+14] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%20 == 6 || i%20 == 7) {
				tmp = coeff[at][l][i+6];
				coeff[at][l][i+6] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%20 == 8 || i%20 == 9) {
				tmp = coeff[at][l][i+10];
				coeff[at][l][i+10] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%20 == 10 || i%20 == 11) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%20 == 12 || i%20 == 13) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			}
			else if (l == 4) {
			    // DIRAC ordering: 400 310 301 220 211 202 130 121 112 103 040 031 022 013 004
			    // QWALK ordering: 400 040 004 310 301 130 031 103 013 220 202 022 211 121 112
			    if (i%30 == 2 || i%30 == 3) {
				tmp = coeff[at][l][i+18];
				coeff[at][l][i+18] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 4 || i%30 == 5) {
				tmp = coeff[at][l][i+24];
				coeff[at][l][i+24] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 6 || i%30 == 7) {
				tmp = coeff[at][l][i+14];
				coeff[at][l][i+14] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 8 || i%30 == 9) {
				tmp = coeff[at][l][i+20];
				coeff[at][l][i+20] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 10 || i%30 == 11) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 12 || i%30 == 13) {
				tmp = coeff[at][l][i+10];
				coeff[at][l][i+10] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 14 || i%30 == 15) {
				tmp = coeff[at][l][i+4];
				coeff[at][l][i+4] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 16 || i%30 == 17) {
				tmp = coeff[at][l][i+10];
				coeff[at][l][i+10] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 18 || i%30 == 19) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 20 || i%30 == 21) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 22 || i%30 == 23) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 24 || i%30 == 25) {
				tmp = coeff[at][l][i+4];
				coeff[at][l][i+4] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			    else if (i%30 == 26 || i%30 == 27) {
				tmp = coeff[at][l][i+2];
				coeff[at][l][i+2] = coeff[at][l][i];
				coeff[at][l][i] = tmp;
			    }
			}
			else {
			    cout << "Error: Only supports up to g functions. Exiting" << endl;
			    exit(1);
			}
		    }
		}
	    }
	}
	void kramers_pair(Orbital & old) {

	}
};

void read_dirac_mol(string & molfilename,
	            vector <Atom> & atoms,
		    vector <Gaussian_pseudo_writer> & arep,
		    vector <Gaussian_pseudo_writer> & sorep,
		    vector <Gaussian_basis_set> & basis);

void convert_arep_sorep_rrep(vector <Gaussian_pseudo_writer> & arep,
	                     vector <Gaussian_pseudo_writer> & sorep,
			     vector <Gaussian_pseudo_writer> & rrep);

void read_orb(vector <string> &orblines,
	      vector <Gaussian_basis_set> &basis,
	      vector <Atom> & atoms,
	      vector <Orbital> &orbs);

void read_dirac_out(string & outfilename,
	               vector <Atom> & atoms,
		       vector <Gaussian_pseudo_writer> & arep,
		       vector <Gaussian_basis_set> & basis,
		       vector <Orbital> & orbs);

void usage(const char * name) {
    cout << "usage: " << name << " <xx.inp> <xx.mol> <output root> " << endl;
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

void test_atoms(vector <Atom> & atoms) {
    for (int i = 0; i < atoms.size(); i++) {
	cout << "Atom: " << atoms[i].name << endl;
	cout << "      " << atoms[i].charge << endl;
	cout << "      " << "(" <<atoms[i].pos[0] << ", " 
	     << atoms[i].pos[1] << ", " << atoms[i].pos[2] << ") " << endl;
    }
}

void test_pseudo(vector <Gaussian_pseudo_writer> & psp) {
    for (int i = 0; i < psp.size(); i++) {
	cout << "ECP for Atom: " << psp[i].atomnum << endl;
	for (int j = 0; j < psp[i].exponents.size(); j++) {
	    cout << "    Channel: " << j << endl;
	    for (int k = 0; k < psp[i].exponents[j].size(); k++) {
		cout << "        " << psp[i].nvalue[j][k] << " " << psp[i].exponents[j][k] << " " << psp[i].coefficients[j][k] << endl;
	    }
	}
    }
}


void write_to_sys(string output_root, vector <Atom> & atoms, vector <Gaussian_pseudo_writer> & rrep);

void write_to_basis(string output_root, vector <Gaussian_basis_set> & basis, vector <Atom> & atoms);

void write_to_jast3(string output_root, vector <Atom> & atoms);

//######################################################################

int main(int argc, char ** argv) {

    if (argc != 4) usage(argv[0]);
    
    ifstream test_out; ifstream test_mol; ifstream test_inp;
    string dirac_out=string(argv[1])+"_"+string(argv[2])+".out";
    string dirac_inp=string(argv[1])+".inp";
    string dirac_mol=string(argv[2])+".mol";
    string output_root=string(argv[3]);
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
    vector <Gaussian_pseudo_writer> arep;
    vector <Gaussian_pseudo_writer> sorep;
    vector <Gaussian_pseudo_writer> rrep;
    vector <Gaussian_basis_set> basis;
    vector <Orbital> orbs;

    read_dirac_mol(dirac_mol,atoms,arep,sorep,basis);
    read_dirac_out(dirac_out,atoms,arep,basis,orbs);

    convert_arep_sorep_rrep(arep,sorep,rrep);

    write_to_sys(output_root, atoms, rrep);
    write_to_basis(output_root, basis, atoms);
    write_to_jast3(output_root, atoms);

    return 0;    
}

//######################################################################

void read_dirac_mol(string & molfilename,
	               vector <Atom> & atoms,
		       vector <Gaussian_pseudo_writer> & arep,
		       vector <Gaussian_pseudo_writer> & sorep,
		       vector <Gaussian_basis_set> & basis) {

    ifstream is(molfilename.c_str());
    if (!is) {
	cout << "Couldn't open " << molfilename << endl;
	exit(1);
    }

    string line;
    vector <string> words;
    int atom_num = 0;
    skiplines(is,3); // 1st 3 lines are comments
    getline(is,line);
    parse(line,words);
    int num_atoms = StringToNumber<int>(words[1]);
    for (int i = 0; i < num_atoms; i++) {
	atoms.push_back(Atom());
	atoms.back().basis = i;
    }
    while(getline(is,line)) {
	parse(line,words);
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
		getline(is,line); words.clear(); parse(line,words);
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
	//Read ECP
	if (line.find("ECP") != line.npos) {
            arep.push_back(Gaussian_pseudo_writer());
            sorep.push_back(Gaussian_pseudo_writer());
	    arep.back().atomnum = atom_num;
	    sorep.back().atomnum = atom_num;
	    int narep, nsorep;
	    narep = StringToNumber<int>(words[2]);
	    nsorep = StringToNumber<int>(words[3]);
	    //Allocate arep & sorep channels
	    for (int i = 0; i < narep; i++) {
                vector <double> tmp;
		vector <int> tmp2;
		arep.back().exponents.push_back(tmp);
		arep.back().coefficients.push_back(tmp);
		arep.back().nvalue.push_back(tmp2);
	    }
	    for (int i = 0; i < nsorep; i++) {
                vector <double> tmp;
		vector <int> tmp2;
		sorep.back().exponents.push_back(tmp);
		sorep.back().coefficients.push_back(tmp);
		sorep.back().nvalue.push_back(tmp2);
	    }
	    //Read AREP
	    for (int i = 0; i < narep; i++) {
                getline(is,line); words.clear(); parse(line,words);
		if (line.find("$") != line.npos || line.find("#") != line.npos) {
		    i -= 1;
		    continue;
		}
		else {
		    int nterms = StringToNumber<int>(words[0]);
		    //allocate number terms for given channel
		    for (int j = 0; j < nterms; j++) {
			getline(is,line); words.clear(); parse(line,words);
			int n = StringToNumber<int>(words[0]);
			double alpha = StringToNumber<double>(words[1]);
			double c = StringToNumber<double>(words[2]);
			arep.back().nvalue[i].push_back(n);
			arep.back().exponents[i].push_back(alpha);
			arep.back().coefficients[i].push_back(c);
		    }
		}
	    }
	    //Read SOREP
	    for (int i = 0; i < nsorep; i++) {
                getline(is,line); words.clear(); parse(line,words);
		if (line.find("$") != line.npos || line.find("#") != line.npos) {
		    i -= 1;
		    continue;
		}
		else {
		    int nterms = StringToNumber<int>(words[0]);
		    //allocate number terms for given channel
		    for (int j = 0; j < nterms; j++) {
			getline(is,line); words.clear(); parse(line,words);
			int n = StringToNumber<int>(words[0]);
			double alpha = StringToNumber<double>(words[1]);
			double c = StringToNumber<double>(words[2]);
			sorep.back().nvalue[i].push_back(n);
			sorep.back().exponents[i].push_back(alpha);
			sorep.back().coefficients[i].push_back(c);
		    }
		}
	    }
	    atom_num+=1;
	}
	words.clear();
    }
    is.close(); is.clear();
}

void convert_arep_sorep_rrep(vector <Gaussian_pseudo_writer> & arep,
	                     vector <Gaussian_pseudo_writer> & sorep,
			     vector <Gaussian_pseudo_writer> & rrep) {

    //Creating RREP for number of atoms. RREP has 2n-2 channels, where
    //n is the number of channels in AREP
    for (int at = 0; at < arep.size(); at++) {
	rrep.push_back(Gaussian_pseudo_writer());
	if (arep[at].nvalue.size() != 1) {
            for (int i = 0; i < (2*arep[at].nvalue.size()-2); i++) {
	        vector <double> tmp;
	        vector <int> tmp2;
	        rrep.back().nvalue.push_back(tmp2);
	        rrep.back().exponents.push_back(tmp);
	        rrep.back().coefficients.push_back(tmp);
	    }
	}
	else {
	     vector <double> tmp;
	     vector <int> tmp2;
	     rrep.back().nvalue.push_back(tmp2);
	     rrep.back().exponents.push_back(tmp);
	     rrep.back().coefficients.push_back(tmp);
	}
    }

    for (int at = 0; at < arep.size(); at++) {
	//Skip the first element in AREP...it is the local channel.
	//Add it last
	rrep[at].atomnum = arep[at].atomnum;
	for (int i = 1; i < arep[at].nvalue.size(); i++) {
	    double l = double(i-1); // because of how arep is stored
	    if (i == 1) {
		for (int j = 0; j < arep[at].coefficients[i].size(); j++) {
		    rrep[at].nvalue[i-1].push_back(arep[at].nvalue[i][j]);
		    rrep[at].exponents[i-1].push_back(arep[at].exponents[i][j]);
		    rrep[at].coefficients[i-1].push_back(arep[at].coefficients[i][j]);
		}
	    }
	    else {
		for (int j = 0; j < arep[at].coefficients[i].size(); j++) {
		    // l j=l+1/2
		    rrep[at].nvalue[i-1].push_back(arep[at].nvalue[i][j]);
		    rrep[at].exponents[i-1].push_back(arep[at].exponents[i][j]);
		    rrep[at].coefficients[i-1].push_back(arep[at].coefficients[i][j]+0.5*l*sorep[at].coefficients[i-2][j]);
		    // l j=l-1/2
		    rrep[at].nvalue[i].push_back(arep[at].nvalue[i][j]);
		    rrep[at].exponents[i].push_back(arep[at].exponents[i][j]);
		    rrep[at].coefficients[i].push_back(arep[at].coefficients[i][j]-0.5*(l+1.0)*sorep[at].coefficients[i-2][j]);
		}
	    }
	}
	//add local channel to end of RREP
	for (int i = 0; i < arep[at].nvalue[0].size(); i++) {
	    rrep[at].nvalue.back().push_back(arep[at].nvalue[0][i]);
	    rrep[at].exponents.back().push_back(arep[at].exponents[0][i]);
	    rrep[at].coefficients.back().push_back(arep[at].coefficients[0][i]);
	}
    }
    for (int at = 0; at < rrep.size(); at++) {
	rrep[at].label = arep[at].label;
	for (int i = 0; i < rrep[at].nvalue.size(); i++)  {
	    for (int j = 0; j < rrep[at].nvalue[i].size(); j++) 
		rrep[at].nvalue[i][j] -= 2;
        }
    }
}

void read_orb(vector <string> & orblines,
	      vector <Gaussian_basis_set> & basis,
	      vector <Atom> & atoms,
	      vector <Orbital> & orbs) {

    double snorm=sqrt(1.0/4.0/3.14159265359);
    double pnorm=snorm*sqrt(3.0);
    double dnorm=snorm*sqrt(15.0);
    double fnorm=snorm*sqrt(105.0);
    double gnorm=snorm*sqrt(315.0);

    vector <string> words;

    int norbs = 0;
    for (int i = 0; i < orblines.size(); i++) {
	if (orblines[i].find("Electronic") != orblines[i].npos) norbs++;
    }
    norbs *= 2;
    for (int mo = 0; mo < norbs; mo++) 
	orbs.push_back(Orbital(basis));
    for (int mo = 0; mo < norbs; mo++) {
	for (int at = 0; at < atoms.size(); at++) {
	    orbs[mo].label[at] = atoms[at].name;
	}
    }
/*
    for (int mo = 0; mo < norbs; mo++) {
	cout << "Orbital: " << mo << endl;
	for (int at = 0; at < orbs[mo].natoms; at++) {
	    cout << "    Atom: " << at << endl;
	    for (int l = 0; l < orbs[mo].coeff[at].size(); l++) {
		cout << "        Types: " << l << endl;
		for (int j = 0; j < orbs[mo].coeff[at][l].size(); j++) 
		    cout << "            " << orbs[mo].coeff[at][l][j] << endl;
	    }
	}
    }
*/

    int mo = 0;
    vector <int> sit;
    vector <int> pit;
    vector <int> dit;
    vector <int> fit;
    vector <int> git;
    vector <double> c(4);
    for (int at = 0; at < orbs[0].natoms; at++) {
        sit.push_back(0);
        pit.push_back(0);
        dit.push_back(0);
        fit.push_back(0);
        git.push_back(0);
    }

    for (int i = 0; i < orblines.size(); i++) {
	if (orblines[i].find("Electronic ") != orblines[i].npos) {
            for (int at = 0; at < orbs[0].natoms; at++) {
		sit[at] = 0;
		pit[at] = 0;
		dit[at] = 0;
		fit[at] = 0;
		git[at] = 0;
	    }
	    i += 2; 
	    while (orblines[i].find("Electronic") == orblines[i].npos) {
		if (orblines[i] == "") {
		    i++; 
		    break;
		}
		words.clear(); parse(orblines[i],words);
		c[0] = StringToNumber<double>(words[5]);
		c[1] = StringToNumber<double>(words[6]);
		c[2] = StringToNumber<double>(words[7]);
		c[3] = StringToNumber<double>(words[8]);
		for (int at = 0; at < orbs[mo].natoms; at++) {
		    if (words[2] == orbs[mo].label[at] && words[4] == "s") {
			orbs[mo].coeff[at][0][sit[at]] = complex<double>(c[0],c[1]);
			sit[at]++;
			orbs[mo].coeff[at][0][sit[at]] = complex<double>(c[2],c[3]);
			sit[at]++;
		    }
		    else if (words[2] == orbs[mo].label[at] && 
			    (words[4] == "px" || words[4] == "py" || words[4] == "pz")) {
			for (int j = 0; j < 4; j++) c[j] *= pnorm;
			if ((words[4] == "px" && pit[at]%6 == 0) ||
			    (words[4] == "py" && pit[at]%6 == 2) ||
			    (words[4] == "pz" && pit[at]%6 == 4)) {

			    orbs[mo].coeff[at][1][pit[at]] = complex<double>(c[0],c[1]);
			    pit[at]++;
			    orbs[mo].coeff[at][1][pit[at]] = complex<double>(c[2],c[3]);
			    pit[at]++;

			}
			else {
			    orbs[mo].coeff[at][1][pit[at]] = complex<double>(0.0,0.0);
			    pit[at]++;
			    orbs[mo].coeff[at][1][pit[at]] = complex<double>(0.0,0.0);
			    pit[at]++;
			    i--;
			}
		    }
		    else if (words[2] == orbs[mo].label[at] && (
			     words[4] == "dxx" || words[4] == "dxy" || words[4] == "dxz" ||
			     words[4] == "dyy" || words[4] == "dyz" || words[4] == "dzz")) {
			for (int j = 0; j < 4; j++) c[j] *= dnorm;
			if ((words[4] == "dxx" && dit[at]%12 == 0) ||
			    (words[4] == "dxy" && dit[at]%12 == 2) ||
			    (words[4] == "dxz" && dit[at]%12 == 4) ||
			    (words[4] == "dyy" && dit[at]%12 == 6) ||
			    (words[4] == "dyz" && dit[at]%12 == 8) ||
			    (words[4] == "dzz" && dit[at]%12 == 10)) {

			    orbs[mo].coeff[at][2][dit[at]] = complex<double>(c[0],c[1]);
			    dit[at]++;
			    orbs[mo].coeff[at][2][dit[at]] = complex<double>(c[2],c[3]);
			    dit[at]++;

			}
			else {
			    orbs[mo].coeff[at][2][dit[at]] = complex<double>(0.0,0.0);
			    dit[at]++;
			    orbs[mo].coeff[at][2][dit[at]] = complex<double>(0.0,0.0);
			    dit[at]++;
			    i--;
			}
		    }
		    else if (words[2] == orbs[mo].label[at] && (
			     words[4] == "fxxx" || words[4] == "fxxy" || words[4] == "fxxz" ||
			     words[4] == "fxyy" || words[4] == "fxyz" || words[4] == "fxzz" ||
			     words[4] == "fyyy" || words[4] == "fyyz" || words[4] == "fyzz" ||
			     words[4] == "fzzz" )) {
			for (int j = 0; j < 4; j++) c[j] *= fnorm;
			if ((words[4] == "fxxx" && fit[at]%20 == 0) ||
			    (words[4] == "fxxy" && fit[at]%20 == 2) ||
			    (words[4] == "fxxz" && fit[at]%20 == 4) ||
			    (words[4] == "fxyy" && fit[at]%20 == 6) ||
			    (words[4] == "fxyz" && fit[at]%20 == 8) ||
			    (words[4] == "fxzz" && fit[at]%20 == 10) ||
			    (words[4] == "fyyy" && fit[at]%20 == 12) ||
			    (words[4] == "fyyz" && fit[at]%20 == 14) ||
			    (words[4] == "fyzz" && fit[at]%20 == 16) ||
			    (words[4] == "fzzz" && fit[at]%20 == 18)) {

			    orbs[mo].coeff[at][3][fit[at]] = complex<double>(c[0],c[1]);
			    fit[at]++;
			    orbs[mo].coeff[at][3][fit[at]] = complex<double>(c[2],c[3]);
			    fit[at]++;
			    fit[at]++;

			}
			else {
			    orbs[mo].coeff[at][3][fit[at]] = complex<double>(0.0,0.0);
			    fit[at]++;
			    orbs[mo].coeff[at][3][fit[at]] = complex<double>(0.0,0.0);
			    fit[at]++;
			    i--;
			}
		    }
		    else if (words[2] == orbs[mo].label[at] && (
			     words[4] == "g400" || words[4] == "g310" || words[4] == "g301" ||
			     words[4] == "g220" || words[4] == "g211" || words[4] == "g202" ||
			     words[4] == "g130" || words[4] == "g121" || words[4] == "g112" ||
			     words[4] == "g103" || words[4] == "g040" || words[4] == "g031" ||
			     words[4] == "g022" || words[4] == "g013" || words[4] == "g004")) {
			for (int j = 0; j < 4; j++) c[j] *= gnorm;
			if ((words[4] == "g400" && git[at]%30 == 0) ||
			    (words[4] == "g310" && git[at]%30 == 2) ||
			    (words[4] == "g301" && git[at]%30 == 4) ||
			    (words[4] == "g220" && git[at]%30 == 6) ||
			    (words[4] == "g211" && git[at]%30 == 8) ||
			    (words[4] == "g202" && git[at]%30 == 10) ||
			    (words[4] == "g130" && git[at]%30 == 12) ||
			    (words[4] == "g121" && git[at]%30 == 14) ||
			    (words[4] == "g112" && git[at]%30 == 16) ||
			    (words[4] == "g103" && git[at]%30 == 18) ||
			    (words[4] == "g040" && git[at]%30 == 20) ||
			    (words[4] == "g031" && git[at]%30 == 22) ||
			    (words[4] == "g022" && git[at]%30 == 24) ||
			    (words[4] == "g013" && git[at]%30 == 26) ||
			    (words[4] == "g004" && git[at]%30 == 28)) {

			    orbs[mo].coeff[at][4][git[at]] = complex<double>(c[0],c[1]);
			    git[at]++;
			    orbs[mo].coeff[at][4][git[at]] = complex<double>(c[2],c[3]);
			    git[at]++;

			}
			else {
			    orbs[mo].coeff[at][4][git[at]] = complex<double>(0.0,0.0);
			    git[at]++;
			    orbs[mo].coeff[at][4][git[at]] = complex<double>(0.0,0.0);
			    git[at]++;
			    i--;
			}
		    }
		    else if (words[2] != orbs[mo].label[at]) continue;
		    else {
			cout << "Error converting orbitals" << endl;
			exit(1);
		    }
		}
		if (i == orblines.size() - 1) break;
		else i++;
	    }
	    i--;
	    orbs[mo].qwalk_ordering();
	    mo++;
	    orbs[mo].kramers_pair(orbs[mo-1]);
	    mo++;
	}
    }
}

void read_dirac_out(string & outfilename,
	               vector <Atom> & atoms,
		       vector <Gaussian_pseudo_writer> & arep,
		       vector < Gaussian_basis_set> & basis,
		       vector <Orbital> & orbs ) {

    ifstream is(outfilename.c_str());
    if (!is) {
	cout << "Couldn't open " << outfilename << endl;
	exit(1);
    }

    string line;
    vector <string> words;
    int num_atoms;

    vector <string> orblines;

    while(getline(is,line)) {
        parse(line,words);
	//Atom names and charges
        if (line.find("Atoms and basis") != line.npos) {
	    skiplines(is,3);
	    words.clear(); getline(is,line); parse(line,words);
	    num_atoms = StringToNumber<int>(words[4]);
            skiplines(is,3);
	    for (int i = 0; i < num_atoms; i++) {
		getline(is,line); words.clear(); parse(line,words);
		atoms[i].name = words[0];
		atoms[i].charge = StringToNumber<double>(words[2]);
		getline(is,line);
	    }
	}
	//Atom Positions
	if (line.find("Cartesian Coord") != line.npos) {
	    skiplines(is,5);
	    for (int i = 0; i < num_atoms; i++) {
                getline(is,line); words.clear(); parse(line,words);
		atoms[i].pos[0] = StringToNumber<double>(words[3]);
                getline(is,line); words.clear(); parse(line,words);
		atoms[i].pos[1] = StringToNumber<double>(words[2]);
                getline(is,line); words.clear(); parse(line,words);
		atoms[i].pos[2] = StringToNumber<double>(words[2]);
                skiplines(is,1);
	    }
	}
	// Orbitals
	if (line.find("***** Vector print *****") != line.npos) {
	    skiplines(is,12);
	    while (getline(is,line)) {
		if (line.find("*****************") != line.npos)
		    break;
		orblines.push_back(line);
	    }
	}
	// Slater Determinants
    }

    is.close(); is.clear();
    for (int at = 0; at < atoms.size(); at++)  
	arep[at].label = atoms[at].name;

    read_orb(orblines,basis,atoms,orbs);

}

void write_to_sys(string output_root, vector <Atom> & atoms, vector <Gaussian_pseudo_writer> & rrep) {

    int nelec = 0;
    for (int at = 0; at < atoms.size(); at++) 
	nelec += atoms[at].charge;

    ofstream sys;
    string sys_name=output_root+".sys";
    sys.open(sys_name.c_str());
    sys << "SYSTEM { MOLECULE" << endl;
    sys << "  SO_NSPIN { " <<  nelec << " 0 }" << endl;
    for (int at = 0; at < atoms.size(); at++) 
	atoms[at].print_atom(sys);
    sys << "}" << endl;
    sys << endl;
    for (int at = 0; at < rrep.size(); at++) 
	rrep[at].print_pseudo(sys);
    sys.close();

}

void write_to_basis(string output_root, vector <Gaussian_basis_set> & basis, vector <Atom> & atoms) {

    ofstream bas;
    string bas_name=output_root+".basis";
    bas.open(bas_name.c_str());
    for (int at = 0; at < basis.size(); at++) {
	bas << "BASIS { " << endl;
	bas << atoms[at].name << endl;
	bas << "SO_spline" << endl;
	bas << endl;
	bas << " GAMESS {" << endl;
	for ( int i = 0; i < basis[at].types.size(); i++) {
	    insertion_sort<double>(basis[at].exponents[i]);
	    for (int j = basis[at].exponents[i].size()-1; j>=0; j--) {
		bas << basis[at].types[i] << "   1 " << endl;
		bas << "  1  " << basis[at].exponents[i][j] << "  1 " << endl;
	    }
	}
        bas << " }"<< endl;
	bas << "} " << endl;
	bas << endl;
    }
    bas.close();
}

void write_to_jast3(string output_root, vector <Atom> & atoms) {

    ofstream jast;
    string jast_name=output_root+".jast3";
    jast.open(jast_name.c_str());
    jast << "JASTROW2" << endl;
    jast << "GROUP {" << endl;
    jast << " OPTIMIZEBASIS " << endl;
    jast << " EEBASIS { EE CUTOFF_CUSP GAMMA 24.0 CUSP 1.0 CUTOFF 15 } " << endl;
    jast << " EEBASIS { EE CUTOFF_CUSP GAMMA 24.0 CUSP 1.0 CUTOFF 15 } " << endl;
    jast << " TWOBODY_SPIN {  FREEZE " << endl;
    jast << "   LIKE_COEFFICIENTS { 0.5  0.0 } " << endl;
    jast << "   UNLIKE_COEFFICIENTS { 0.0 0.5 } " << endl;
    jast << " }" << endl;
    jast << "}" << endl;
    jast << "GROUP { " << endl;
    jast << " OPTIMIZEBASIS" << endl;
    jast << " EEBASIS { EE POLYPADE BETA0 0.5 NFUNC 4 RCUT 15 } " << endl;
    for (int at = 0; at < atoms.size(); at++) 
	jast << " EIBASIS { " << atoms[at].name << " POLYPADE BETA0 0.2 NFUNC 4 RCUT 15 }" << endl;
    jast << " ONEBODY {" << endl;
    for (int at = 0; at < atoms.size(); at++) 
	jast << "  COEFFICIENTS { " << atoms[at].name << " 0 0 0 0 }" << endl;
    jast << " }" << endl;
    jast << " TWOBODY {" << endl;
    jast << "  COEFFICIENTS { 0 0 0 0 } " << endl;
    jast << " }" << endl;
    jast << " THREEBODY {" << endl;
    for (int at = 0; at < atoms.size(); at++) 
	jast << "  COEFFICIENTS { " << atoms[at].name << " 0 0 0 0 0 0 0 0 0 0 0 0 } " << endl;
    jast << " }" << endl;
    jast << "}" << endl;

    jast.close();
}
