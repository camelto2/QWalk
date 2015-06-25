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
#include <vector>
#include <fstream>
#include <cmath>
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include "wf_writer.h"
using namespace std;

int main(int argc, char ** argv) {
    string diracout_filename;
    if (argc != 2) {
        cout << "Error: Need to provide Dirac output filename" << endl;
	exit(1);
    }

    return 0;    
}
