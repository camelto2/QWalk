#! /usr/bin/env python

'''
Write QWalk files
'''
from functools import reduce
from pyscf import gto
from pyscf import lib
import numpy

def write_basis(basename,mol):
    import sys
    fout = open(basename+'.basis','w')
    for ia in mol._basis.keys():
        fout.write('BASIS {\n')
        fout.write(ia+'\n')
        fout.write('AOSPLINE\n\n')
        fout.write('GAMESS {\n')
        for ib in range(len(mol._basis[ia])):
            l = mol._basis[ia][ib][0]
            np = len(mol._basis[ia][ib])-1
            nc = len(mol._basis[ia][ib][1])-1
            for c in range(nc):
                if (l == 0):
                    fout.write('S  '+str(np)+'\n')
                elif (l == 1):
                    fout.write('P  '+str(np)+'\n')
                elif (l == 2):
                    fout.write('5D '+str(np)+'\n')
                elif (l == 3):
                    fout.write('7F_crystal '+str(np)+'\n')
                elif (l == 4):
                    fout.write('9G '+str(np)+'\n')
                else:
                    sys.stderr.write('QWalk cannot handle higher than G basis functions\n')
                    sys.stderr.write('Exiting!')
                    exit()
                for p in range(np):
                    exp = mol._basis[ia][ib][p+1][0]
                    coeff = mol._basis[ia][ib][p+1][c+1]
                    fout.write('  {}  {}  {}\n'.format(p+1,exp,coeff))
                    fout.write('}\n}\n\n')

def write_jast3(basename,mol):
    fout = open(basename+'.jast3','w')
    unique_atoms = []
    for ia in range(mol.natm):
        if mol.atom_pure_symbol(ia) not in unique_atoms:
            unique_atoms.append(mol.atom_pure_symbol(ia))
    fout.write('''JASTROW2
GROUP {
 OPTIMIZEBASIS
 EEBASIS { EE CUTOFF_CUSP GAMMA 24.0 CUSP 1.0 CUTOFF 7.5 }
 EEBASIS { EE CUTOFF_CUSP GAMMA 24.0 CUSP 1.0 CUTOFF 7.5 }
 TWOBODY_SPIN {  FREEZE
   LIKE_COEFFICIENTS { 0.25  0.0 }
   UNLIKE_COEFFICIENTS { 0.0 0.5 }
  }
}''')
    fout.write('\nGROUP {\n OPTIMIZEBASIS\n EEBASIS { EE POLYPADE BETA0 0.5 NFUNC 4 RCUT 7.5 }\n''')
    for ia in range(len(unique_atoms)):
        fout.write(' EIBASIS { '+unique_atoms[ia]+' POLYPADE BETA0 0.2 NFUNC 4 RCUT 7.5 }\n')
    fout.write(' ONEBODY {\n')
    for ia in range(len(unique_atoms)):
        fout.write('  COEFFICIENTS { '+unique_atoms[ia]+' 0 0 0 0 }\n')
    fout.write(' }\n TWOBODY {\n  COEFFICIENTS { 0 0 0 0 }\n }\n THREEBODY {\n')
    for ia in range(len(unique_atoms)):
        fout.write('COEFFICIENTS { '+unique_atoms[ia]+' 0 0 0 0 0 0 0 0 0 0 0 0 }\n')
    fout.write(' }\n}')


def count_terms(channel):
    tot=0
    for n in range(len(channel)):
        tot += len(channel[n])
    return tot


def print_pp_terms(channel,pp,fout):
    n = {}
    assert(len(pp[channel]) == 4)
    for i in range(4):
        n.update({i:pp[channel][i]})
    for i in range(4):
        for j in range(len(n[i])):
            fout.write(' '+str(i)+' '+str(n[i][j][0])+' '+str(n[i][j][1])+'\n')
    

def write_sys(basename,mol):
    import sys
    fout = open(basename+'.sys','w')
    fout.write('SYSTEM { \n MOLECULE\n')
    fout.write(' NSPIN { '+str(mol.nelec[0])+' '+str(mol.nelec[1])+' }\n')
    for ia in range(mol.natm):
        x,y,z = mol.atom_coord(ia)
        fout.write(' ATOM { '+mol.atom_pure_symbol(ia)+' '+str(mol.atom_charge(ia))+' COOR '+str(x)+' '+str(y)+' '+str(z)+' }\n')
    fout.write('}\n\n')
    ecp = mol.format_ecp(mol.ecp)
    pp = {}
    for at in ecp:
        fout.write('PSEUDO {\n '+at+'\n AIP 12\n')
        fout.write(' BASIS { '+at+'\n  RGAUSSIAN\n  OLDQMC {\n')
        nn = len(ecp[at][1])
        fout.write('0.0 '+str(nn)+'\n')
        for l in range(nn):
            if ecp[at][1][l][0] == -1:
                pp.update({'loc':ecp[at][1][l][1]})
            elif ecp[at][1][l][0] == 0:
                pp.update({'s':ecp[at][1][l][1]})
            elif ecp[at][1][l][0] == 1:
                pp.update({'p':ecp[at][1][l][1]})
            elif ecp[at][1][l][0] == 2:
                pp.update({'d':ecp[at][1][l][1]})
            elif ecp[at][1][l][0] == 3:
                pp.update({'f':ecp[at][1][l][1]})
            elif ecp[at][1][l][0] == 4:
                pp.update({'g':ecp[at][1][l][1]})
            else:
                sys.stderr.write("ECPs only up to G channel\n")
                sys.stderr.write("Exiting\n")
                exit()
        if 's' in pp.keys():
            fout.write(str(count_terms(pp['s']))+' ')
        if 'p' in pp.keys():
            fout.write(str(count_terms(pp['p']))+' ')
        if 'd' in pp.keys():
            fout.write(str(count_terms(pp['d']))+' ')
        if 'f' in pp.keys():
            fout.write(str(count_terms(pp['f']))+' ')
        if 'g' in pp.keys():
            fout.write(str(count_terms(pp['g']))+' ')
        if 'loc' in pp.keys():
            fout.write(str(count_terms(pp['loc']))+' ')

        fout.write('\n')

        if 's' in pp.keys():
            print_pp_terms('s',pp,fout)
        if 'p' in pp.keys():
            print_pp_terms('p',pp,fout)
        if 'd' in pp.keys():
            print_pp_terms('d',pp,fout)
        if 'f' in pp.keys():
            print_pp_terms('f',pp,fout)
        if 'g' in pp.keys():
            print_pp_terms('g',pp,fout)
        if 'loc' in pp.keys():
            print_pp_terms('loc',pp,fout)

        fout.write('  }\n }\n}\n\n')


def write_slater(basename,mol,mf):
    import sys
    import pyscf.scf
    import pyscf.dft
    fout = open(basename+'.slater','w')
    fout.write('SLATER\nORBITALS {\n CUTOFF_MO\n')
    nelec = mol.nelec[0]+mol.nelec[1]
    mag = max(nelec/20.0,1.0)
    fout.write(' MAGNIFY '+str(mag)+'\n')
    occa = []
    occb = []
    if isinstance(mf,pyscf.scf.uhf.UHF) or 'UHF' == mf.__class__.__name__\
    or isinstance(mf,pyscf.dft.uks.UKS) or 'UKS' == mf.__class__.__name__:
        norbs = len(mf.mo_occ[0])+len(mf.mo_occ[1])
        for i in range(len(mf.mo_occ[0])):
            if mf.mo_occ[0][i] > 0:
                occa.append(i+1)
        for i in range(len(mf.mo_occ[1])):
            if mf.mo_occ[1][i] > 0:
                occb.append(i+1+len(mf.mo_occ[0]))
    elif isinstance(mf,pyscf.scf.rhf.RHF) or 'RHF' == mf.__class__.__name__ \
    or isinstance(mf,pyscf.scf.rohf.ROHF) or 'ROHF' == mf.__class__.__name__\
    or isinstance(mf,pyscf.dft.rks.RKS) or 'RKS' == mf.__class__.__name__\
    or isinstance(mf,pyscf.dft.roks.ROKS) or 'ROKS' == mf.__class__.__name__:
        norbs = len(mf.mo_occ)
        for i in range(len(mf.mo_occ)):
            if mf.mo_occ[i] > 0:
                occa.append(i+1)
                if mf.mo_occ[i] == 2:
                    occb.append(i+1)
    else:
        sys.stderr.write('Only handles HF/KS calculations at the moment')
        sys.stderr.write('Exiting!')
        exit()
    fout.write(' NMO '+str(norbs)+'\n')
    fout.write(' ORBFILE '+basename+'.orb\n')
    fout.write(' INCLUDE '+basename+'.basis\n')
    fout.write(' CENTERS { USEATOMS }\n')
    fout.write('}\n\n')
    fout.write('DETWT { 1.0 }\n')
    fout.write('STATES {\n')
    for i in range(len(occa)):
        fout.write(' '+str(occa[i]))
    fout.write('\n')
    for i in range(len(occb)):
        fout.write(' '+str(occb[i]))
    fout.write('\n}')

def order_ao_index(mol):
# reorders the d, f, to qwalk orderings, i.e. 5D, 7D_crystal, 9G
# reorder d,f,g fucntion to
#  5D: D 0, D+1, D-1, D+2, D-2
#
#  7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
#
#  9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
    idx = []
    off = 0
    for ib in range(mol.nbas):
        l = mol.bas_angular(ib)
        for n in range(mol.bas_nctr(ib)):
            if l == 2:
                idx.extend([off+2,off+3,off+1,off+4,off+0])
            elif l == 3:
                idx.extend([off+3,off+4,off+2,off+5,off+1,off+6,off+0])
            elif l == 4:
                idx.extend([off+4,off+5,off+3,off+6,off+2,
                            off+7,off+1,off+8,off+0])
            elif l > 4:
                raise RuntimeError('l=5 is not supported')
            else:
                idx.extend(range(off,off+l*2+1))
            off += l * 2 + 1
    return idx

#just in case the norms are not right, we would need this. But I think they are fine as is
def ao_norm(mol):
    norm = []
    snorm = 1.0
    for ib in range(mol.nbas):
        l = mol.bas_angular(ib)
        for n in range(mol.bas_nctr(ib)):
            if l == 0:
                norm.extend([snorm])
            elif l == 1:
                norm.extend([snorm,snorm,snorm])
            elif l == 2:
                norm.extend([snorm,snorm,snorm,snorm,snorm])
            elif l == 3:
                norm.extend([snorm,snorm,snorm,snorm,snorm,snorm,snorm])
            elif l == 4:
                norm.extend([snorm,snorm,snorm,snorm,snorm,snorm,snorm,snorm,snorm])
            elif l > 4:
                raise RuntimeError('l=5 is not supported')
    return norm

def write_orb(basename,mol,mf):
    import pyscf.scf
    import pyscf.dft
    import sys
    fout = open(basename+'.orb','w')
    centers = []
    mo_coeff = []
    if isinstance(mf,pyscf.scf.uhf.UHF) or 'UHF' == mf.__class__.__name__\
    or isinstance(mf,pyscf.dft.uks.UKS) or 'UKS' == mf.__class__.__name__:
        mo_coeff.append(mf.mo_coeff[0])
        mo_coeff.append(mf.mo_coeff[1])
        mo_coeff = numpy.hstack(mo_coeff)
    elif isinstance(mf,pyscf.scf.rhf.RHF) or 'RHF' == mf.__class__.__name__ \
    or isinstance(mf,pyscf.scf.rohf.ROHF) or 'ROHF' == mf.__class__.__name__\
    or isinstance(mf,pyscf.dft.rks.RKS) or 'RKS' == mf.__class__.__name__\
    or isinstance(mf,pyscf.dft.roks.ROKS) or 'ROKS' == mf.__class__.__name__:
        mo_coeff = mf.mo_coeff
    else:
        sys.stderr.write('Only handles HF/KS calculations at the moment')
        sys.stderr.write('Exiting!')
        exit()

    for ib in range(mol.nbas):
        ia = mol.bas_atom(ib)
        l = mol.bas_angular(ib)
        c = mol._libcint_ctr_coeff(ib)
        np,nc = c.shape
        nsph = 2*l+1
        nd = nc*nsph
        centers.extend([ia+1]*(np*nsph))
    centers = numpy.hstack(centers)
    nprim, nmo = mo_coeff.shape

    count = 0
    for mo in range(nmo):
        for prim in range(nprim):
            fout.write('  {}  {}  {}  {}\n'.format(mo+1,prim+1,centers[prim],count+1))
            count += 1
    fout.write(' COEFFICIENTS')

    aoidx = order_ao_index(mol)
    aonorm = ao_norm(mol)
    assert(len(aoidx) == len(aonorm))
    count = 0
    for imo in range(nmo):
        for i,j in enumerate(aoidx):
            if count%5 == 0:
                fout.write('\n')
            fout.write('{} '.format(mo_coeff[j,imo]*aonorm[i]))
            count += 1

def write_qwalk(basename,mol,mf):
    write_basis(basename,mol)
    write_jast3(basename,mol)
    write_sys(basename,mol)
    write_slater(basename,mol,mf)
    write_orb(basename,mol,mf)


if __name__ == '__main__':
    import sys
    from pyscf import gto, dft

    mol = gto.Mole()
    mol.atom = ''' O 0. 0. 0. '''
    mol.basis = {'O': 'bfd-vdz'}
    mol.ecp = {'O': 'bfd-pp'}
    mol.symmetry = 'D2h'
    mol.spin = 2
    mol.charge = 0
    mol.build()

    mf = dft.ROKS(mol)
    mf.irrep_nelec = {'Ag': (1,1), 'B3u':(1,1),'B2u':(1,0),'B1u':(1,0)}
    mf.kernel()

    write_qwalk('O',mol,mf)
