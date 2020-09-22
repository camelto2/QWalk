#!/usr/bin/python -tt
import sys
import string
import math
import os
from decimal import *

# Original by minyi, modified and expanded by Cody Melton
# Need dirac orbital information from output file by *PRIVEC
# just copy the relevant spinors to another file and use this to create 
# Qwalk orb file

# The printed spinor is only one of the kramer pairs, and it uses quaternion coefficients
# This generates the kramers pair by using the quaternion conjugate of the coefficient for
# the basis function
# also normalizes correctly for qwalk

def func(nIon,name,ns_list,np_list,nd_list,nf_list,ng_list,orbital,table):
        basis=list()
        nbasis=list()
        totbas=list()
        for i in range(nIon):
          ns=ns_list[i]
          np=np_list[i]
          nd=nd_list[i]
          nf=nf_list[i]
          ng=ng_list[i]
          nba=[0,ns,ns+np,ns+np*2,ns+np*3,ns+np*3+nd,ns+np*3+nd*2,ns+np*3+nd*3,ns+np*3+nd*4,ns+np*3+nd*5,ns+np*3+nd*6,ns+np*3+nd*6+nf,ns+np*3+nd*6+nf*2,ns+np*3+nd*6+nf*3,ns+np*3+nd*6+nf*4,ns+np*3+nd*6+nf*5,ns+np*3+nd*6+nf*6,ns+np*3+nd*6+nf*7,ns+np*3+nd*6+nf*8,ns+np*3+nd*6+nf*9,ns+np*3+nd*6+nf*10,ns+np*3+nd*6+nf*10+ng,ns+np*3+nd*6+nf*10+ng*2,ns+np*3+nd*6+nf*10+ng*3,ns+np*3+nd*6+nf*10+ng*4,ns+np*3+nd*6+nf*10+ng*5,ns+np*3+nd*6+nf*10+ng*6,ns+np*3+nd*6+nf*10+ng*7,ns+np*3+nd*6+nf*10+ng*8,ns+np*3+nd*6+nf*10+ng*9,ns+np*3+nd*6+nf*10+ng*10,ns+np*3+nd*6+nf*10+ng*11,ns+np*3+nd*6+nf*10+ng*12,ns+np*3+nd*6+nf*10+ng*13,ns+np*3+nd*6+nf*10+ng*14,ns+np*3+nd*6+nf*10+ng*15]
          nbasis.append(nba)
          totbas.append([ns,np,nd,nf,ng])
          bs=list()
          for j in range(ns+np*3+nd*6+nf*10+ng*15):
            bs.append([0.0,0.0,0.0,0.0])
          basis.append(bs)
        index=list()
        for i in range(nIon):
            tmp=list(nbasis[i])
            index.append(tmp)
        snorm=math.sqrt(1.0/4.0/3.14159265359)
        pnorm=snorm*math.sqrt(3.0)
        dnorm=snorm*math.sqrt(15.0)
        dnorm2=dnorm
        fnorm=snorm*math.sqrt(105.0)
        fnorm2=fnorm
        fnorm3=fnorm
        gnorm=snorm*math.sqrt(945.0)
        gnorm2=gnorm
        gnorm3=gnorm
        gnorm4=gnorm
        for element in orbital:
          temp=element.split()
          for i in range(nIon):
            if temp[2]==name[i]:#s px py pz dxx dyy dzz dxy dxz dyz fxxx fyyy fzzz fxxy fxxz 
                # fyyx fyyz fzzx fzzy fxyz gxxxx gyyyy gzzzz gxxxy gxxxz gyyyx gyyyz 
                # gzzzx gzzzy gxxyy gxxzz gyyzz gxxyz gyyxz gzzxy
              a=float(temp[5])
              b=float(temp[6])
              c=float(temp[7])
              d=float(temp[8])
              if temp[4]=='s':
                basis[i][index[i][0]][0]=a*snorm
                basis[i][index[i][0]][1]=b*snorm
                basis[i][index[i][0]][2]=c*snorm
                basis[i][index[i][0]][3]=d*snorm
                index[i][0]+=1
              if temp[4]=='px':
                basis[i][index[i][1]][0]=a*pnorm
                basis[i][index[i][1]][1]=b*pnorm
                basis[i][index[i][1]][2]=c*pnorm
                basis[i][index[i][1]][3]=d*pnorm
                index[i][1]+=1
              if temp[4]=='py':
                basis[i][index[i][2]][0]=a*pnorm
                basis[i][index[i][2]][1]=b*pnorm
                basis[i][index[i][2]][2]=c*pnorm
                basis[i][index[i][2]][3]=d*pnorm
                index[i][2]+=1
              if temp[4]=='pz':
                basis[i][index[i][3]][0]=a*pnorm
                basis[i][index[i][3]][1]=b*pnorm
                basis[i][index[i][3]][2]=c*pnorm
                basis[i][index[i][3]][3]=d*pnorm
                index[i][3]+=1
              if temp[4]=='dxx': 
                basis[i][index[i][4]][0]=a*dnorm
                basis[i][index[i][4]][1]=b*dnorm
                basis[i][index[i][4]][2]=c*dnorm
                basis[i][index[i][4]][3]=d*dnorm
                index[i][4]+=1
              if temp[4]=='dyy':
                basis[i][index[i][5]][0]=a*dnorm
                basis[i][index[i][5]][1]=b*dnorm
                basis[i][index[i][5]][2]=c*dnorm
                basis[i][index[i][5]][3]=d*dnorm
                index[i][5]+=1
              if temp[4]=='dzz':
                basis[i][index[i][6]][0]=a*dnorm
                basis[i][index[i][6]][1]=b*dnorm
                basis[i][index[i][6]][2]=c*dnorm
                basis[i][index[i][6]][3]=d*dnorm
                index[i][6]+=1
              if temp[4]=='dxy':
                basis[i][index[i][7]][0]=a*dnorm2
                basis[i][index[i][7]][1]=b*dnorm2
                basis[i][index[i][7]][2]=c*dnorm2
                basis[i][index[i][7]][3]=d*dnorm2
                index[i][7]+=1     
              if temp[4]=='dxz':
                basis[i][index[i][8]][0]=a*dnorm2
                basis[i][index[i][8]][1]=b*dnorm2
                basis[i][index[i][8]][2]=c*dnorm2
                basis[i][index[i][8]][3]=d*dnorm2
                index[i][8]+=1
              if temp[4]=='dyz':
                basis[i][index[i][9]][0]=a*dnorm2
                basis[i][index[i][9]][1]=b*dnorm2
                basis[i][index[i][9]][2]=c*dnorm2
                basis[i][index[i][9]][3]=d*dnorm2
                index[i][9]+=1
              if temp[4]=='fxxx':
                  basis[i][index[i][10]][0]=a*fnorm
                  basis[i][index[i][10]][1]=b*fnorm
                  basis[i][index[i][10]][2]=c*fnorm
                  basis[i][index[i][10]][3]=d*fnorm
                  index[i][10]+=1
              if temp[4]=='fyyy':
                  basis[i][index[i][11]][0]=a*fnorm
                  basis[i][index[i][11]][1]=b*fnorm
                  basis[i][index[i][11]][2]=c*fnorm
                  basis[i][index[i][11]][3]=d*fnorm
                  index[i][11]+=1
              if temp[4]=='fzzz':
                  basis[i][index[i][12]][0]=a*fnorm
                  basis[i][index[i][12]][1]=b*fnorm
                  basis[i][index[i][12]][2]=c*fnorm
                  basis[i][index[i][12]][3]=d*fnorm
                  index[i][12]+=1
              if temp[4]=='fxxy':
                  basis[i][index[i][13]][0]=a*fnorm2
                  basis[i][index[i][13]][1]=b*fnorm2
                  basis[i][index[i][13]][2]=c*fnorm2
                  basis[i][index[i][13]][3]=d*fnorm2
                  index[i][13]+=1
              if temp[4]=='fxxz':
                  basis[i][index[i][14]][0]=a*fnorm2
                  basis[i][index[i][14]][1]=b*fnorm2
                  basis[i][index[i][14]][2]=c*fnorm2
                  basis[i][index[i][14]][3]=d*fnorm2
                  index[i][14]+=1
              if temp[4]=='fxyy':
                  basis[i][index[i][15]][0]=a*fnorm2
                  basis[i][index[i][15]][1]=b*fnorm2
                  basis[i][index[i][15]][2]=c*fnorm2
                  basis[i][index[i][15]][3]=d*fnorm2
                  index[i][15]+=1
              if temp[4]=='fyyz':
                  basis[i][index[i][16]][0]=a*fnorm2
                  basis[i][index[i][16]][1]=b*fnorm2
                  basis[i][index[i][16]][2]=c*fnorm2
                  basis[i][index[i][16]][3]=d*fnorm2
                  index[i][16]+=1
              if temp[4]=='fxzz':
                  basis[i][index[i][17]][0]=a*fnorm2
                  basis[i][index[i][17]][1]=b*fnorm2
                  basis[i][index[i][17]][2]=c*fnorm2
                  basis[i][index[i][17]][3]=d*fnorm2
                  index[i][17]+=1
              if temp[4]=='fyzz':
                  basis[i][index[i][18]][0]=a*fnorm2
                  basis[i][index[i][18]][1]=b*fnorm2
                  basis[i][index[i][18]][2]=c*fnorm2
                  basis[i][index[i][18]][3]=d*fnorm2
                  index[i][18]+=1
              if temp[4]=='fxyz':
                  basis[i][index[i][19]][0]=a*fnorm3
                  basis[i][index[i][19]][1]=b*fnorm3
                  basis[i][index[i][19]][2]=c*fnorm3
                  basis[i][index[i][19]][3]=d*fnorm3
                  index[i][19]+=1
              if temp[4]=='g400':
                  basis[i][index[i][20]][0]=a*gnorm
                  basis[i][index[i][20]][1]=b*gnorm
                  basis[i][index[i][20]][2]=c*gnorm
                  basis[i][index[i][20]][3]=d*gnorm
                  index[i][20]+=1
              if temp[4]=='g040':
                  basis[i][index[i][21]][0]=a*gnorm
                  basis[i][index[i][21]][1]=b*gnorm
                  basis[i][index[i][21]][2]=c*gnorm
                  basis[i][index[i][21]][3]=d*gnorm
                  index[i][21]+=1
              if temp[4]=='g004':
                  basis[i][index[i][22]][0]=a*gnorm
                  basis[i][index[i][22]][1]=b*gnorm
                  basis[i][index[i][22]][2]=c*gnorm
                  basis[i][index[i][22]][3]=d*gnorm
                  index[i][22]+=1
              if temp[4]=='g310':
                  basis[i][index[i][23]][0]=a*gnorm2
                  basis[i][index[i][23]][1]=b*gnorm2
                  basis[i][index[i][23]][2]=c*gnorm2
                  basis[i][index[i][23]][3]=d*gnorm2
                  index[i][23]+=1
              if temp[4]=='g301':
                  basis[i][index[i][24]][0]=a*gnorm2
                  basis[i][index[i][24]][1]=b*gnorm2
                  basis[i][index[i][24]][2]=c*gnorm2
                  basis[i][index[i][24]][3]=d*gnorm2
                  index[i][24]+=1
              if temp[4]=='g130':
                  basis[i][index[i][25]][0]=a*gnorm2
                  basis[i][index[i][25]][1]=b*gnorm2
                  basis[i][index[i][25]][2]=c*gnorm2
                  basis[i][index[i][25]][3]=d*gnorm2
                  index[i][25]+=1
              if temp[4]=='g031':
                  basis[i][index[i][26]][0]=a*gnorm2
                  basis[i][index[i][26]][1]=b*gnorm2
                  basis[i][index[i][26]][2]=c*gnorm2
                  basis[i][index[i][26]][3]=d*gnorm2
                  index[i][26]+=1
              if temp[4]=='g103':
                  basis[i][index[i][27]][0]=a*gnorm2
                  basis[i][index[i][27]][1]=b*gnorm2
                  basis[i][index[i][27]][2]=c*gnorm2
                  basis[i][index[i][27]][3]=d*gnorm2
                  index[i][27]+=1
              if temp[4]=='g013':
                  basis[i][index[i][28]][0]=a*gnorm2
                  basis[i][index[i][28]][1]=b*gnorm2
                  basis[i][index[i][28]][2]=c*gnorm2
                  basis[i][index[i][28]][3]=d*gnorm2
                  index[i][28]+=1
              if temp[4]=='g220':
                  basis[i][index[i][29]][0]=a*gnorm3
                  basis[i][index[i][29]][1]=b*gnorm3
                  basis[i][index[i][29]][2]=c*gnorm3
                  basis[i][index[i][29]][3]=d*gnorm3
                  index[i][29]+=1
              if temp[4]=='g202':
                  basis[i][index[i][30]][0]=a*gnorm3
                  basis[i][index[i][30]][1]=b*gnorm3
                  basis[i][index[i][30]][2]=c*gnorm3
                  basis[i][index[i][30]][3]=d*gnorm3
                  index[i][30]+=1
              if temp[4]=='g022':
                  basis[i][index[i][31]][0]=a*gnorm3
                  basis[i][index[i][31]][1]=b*gnorm3
                  basis[i][index[i][31]][2]=c*gnorm3
                  basis[i][index[i][31]][3]=d*gnorm3
                  index[i][31]+=1
              if temp[4]=='g211':
                  basis[i][index[i][32]][0]=a*gnorm4
                  basis[i][index[i][32]][1]=b*gnorm4
                  basis[i][index[i][32]][2]=c*gnorm4
                  basis[i][index[i][32]][3]=d*gnorm4
                  index[i][32]+=1
              if temp[4]=='g121':
                  basis[i][index[i][33]][0]=a*gnorm4
                  basis[i][index[i][33]][1]=b*gnorm4
                  basis[i][index[i][33]][2]=c*gnorm4
                  basis[i][index[i][33]][3]=d*gnorm4
                  index[i][33]+=1
              if temp[4]=='g112':
                  basis[i][index[i][34]][0]=a*gnorm4
                  basis[i][index[i][34]][1]=b*gnorm4
                  basis[i][index[i][34]][2]=c*gnorm4
                  basis[i][index[i][34]][3]=d*gnorm4
                  index[i][34]+=1

        index=list()
        for i in range(nIon):
            tmp=list(nbasis[i])
            index.append(tmp)
        for i in range(nIon):        
           for j in range(totbas[i][0]):
                table.append([basis[i][index[i][0]+j][0],basis[i][index[i][0]+j][1] ])
                table.append([basis[i][index[i][0]+j][2],basis[i][index[i][0]+j][3] ])
           for j in range(totbas[i][1]):
                table.append([basis[i][index[i][1]+j][0],basis[i][index[i][1]+j][1] ])
                table.append([basis[i][index[i][1]+j][2],basis[i][index[i][1]+j][3] ])
                table.append([basis[i][index[i][2]+j][0],basis[i][index[i][2]+j][1] ])
                table.append([basis[i][index[i][2]+j][2],basis[i][index[i][2]+j][3] ])
                table.append([basis[i][index[i][3]+j][0],basis[i][index[i][3]+j][1] ])
                table.append([basis[i][index[i][3]+j][2],basis[i][index[i][3]+j][3] ])
           for j in range(totbas[i][2]):
                table.append([basis[i][index[i][4]+j][0],basis[i][index[i][4]+j][1] ])
                table.append([basis[i][index[i][4]+j][2],basis[i][index[i][4]+j][3] ])
                table.append([basis[i][index[i][5]+j][0],basis[i][index[i][5]+j][1] ])
                table.append([basis[i][index[i][5]+j][2],basis[i][index[i][5]+j][3] ])
                table.append([basis[i][index[i][6]+j][0],basis[i][index[i][6]+j][1] ])
                table.append([basis[i][index[i][6]+j][2],basis[i][index[i][6]+j][3] ])
                table.append([basis[i][index[i][7]+j][0],basis[i][index[i][7]+j][1] ])
                table.append([basis[i][index[i][7]+j][2],basis[i][index[i][7]+j][3] ])
                table.append([basis[i][index[i][8]+j][0],basis[i][index[i][8]+j][1] ])
                table.append([basis[i][index[i][8]+j][2],basis[i][index[i][8]+j][3] ])
                table.append([basis[i][index[i][9]+j][0],basis[i][index[i][9]+j][1] ])
                table.append([basis[i][index[i][9]+j][2],basis[i][index[i][9]+j][3] ])
           for j in range(totbas[i][3]):
                table.append([basis[i][index[i][10]+j][0],basis[i][index[i][10]+j][1] ])
                table.append([basis[i][index[i][10]+j][2],basis[i][index[i][10]+j][3] ])
                table.append([basis[i][index[i][11]+j][0],basis[i][index[i][11]+j][1] ])
                table.append([basis[i][index[i][11]+j][2],basis[i][index[i][11]+j][3] ])
                table.append([basis[i][index[i][12]+j][0],basis[i][index[i][12]+j][1] ])
                table.append([basis[i][index[i][12]+j][2],basis[i][index[i][12]+j][3] ])
                table.append([basis[i][index[i][13]+j][0],basis[i][index[i][13]+j][1] ])
                table.append([basis[i][index[i][13]+j][2],basis[i][index[i][13]+j][3] ])
                table.append([basis[i][index[i][14]+j][0],basis[i][index[i][14]+j][1] ])
                table.append([basis[i][index[i][14]+j][2],basis[i][index[i][14]+j][3] ])
                table.append([basis[i][index[i][15]+j][0],basis[i][index[i][15]+j][1] ])
                table.append([basis[i][index[i][15]+j][2],basis[i][index[i][15]+j][3] ])
                table.append([basis[i][index[i][16]+j][0],basis[i][index[i][16]+j][1] ])
                table.append([basis[i][index[i][16]+j][2],basis[i][index[i][16]+j][3] ])
                table.append([basis[i][index[i][17]+j][0],basis[i][index[i][17]+j][1] ])
                table.append([basis[i][index[i][17]+j][2],basis[i][index[i][17]+j][3] ])
                table.append([basis[i][index[i][18]+j][0],basis[i][index[i][18]+j][1] ])
                table.append([basis[i][index[i][18]+j][2],basis[i][index[i][18]+j][3] ])
                table.append([basis[i][index[i][19]+j][0],basis[i][index[i][19]+j][1] ])
                table.append([basis[i][index[i][19]+j][2],basis[i][index[i][19]+j][3] ])
           for j in range(totbas[i][4]):
                table.append([basis[i][index[i][20]+j][0],basis[i][index[i][20]+j][1] ])
                table.append([basis[i][index[i][20]+j][2],basis[i][index[i][20]+j][3] ])
                table.append([basis[i][index[i][21]+j][0],basis[i][index[i][21]+j][1] ])
                table.append([basis[i][index[i][21]+j][2],basis[i][index[i][21]+j][3] ])
                table.append([basis[i][index[i][22]+j][0],basis[i][index[i][22]+j][1] ])
                table.append([basis[i][index[i][22]+j][2],basis[i][index[i][22]+j][3] ])
                table.append([basis[i][index[i][23]+j][0],basis[i][index[i][23]+j][1] ])
                table.append([basis[i][index[i][23]+j][2],basis[i][index[i][23]+j][3] ])
                table.append([basis[i][index[i][24]+j][0],basis[i][index[i][24]+j][1] ])
                table.append([basis[i][index[i][24]+j][2],basis[i][index[i][24]+j][3] ])
                table.append([basis[i][index[i][25]+j][0],basis[i][index[i][25]+j][1] ])
                table.append([basis[i][index[i][25]+j][2],basis[i][index[i][25]+j][3] ])
                table.append([basis[i][index[i][26]+j][0],basis[i][index[i][26]+j][1] ])
                table.append([basis[i][index[i][26]+j][2],basis[i][index[i][26]+j][3] ])
                table.append([basis[i][index[i][27]+j][0],basis[i][index[i][27]+j][1] ])
                table.append([basis[i][index[i][27]+j][2],basis[i][index[i][27]+j][3] ])
                table.append([basis[i][index[i][28]+j][0],basis[i][index[i][28]+j][1] ])
                table.append([basis[i][index[i][28]+j][2],basis[i][index[i][28]+j][3] ])
                table.append([basis[i][index[i][29]+j][0],basis[i][index[i][29]+j][1] ])
                table.append([basis[i][index[i][29]+j][2],basis[i][index[i][29]+j][3] ])
                table.append([basis[i][index[i][30]+j][0],basis[i][index[i][30]+j][1] ])
                table.append([basis[i][index[i][30]+j][2],basis[i][index[i][30]+j][3] ])
                table.append([basis[i][index[i][31]+j][0],basis[i][index[i][31]+j][1] ])
                table.append([basis[i][index[i][31]+j][2],basis[i][index[i][31]+j][3] ])
                table.append([basis[i][index[i][32]+j][0],basis[i][index[i][32]+j][1] ])
                table.append([basis[i][index[i][32]+j][2],basis[i][index[i][32]+j][3] ])
                table.append([basis[i][index[i][33]+j][0],basis[i][index[i][33]+j][1] ])
                table.append([basis[i][index[i][33]+j][2],basis[i][index[i][33]+j][3] ])
                table.append([basis[i][index[i][34]+j][0],basis[i][index[i][34]+j][1] ])
                table.append([basis[i][index[i][34]+j][2],basis[i][index[i][34]+j][3] ])
        for i in range(nIon):
           for j in range(totbas[i][0]):
                table.append([-1.0*basis[i][index[i][0]+j][2],basis[i][index[i][0]+j][3] ])
                table.append([basis[i][index[i][0]+j][0],-1.0*basis[i][index[i][0]+j][1] ])
           for j in range(totbas[i][1]):
                table.append([-1.0*basis[i][index[i][1]+j][2],basis[i][index[i][1]+j][3] ])
                table.append([basis[i][index[i][1]+j][0],-1.0*basis[i][index[i][1]+j][1] ])
                table.append([-1.0*basis[i][index[i][2]+j][2],basis[i][index[i][2]+j][3] ])
                table.append([basis[i][index[i][2]+j][0],-1.0*basis[i][index[i][2]+j][1] ])
                table.append([-1.0*basis[i][index[i][3]+j][2],basis[i][index[i][3]+j][3] ])
                table.append([basis[i][index[i][3]+j][0],-1.0*basis[i][index[i][3]+j][1] ])
           for j in range(totbas[i][2]):
                table.append([-1.0*basis[i][index[i][4]+j][2],basis[i][index[i][4]+j][3] ])
                table.append([basis[i][index[i][4]+j][0],-1.0*basis[i][index[i][4]+j][1] ])
                table.append([-1.0*basis[i][index[i][5]+j][2],basis[i][index[i][5]+j][3] ])
                table.append([basis[i][index[i][5]+j][0],-1.0*basis[i][index[i][5]+j][1] ])
                table.append([-1.0*basis[i][index[i][6]+j][2],basis[i][index[i][6]+j][3] ])
                table.append([basis[i][index[i][6]+j][0],-1.0*basis[i][index[i][6]+j][1] ])
                table.append([-1.0*basis[i][index[i][7]+j][2],basis[i][index[i][7]+j][3] ])
                table.append([basis[i][index[i][7]+j][0],-1.0*basis[i][index[i][7]+j][1] ])
                table.append([-1.0*basis[i][index[i][8]+j][2],basis[i][index[i][8]+j][3] ])
                table.append([basis[i][index[i][8]+j][0],-1.0*basis[i][index[i][8]+j][1] ])
                table.append([-1.0*basis[i][index[i][9]+j][2],basis[i][index[i][9]+j][3] ])
                table.append([basis[i][index[i][9]+j][0],-1.0*basis[i][index[i][9]+j][1] ])
           for j in range(totbas[i][3]):
                table.append([-1.0*basis[i][index[i][10]+j][2],basis[i][index[i][10]+j][3] ])
                table.append([basis[i][index[i][10]+j][0],-1.0*basis[i][index[i][10]+j][1] ])
                table.append([-1.0*basis[i][index[i][11]+j][2],basis[i][index[i][11]+j][3] ])
                table.append([basis[i][index[i][11]+j][0],-1.0*basis[i][index[i][11]+j][1] ])
                table.append([-1.0*basis[i][index[i][12]+j][2],basis[i][index[i][12]+j][3] ])
                table.append([basis[i][index[i][12]+j][0],-1.0*basis[i][index[i][12]+j][1] ])
                table.append([-1.0*basis[i][index[i][13]+j][2],basis[i][index[i][13]+j][3] ])
                table.append([basis[i][index[i][13]+j][0],-1.0*basis[i][index[i][13]+j][1] ])
                table.append([-1.0*basis[i][index[i][14]+j][2],basis[i][index[i][14]+j][3] ])
                table.append([basis[i][index[i][14]+j][0],-1.0*basis[i][index[i][14]+j][1] ])
                table.append([-1.0*basis[i][index[i][15]+j][2],basis[i][index[i][15]+j][3] ])
                table.append([basis[i][index[i][15]+j][0],-1.0*basis[i][index[i][15]+j][1] ])
                table.append([-1.0*basis[i][index[i][16]+j][2],basis[i][index[i][16]+j][3] ])
                table.append([basis[i][index[i][16]+j][0],-1.0*basis[i][index[i][16]+j][1] ])
                table.append([-1.0*basis[i][index[i][17]+j][2],basis[i][index[i][17]+j][3] ])
                table.append([basis[i][index[i][17]+j][0],-1.0*basis[i][index[i][17]+j][1] ])
                table.append([-1.0*basis[i][index[i][18]+j][2],basis[i][index[i][18]+j][3] ])
                table.append([basis[i][index[i][18]+j][0],-1.0*basis[i][index[i][18]+j][1] ])
                table.append([-1.0*basis[i][index[i][19]+j][2],basis[i][index[i][19]+j][3] ])
                table.append([basis[i][index[i][19]+j][0],-1.0*basis[i][index[i][19]+j][1] ])
           for j in range(totbas[i][4]):
                table.append([-1.0*basis[i][index[i][20]+j][2],basis[i][index[i][20]+j][3] ])
                table.append([basis[i][index[i][20]+j][0],-1.0*basis[i][index[i][20]+j][1] ])
                table.append([-1.0*basis[i][index[i][21]+j][2],basis[i][index[i][21]+j][3] ])
                table.append([basis[i][index[i][21]+j][0],-1.0*basis[i][index[i][21]+j][1] ])
                table.append([-1.0*basis[i][index[i][22]+j][2],basis[i][index[i][22]+j][3] ])
                table.append([basis[i][index[i][22]+j][0],-1.0*basis[i][index[i][22]+j][1] ])
                table.append([-1.0*basis[i][index[i][23]+j][2],basis[i][index[i][23]+j][3] ])
                table.append([basis[i][index[i][23]+j][0],-1.0*basis[i][index[i][23]+j][1] ])
                table.append([-1.0*basis[i][index[i][24]+j][2],basis[i][index[i][24]+j][3] ])
                table.append([basis[i][index[i][24]+j][0],-1.0*basis[i][index[i][24]+j][1] ])
                table.append([-1.0*basis[i][index[i][25]+j][2],basis[i][index[i][25]+j][3] ])
                table.append([basis[i][index[i][25]+j][0],-1.0*basis[i][index[i][25]+j][1] ])
                table.append([-1.0*basis[i][index[i][26]+j][2],basis[i][index[i][26]+j][3] ])
                table.append([basis[i][index[i][26]+j][0],-1.0*basis[i][index[i][26]+j][1] ])
                table.append([-1.0*basis[i][index[i][27]+j][2],basis[i][index[i][27]+j][3] ])
                table.append([basis[i][index[i][27]+j][0],-1.0*basis[i][index[i][27]+j][1] ])
                table.append([-1.0*basis[i][index[i][28]+j][2],basis[i][index[i][28]+j][3] ])
                table.append([basis[i][index[i][28]+j][0],-1.0*basis[i][index[i][28]+j][1] ])
                table.append([-1.0*basis[i][index[i][29]+j][2],basis[i][index[i][29]+j][3] ])
                table.append([basis[i][index[i][29]+j][0],-1.0*basis[i][index[i][29]+j][1] ])
                table.append([-1.0*basis[i][index[i][30]+j][2],basis[i][index[i][30]+j][3] ])
                table.append([basis[i][index[i][30]+j][0],-1.0*basis[i][index[i][30]+j][1] ])
                table.append([-1.0*basis[i][index[i][31]+j][2],basis[i][index[i][31]+j][3] ])
                table.append([basis[i][index[i][31]+j][0],-1.0*basis[i][index[i][31]+j][1] ])
                table.append([-1.0*basis[i][index[i][32]+j][2],basis[i][index[i][32]+j][3] ])
                table.append([basis[i][index[i][32]+j][0],-1.0*basis[i][index[i][32]+j][1] ])
                table.append([-1.0*basis[i][index[i][33]+j][2],basis[i][index[i][33]+j][3] ])
                table.append([basis[i][index[i][33]+j][0],-1.0*basis[i][index[i][33]+j][1] ])
                table.append([-1.0*basis[i][index[i][34]+j][2],basis[i][index[i][34]+j][3] ])
                table.append([basis[i][index[i][34]+j][0],-1.0*basis[i][index[i][34]+j][1] ])

orbfilename=raw_input('Please input  orb filename: ')
my_file=open(orbfilename,"r+")
nIon=int(raw_input('Please input number of ions: '))
nMO=int(raw_input('Please input number of original MOs: '))
ns_list=list()
np_list=list()
nd_list=list()
nf_list=list()
ng_list=list()
name=list()
for i in range(nIon):
  print('This is the ion: ',i)
  name.append(raw_input('Please input the name of ion: '))
  ns_list.append(int(raw_input('Please input the number of s functions: ')))
  np_list.append(int(raw_input('Please input the number of p functions: ')))
  nd_list.append(int(raw_input('Please input the number of d functions: ')))
  nf_list.append(int(raw_input('Please input the number of f functions: ')))
  ng_list.append(int(raw_input('Please input the number of g functinos: ')))

with my_file as f_in:
    lines = list(line for line in (l.strip() for l in f_in) if line and line.split()[0]!='*')
orbitals=list()
xx=list()
for element in lines:
   if element=='====================================================':
      if xx:
         yy=xx[:]
         orbitals.append(yy)
         xx=[]
   else: 
      xx.append(element) 
if xx:
   yy=xx[:]
   orbitals.append(yy)
table=list()
for mo in range(nMO):
  func(nIon,name,ns_list,np_list,nd_list,nf_list,ng_list,orbitals[mo],table)
k=1
f=open("new.orb","w")
for i in range(2*nMO):
  for ion in range(nIon):
   for j in range(ns_list[ion]*1+np_list[ion]*3+nd_list[ion]*6+nf_list[ion]*10+ng_list[ion]*15):
     f.write("{:>8} {:>8} {:>8} {:>8}".format(str(i+1),str(j+1),str(ion+1),str(k)))
     k+=1
     f.write ( '\n ')
f.write('COEFFICIENTS')
f.write ( '\n ')
i=0
while(i<len(table)):
           if(i<len(table)):
             f.write('({}, {}, {}, {})\t'.format(str(table[i][0]),str(table[i][1]),str(table[i+1][0]),str(table[i+1][1])))
             i+=2
           if(i<len(table)):
             f.write('({}, {}, {}, {})\t'.format(str(table[i][0]),str(table[i][1]),str(table[i+1][0]),str(table[i+1][1])))
             i+=2
           if(i<len(table)):
             f.write('({}, {}, {}, {})\t'.format(str(table[i][0]),str(table[i][1]),str(table[i+1][0]),str(table[i+1][1])))
             i+=2
           if(i<len(table)):
             f.write('({}, {}, {}, {})\t'.format(str(table[i][0]),str(table[i][1]),str(table[i+1][0]),str(table[i+1][1])))
             i+=2
           if(i<len(table)):
             f.write('({}, {}, {}, {})\t'.format(str(table[i][0]),str(table[i][1]),str(table[i+1][0]),str(table[i+1][1])))
             i+=2
           f.write('\n')


