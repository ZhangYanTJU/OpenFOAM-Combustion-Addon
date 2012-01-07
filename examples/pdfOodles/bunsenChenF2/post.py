#!/usr/bin/python

import os
from numpy import *

UMEAN='UMean'
RMEAN='UPrime2Mean'
TMEAN='TMean'
TVAR='TPrime2Mean'
U0=50.0
Tb=2248.0
Tu=300.0
D=0.012
k0=10.8

sample=False

sampleDictHeader="""
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    \"system\";
    object      sampleDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
interpolationScheme cellPointFace;
setFormat     raw;
surfaceFormat vtk;
surfaces
(
);
"""


if (sample):
 sdfile=open('system/sampleDict', 'w')
 sdfile.write(sampleDictHeader)
 sdfile.write("sets (\n")
 for XH in [25, 45, 65, 85, 105, 125]:
  x=0.1*XH*D
  sdfile.write("xToD_%d { type uniform; axis y;\
 start (%f 0.0 0.000); end (%f 0.024 0.000); nPoints 50; }\n" % (XH, x, x))
 sdfile.write("ax { type uniform; axis x; start (0 0 0); end (0.144 0 0); nPoints 50; }\n")
 sdfile.write(");\n")
 sdfile.write("fields (%s %s %s %s);\n" % (UMEAN, TMEAN, RMEAN, TVAR));
 sdfile.close()
 
 os.system('sample -latestTime')

SDIR=map(float, os.listdir('sets'))
SDIR.sort()
DIR=SDIR[-1]

import matplotlib
matplotlib.rcParams['font.size'] = 8
import matplotlib.pyplot as plt

DATA=os.environ['VGLDATEN']+"/bunsenChen"
print "Using data from "+DATA

plt.figure()
plt.subplots_adjust(wspace=0.5, hspace=0.6)
i=1
heights=[25, 45, 65, 85, 105]
rows=ceil(len(heights)/2.0)
for XH in heights:
 xd=0.1*XH
 T=loadtxt("sets/%s/xToD_%d_%s_%s.xy"%(DIR, XH, TMEAN, TVAR))
 c=loadtxt(DATA+"/c/cF2%d"%XH, skiprows=1)
 cp=c.take(where(c[:,0]>=0)[0], axis=0)
 crms=loadtxt(DATA+"/crms/crmsF2%d"%XH, skiprows=1)
 crmsp=crms.take(where(crms[:,0]>=0)[0], axis=0)
 plt.subplot(rows,2,i)
 plt.title('$x/D=%1.1f$'%xd)
 i+=1
 plt.grid(True)
 plt.xlabel('$r/D$')
 plt.ylabel('$\\langle c \\rangle$')
 l1=plt.plot(-T[:,0]/D, (T[:,1]-Tu)/(Tb-Tu), "r-")
 l2=plt.plot(-cp[:,0], cp[:,1], 'ro')
 plt.twinx()
 plt.ylabel('$\\sqrt{\\langle c\'^2 \\rangle}$')
 l3=plt.plot(T[:,0]/D, sqrt(T[:,2])/(Tb-Tu), "b-")
 l4=plt.plot(crmsp[:,0], crmsp[:,1], 'b+')
 if i==len(heights):
  plt.legend((l1, l2, l3, l4), 
  ('$\\langle c \\rangle$ (LES)', '$\\langle c \\rangle$ (Experiment)', 
   '$\\sqrt{\\langle c\'^2 \\rangle}$ (LES)', '$\\sqrt{\\langle c\'^2 \\rangle}$ (Experiment)'),
   loc=(0.1, -1.6))

plt.savefig('temp.pdf')



plt.figure()
plt.subplots_adjust(wspace=0.5, hspace=0.6)
i=1
heights=[25, 45, 65, 85, 105]
rows=ceil(len(heights)/2.0)
for XH in heights:
 xd=0.1*XH
 U=loadtxt("sets/%s/xToD_%d_%s.xy"%(DIR, XH, UMEAN))
 R=loadtxt("sets/%s/xToD_%d_%s.xy"%(DIR, XH, RMEAN))
 uexp=loadtxt(DATA+"/U/U%dF2H"%XH, skiprows=1)
 uexpp=uexp.take(where(uexp[:,0]>=0)[0], axis=0)
 Uc=loadtxt("cold/sets/1/xToD_%d_Umean.xy"%(XH))
 Rc=loadtxt("cold/sets/1/xToD_%d_R.xy"%(XH))
 uexpc=loadtxt(DATA+"/U/U%dF2C"%XH, skiprows=1)
 uexppc=uexpc.take(where(uexpc[:,0]>=0)[0], axis=0)
 kexp=loadtxt(DATA+"/K/K%dF2H"%XH, skiprows=1)
 kexpp=kexp.take(where(kexp[:,0]>=0)[0], axis=0)
 plt.subplot(rows,2,i)
 plt.title('$x/D=%1.1f$'%xd)
 i+=1
 plt.grid(True)
 plt.xlabel('$r/D$')
 plt.ylabel('$\\langle U_{ax} \\rangle / U_0$')
 l1=plt.plot(-U[:,0]/D, U[:,1]/U0, "b-")
 l2=plt.plot(-uexp[:,0], uexp[:,1], 'b+')
 l5=plt.plot(-Uc[:,0]/D, Uc[:,1]/U0, "g-")
 l6=plt.plot(-uexpc[:,0], uexpc[:,1], 'gx')
 plt.twinx()
 plt.ylabel('$\\langle k \\rangle / k_0$')
 kLES=0.5*(R[:,1]+R[:,5]+R[:,9])
 #kLES=0.5*(R[:,1]+2.*R[:,9])
 l3=plt.plot(R[:,0]/D, kLES/k0, "r-")
 l4=plt.plot(kexpp[:,0], kexpp[:,1], 'ro')
 if i==len(heights):
  plt.legend((l1, l2, l5, l6, l3, l4), 
  ('$\\langle U_{ax} \\rangle / U_0$ (LES, reacting)', 
   '$\\langle U_{ax} \\rangle / U_0$ (Experiment, reacting)', 
   '$\\langle U_{ax} \\rangle / U_0$ (LES, isothermal)', 
   '$\\langle U_{ax} \\rangle / U_0$ (Experiment, isothermal)', 
   '$\\langle k \\rangle / k_0$ (LES, reacting)', 
   '$\\langle k \\rangle / k_0$ (Experiment, reacting)'),
   loc=(-0.1, -1.7))

plt.savefig('velocity.pdf')
