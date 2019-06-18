#!/usr/bin/env python3
#
# Just test that the code runs and produces sensible output!
#
import numpy as np                                                      
import zeftr as Z                                                       
import matplotlib.pyplot as plt
# These are taken from the DarkSky simulation, but it could
# be almost anything for our purposes.
pk="pk_test.txt"
ff=0.50765
b1=1.0
b2=1.0
bs=1.0
Rf=15.00
# Call the code for pre- and post-reconstruction xi.
zel=Z.ZeftRecon()                                                       
raw=zel(pk,ff,b1,b2,bs,-1,1.0,1.0,0.0,0.0)
rec=zel(pk,ff,b1,b2,bs,Rf,1.0,1.0,0.0,0.0)
# and compare them in a figure.
fig,ax = plt.subplots(1,1)
ax.plot(raw[1:,0], raw[1:,1],'b:' ,label=r'Raw$(\ell=0)$')
ax.plot(raw[1:,0],-raw[1:,2],'b--',label=r'Raw$(\ell=2)$')
ax.plot(rec[1:,0], rec[1:,1],'r:' ,label=r'Rec$(\ell=0)$')
ax.plot(rec[1:,0],-rec[1:,2],'r--',label=r'Rec$(\ell=2)$')
ax.set_xlabel(r'$s\quad[h^{-1}{\rm Mpc}]$')
ax.set_ylabel(r'$i^\ell s^2\xi(s)\quad[h^{-2}{\rm Mpc}^2]$')
# tidy up and save the figure.
plt.tight_layout()
plt.savefig('test.pdf')
#
