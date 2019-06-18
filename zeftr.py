# Python wrapper class for zeft_recon code.
# This is just like zelrecon, but with an extra parameter.


import numpy  as np
import socket
import ctypes
import os

# Find out where we are.
fullpath = os.path.dirname(__file__)

# Some global variables, can be used by other codes
# calling ZelRecon as "defaults".
Rf = 15.0
ff = 0.76



class ZeftRecon:
    """
    A Python class to conviently "wrap" calls to the
    Zeldovich+EFT reconstruction code.
    This uses temporary files to pass the information around.
    """
    __author__ = "Martin White"
    __version__ = "1.0"
    __email__  = "mwhite@berkeley.edu"
    #
    def __call__(self,pkfile,ff=0,b1=0,b2=0,bs=0,Rf=-1,\
                 Apar=1,Aperp=1,alphad=0,alphan=0):
        """
        Runs the code and returns a NumPy array containing
        the data returned by the code.
        Note: as currently configured returns s^2 xi_ell, not  xi_ell.
        """
        ret = self.mylib.call_zeft_recon(\
          ctypes.c_char_p(pkfile.encode('utf-8')),\
          ctypes.c_double(ff),ctypes.c_double(b1),ctypes.c_double(b2),\
          ctypes.c_double(bs),ctypes.c_double(Rf),\
          ctypes.c_double(Apar),ctypes.c_double(Aperp),\
          ctypes.c_double(alphad),ctypes.c_double(alphan),\
          ctypes.c_char_p(self.tmpfn.encode('utf-8')))
        if (ret==0)&(os.path.isfile(self.tmpfn)):
            dd = np.loadtxt(self.tmpfn)
            os.remove(self.tmpfn)
        else:
            outstr = "ZeftRecon call failed with: "+pkfile+","+str(ff)+\
                     ","+str(b1)+","+str(b2)+","+str(bs)+","+str(Rf)+\
                     ","+str(Apar)+","+str(Aperp)+\
                     ","+str(alphad)+","+str(alphan)
            raise(RuntimeError,outstr)
            dd = None
        return(dd)
        #
    def __init__(self):
        """
        Initialize the class...very lightweight.
        """
        # Basic initialization, including a temporary file
        # whose name is based on the current host, PPID and PID.
        self.tmpfn = "zeft_recon_{:s}_{:d}_{:d}.txt".\
          format(socket.gethostname(),os.getppid(),os.getpid())
        #self.mylib = ctypes.CDLL(os.getcwd()+"/zeft_recon_ctypes.so")
        self.mylib = ctypes.CDLL(fullpath+"/zeft_recon_ctypes.so")
    #




def peak_background_bias(nu):
    """
    Returns the Lagrangian biases, (b1,b2), given nu.
    This is helpful if we want to make our basis set f, nu.
    """
    delc = 1.686
    a    = 0.707
    p    = 0.30
    anu2 = a*nu**2
    b1   = (anu2-1+2*p/(1+anu2**p))/delc
    b2   = (anu2**2-3*anu2+2*p*(2*anu2+2*p-1)/(1+anu2**p))/delc**2
    return( (b1,b2) )
    #
