# Python wrapper class for zelrecon code.


import numpy  as np
import socket
import ctypes
import os

# Some global variables, can be used by other codes
# calling ZelRecon as "defaults".
Rf = 15.0
ff = 0.76



class ZelRecon:
    """
    ZelRecon:
    A Python class to conviently "wrap" calls to the
    Zeldovich reconstruction code.
    This uses temporary files to pass the information around.
    """
    __author__ = "Martin White"
    __version__ = "1.0"
    __email__  = "mwhite@berkeley.edu"
    #
    def __call__(self,pkfile,ff,F1,F2,Rf,Apar,Aperp):
        """
        __call__(self,pkfile,ff,F1,F2,Rf,Apar,Aperp):
        Runs the code and returns a NumPy array containing
        the data returned by the code.
        Note: as currently configured returns s^2 xi_ell, not  xi_ell.
        """
        ret = self.mylib.call_zelrecon(ctypes.c_char_p(pkfile),\
          ctypes.c_double(ff),ctypes.c_double(F1),ctypes.c_double(F2),\
          ctypes.c_double(Rf),ctypes.c_double(Apar),ctypes.c_double(Aperp),\
          ctypes.c_char_p(self.tmpfn))
        if (ret==0)&(os.path.isfile(self.tmpfn)):
            dd = np.loadtxt(self.tmpfn)
            os.remove(self.tmpfn)
        else:
            outstr = "ZelRecon call failed with: "+pkfile+","+str(ff)+\
                     ","+str(F1)+","+str(F2)+","+str(Rf)+\
                     ","+str(Apar)+","+str(Aperp)
            raise RuntimeError,outstr
            dd = None
        return(dd)
        #
    def __init__(self):
        """
        __init__(self):
        """
        # Basic initialization, including a temporary file
        # whose name is based on the current host, PPID and PID.
        self.tmpfn = "zelrecon_%s_%d_%d.txt"%\
          (socket.gethostname(),os.getppid(),os.getpid())
        self.mylib = ctypes.CDLL(os.getcwd()+"/zelrecon_ctypes.so")
    #




def peak_background_bias(nu):
    """
    peak_background_bias(nu):
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
