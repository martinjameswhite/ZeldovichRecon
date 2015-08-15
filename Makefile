CC       = g++
OPTIONS  = -O -funroll-loops -fopenmp -shared -fPIC

all: zeft_recon_ctypes


zeft_recon_ctypes: zeft_recon_ctypes.cpp
        g++ $(OPTIONS) -o zeft_recon_ctypes.so zeft_recon_ctypes.cpp -lm


.PHONY : clean
clean:
        rm -f zeft_recon_ctypes.so
