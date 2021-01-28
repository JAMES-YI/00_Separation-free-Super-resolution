
This file provides overall picture of separation-free super-resolution project.

Created by JYI, 01/27/2021

------------------------------
ToDos
- clean codes and data
- experimental results demonstration of HMC and ANM
  with noiseless observations: (1) under different frequency separations; 
  (2) different sampling patterns; (3) different number of atoms; 
  (4) different sampling complexity; (5) different types of atoms
- experimental results demonstration of HMC and ANM 
  with noisy observations: (0) under different noise levels; 
  (1) under different frequency separations; 
  (2) different sampling patterns; (3) different number of atoms; 
  (4) different sampling complexity; (5) different types of atoms
- experimental results demonstration of HMC and ANM in phase transition figures:
  (0) under different noise levels; 
  (1) under different frequency separations; 
  (2) different sampling patterns; (3) different number of atoms; 
  (4) different sampling complexity; (5) different types of atoms
- 
- experimental results demonstration for successful recovery rate wrt separation
- experimental results demonstration for theorems

experimental results demonstration for atom separation lower bound 
  in [Xu et al., 2017] and [Tang et al., 2015]
- generate a sequence of frequencies ==> get atoms
- directly calculate atom distance by using the atoms corresponding 
  to the closest two frequencies
- directly calculate atom distance by using formula in Theorem 4


Reviewer's comments 1st round
- frequency separation in [Tang et al., 2015] is 2/((2N-2)*2pi) or 1/(N-1)
- In uniform sample at integer indices case, 
  why the frequency separation is 0.167 instead of 1/3 when 2N-2=2?
- The computation for frequency from atom separation can be problematic
- For calculating the bound from [Tang et al., 2015]: (1) directly set the 
  two frequencies to be apart by frequency separation df; 
  (2) construct two atoms, and calculate the atom distance
- For calculating the bound from [Xu et al., 2017]: (1) set the 
  frequencies to be on the grid; (2) for the atom set and construct the matrix;
  (3) calculate the atom bound from Theorem 4; (4) use interpolation relation between
  atom separation and frequency separation to estimate the frequency separation;
- When obtaining the interpolation relation: (1) one of the frequency should be fixed as 0
- why does Theorem 4 always give atom lower bound of 2?
- In non-uniform sampling (with missing samples) at integer indices, is the frequency separation
  lower bound 2/(N_max*2*pi) or 2/N_max?
- give details about how the bounds are computed


software development
- 

------------------------------


------------------------------


------------------------------

