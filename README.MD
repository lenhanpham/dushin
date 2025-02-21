### Dushin 2.0 supports all gaussian versions from g94 to g16

- This is the new version of dushin program for decomposition of reorganization energy. 
- The program supports all gaussian versions from 94 to 16: both fchk and log files
- This new version can read in any input names e.g.  input.in or input.inp or acd.input 

#### Compilation:

* The current version supports both gfortran and ifort compilers.  
* cd to the dushin dir src 
* module load intel-compiler/2021.1.1 or load gcc
* make changes to config.mk to select corresponding compiler. The default is gfortran 
* make all 
* make clean to clean up the src if errors happen

#### Run the program

Please have a look at the example to prepare input files

* The program needs: log outputs from Gaussian, fchk obtained using formchk input.chk 
* An example input is provided in the examples directory 

#### Leave comments in issues to have more support from me.

I do not have free time to explain how to use Dushin program to calculate internal reorganization energy from harmonic vibrational modes. Leave comments in Issues, and there will be support.

