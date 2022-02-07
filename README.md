Barycentric correction for the BL SETI data products SIGPROC filterbank files

Requirements:
-----------------
1. SIGPROC
2. TEMPO
3. gsl


Install:
------------------
1. Edit Makefile to connect SIGPROC and GSL
2. make clean
3. make

Use:
-----------------
> barycentre_seti "Input file" -verbose > "Output file"

Limitations:
-----------------
#1. At this point, only works with 32-bit SIGPROC filterbank file with positive foff
#2. Negative and positive relative velocity both fixes spectra-to-spectra 
#3. Only negative relative velocity corrects the individual spectra 
 
Testing:
-----------------
Testing for this code is inside directory "pipeline_tests"

Python notebook to create four filterbank files using setigen. 
nVelDec.fil : Negative Velocity which is decreasing spectra-to-spectra
nVelInc.fil : Negative Velocity which is increasing spectra-to-spectra
pVelDec.fil : Positive Velocity which is decreasing spectra-to-spectra
pVelInc.fil : Positive Velocity which is increasing spectra-to-spectra 

Ones these files are created using the python code, use following code to create four plots.

python3 barycentric_correction_code_check.py 


