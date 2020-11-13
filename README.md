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
1. At this point, only works with 32-bit SIGPROC filterbank file with positive foff
2. Negative and positive relative velocity both fixes spectra-to-spectra 
3. Only negative relative velocity corrects the individual spectra 
 

