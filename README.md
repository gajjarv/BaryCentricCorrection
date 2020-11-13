Barycentric correction for the BL SETI data products SIGPROC filterbank files

Requirements:
-----------------
SIGPROC
gsl

Install:
------------------
Connect installed version of SIGPROC 
> export SIGPROC="<path to libsigproc.a file>"
> make clean
> make

Use:
-----------------
> barycentre_seti <Input file> -verbose > <Output file>

Limitations:
-----------------
1. At this point, only works with 32-bit SIGPROC filterbank file
2. Negative and positive relative velocity both fixes spectra-to-spectra 
3. Only negative relative velocity corrects the individual spectra 
 

