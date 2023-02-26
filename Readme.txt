This code is developed based on the calculation code of seismic deformation by Professor Wenke Sun.
The main purpose of the current code is to calculate the surface deformation caused by a 1 kg point source inside the PREM Earth Model, including displacement and gravity.

The main way to use it is to run "run_code.bat". This is a batch file.
It will call the intel Fortran compiler in the system and compile the Fortran source files to get a.out.
Then it runs a.out to calculate the Love number and Green's function due to some source.

After about 5~10 minutes, the program will finish running and you can get several folders.
The description of the files in a folder nambed by source depth is as follows:

LOVE33.DAT contains y-numbers of all degrees.
The unit of each physical quantity is consistent with that in the differential equation.
# n y1 y2 y3 y4 y5 y6 y1t,y2t.
# Note y1t=0 and y2t=0 here.

DSLOV33 is similar as LOVE33.DAT but used just a linear interpolation.
# n y1 y2 y3 y4 y5 y6 y1t,y2t.
# Note y1t=0 and y2t=0 here (Ignore them).

GRNFN33.DAT is the Green's function as the summation of Love-numbers due to 1 kg loading. 
# Calculate each physical quantity in the spherical coordinate system (r,theta,phi).
# Four columns are
# theta(°) u_r(m/kg)*a*theta*1e12 u_theta(m/kg)*a*theta*1e12 gE(m/s^2/kg)*a*theta*1e18

input.dat
# This this the input file for a.out
# The second row is a value of source depth.

model.dat
# This is the parameters of the PREM model and it contains the explanation of each column.


#########################################################################################
Contact information:
He Tang, hetang@ucas.ac.cn
University of Chinese Academy of sciences.
