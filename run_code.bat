@ echo off
REM :: Compile and get executable program
ifort module*.f90 -c /extend-source=132
ifort *.for *.f90 /exe:a.out /extend-source=132
ifort *.for *.f90 /exe:a.out /extend-source=132

REM :: Calculate the Love number and deformed Green's function of several depths through the For loop.
FOR %%d IN ( 01 03 ) DO (
  echo depth of source       > input.dat
  echo %%d                  >> input.dat
  echo depth of deformation >> input.dat
  echo 0.0                  >> input.dat
  echo filename of love     >> input.dat
  echo 'LOVE12.DAT'         >> input.dat 
  echo 'LOVE32.DAT'         >> input.dat
  echo 'LOVE220.DAT'        >> input.dat
  echo 'LOVE33.DAT'         >> input.dat
  echo 'model.dat'          >> input.dat
  
  start /wait a.out
  
REM :: Delete some unnecessary files.
  if exist *.mod    del /q *.mod
  if exist *.obj    del /q *.obj
  if exist *12.DAT  del /q *12.DAT
  if exist *32.DAT  del /q *32.DAT
  if exist *220.DAT del /q *220.DAT
  if exist model1_*.DAT del /q model1_*.DAT
  if exist model2_*.DAT del /q model2_*.DAT
  if exist model3_*.DAT del /q model3_*.DAT
  
  if not exist %%d (
    mkdir %%d
  )
  copy *.DAT %%d /Y 
)