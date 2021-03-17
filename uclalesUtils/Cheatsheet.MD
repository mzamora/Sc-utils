# Contents
1. [Compiling UCLALES](#compiling-uclales)
2. [Run DYCOMS RF-01](run-dycoms-rf-01)

# Compiling UCLALES

Get uclales

`git clone https://github.com/uclales/uclales.git`

Create config/default.cmake with the following content

` cd uclales`

` nano config/default.cmake`

Copy this into the editor and save it

```
set(CMAKE_Fortran_COMPILER "gfortran")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_1       "/usr/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/usr/lib/libnetcdf.a")
set(HDF5_LIB_1         "/usr/lib/x86_64-linux-gnu/libhdf5_hl.a")
set(HDF5_LIB_2         "/usr/lib/x86_64-linux-gnu/libhdf5.a")
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} -ldl m z curl)
```
Note that the previous config worked in our ubuntu 18.04(?) servers. My current config working on ubuntu 20.04 on my laptop uses mpifort from mpich, because the one from openmpi keeps throwing me an error (although it works ok with other software), and hdf5 libraries switched to serial versions. As always, be careful with paths:
```
set(CMAKE_Fortran_COMPILER "gfortran")
set(Fortran_COMPILER_WRAPPER "/usr/bin/mpifort.mpich")

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

set(NETCDF_INCLUDE_DIR "/usr/include")
set(NETCDF_LIB_1       "/usr/lib/x86_64-linux-gnu/libnetcdff.a")
set(NETCDF_LIB_2       "/usr/lib/x86_64-linux-gnu/libnetcdf.so")
set(HDF5_LIB_1         "/usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.a")
set(HDF5_LIB_2         "/usr/lib/x86_64-linux-gnu/libhdf5_serial.a")
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
```


Make a build directory and compile 
```
mkdir build
cd build
cmake -D MPI=TRUE ..
make -j4
```

Create a run folder
```
mkdir ../run
cd ../run
```

Link the uclales executable and datafiles folder (for full rad scheme)
```
ln -s ../build/uclales .
ln -s ../bin/datafiles/ .
```

Compile  and link programs for stitching
```
cd ../misc/synthesis/
gfortran reduceps.f90 -I /usr/include/ -lnetcdff -lnetcdf -lcurl -lz -lhdf5_fortran -lhdf5_hl -lhdf5 -o reduceps
gfortran reducets.f90 -I /usr/include/ -lnetcdff -lnetcdf -lcurl -lz -lhdf5_fortran -lhdf5_hl -lhdf5 -o reducets
cd ../../run/
ln -s ../misc/synthesis/reducets
ln -s ../misc/synthesis/reduceps
```
Similarly to the config file, if you have issues in newer systems, try these:
```
gfortran reduceps.f90 -I /usr/include/ -lnetcdff -lnetcdf -lcurl -lz -lhdf5_serial_hl -lhdf5_serial -o reduceps
gfortran reducets.f90 -I /usr/include/ -lnetcdff -lnetcdf -lcurl -lz -lhdf5_serial_hl -lhdf5_serial -o reducets

```


# Run DYCOMS-RF01
(Currently in the run folder) Copy the NAMELIST

`cp misc/initfiles/namelist_dycm01 NAMELIST`

To create the input sounding, first compile bstate.f95
```
cd ../misc/initfiles
gfortran bstate.f95 -I /usr/include/ -lnetcdff -lnetcdf -lcurl -lz -lhdf5_fortran -lhdf5_hl -lhdf5
mv a.out bstate
./bstate
mv sound_in ../../run/
```

Everything is ready to run stuff. Make sure that your mpiexec or mpirun are consistent with what you put in the config file (openMPI or mpich choice)
`mpiexec -n 8 ./uclales > log`

Stitch the files
```
./reducets $casename 4 2
./reduceps $casename 4 2
```
