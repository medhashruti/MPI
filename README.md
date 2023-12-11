# MPI

Configuration Properties : configuration 'Active(Debug)' ; platform 'x64'

Debugging>Command...
$(I_MPI_ONEAPI_ROOT)\bin\mpiexec.exe

Debugging>Command arguments...
-np 12 "$(TargetPath)"

Fortran>Additional Include Directories....
$(I_MPI_ONEAPI_ROOT)\include

Linker>Additional Library Directories....
$(I_MPI_ONEAPI_ROOT)\lib\debug;$(I_MPI_ONEAPI_ROOT)\lib\

Linker>Input>Additional Dependancies....
impi.lib libmpi_ilp64.lib
