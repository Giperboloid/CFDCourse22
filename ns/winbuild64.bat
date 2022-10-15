SETLOCAL

REM =============== System options
REM Define generator
SET CMGenerator="Visual Studio 16 2019"

REM Architechture: x64 or Win32
SET CMArch=x64

REM root boost library
SET CMBOOST_ROOT="c:/libs/boost_1_78_0"

REM root amgc library
SET CMAMGCL_ROOT="c:/libs/amgcl-1.4.2"

cmake -G %CMGenerator% -A %CMArch% ..^
	-DBOOST_ROOT=%CMBOOST_ROOT% ^
	-DAMGCL_ROOT=%CMAMGCL_ROOT%

ENDLOCAL

