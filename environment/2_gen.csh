#!/bin/tcsh

########################################################################
# PREFIX and non-system programs/libraries
########################################################################

### prefix area
setenv PREFIX /group/clas/builds/centos7/trunk

### non-system builds of programs and libraries
setenv GCC /apps/gcc/5.3.0
setenv ROOT /apps/root/5.34.36
setenv CERN /apps/cernlib/x86_64_rhel7/2005
setenv PYTHON /apps/python/2.7.12
setenv SCONS /apps/scons
setenv BOOST /group/clas/boost/boost-1.53.0

########################################################################
# PATH
########################################################################

setenv PATH .:${PREFIX}/build/bin
setenv PATH ${PATH}:${PREFIX}/scripts

setenv PATH ${PATH}:${GCC}/bin
setenv PATH ${PATH}:${ROOT}/root/bin
setenv PATH ${PATH}:${PYTHON}/bin
setenv PATH ${PATH}:${SCONS}/bin

### standard system paths
setenv PATH ${PATH}:/site/bin:/apps/bin
setenv PATH ${PATH}:/usr/bin:/bin:/usr/sbin:/sbin

setenv PATH ${PATH}:./bin:./build/bin

########################################################################
# LD_LIBRARY_PATH
########################################################################

### run-time library loading path
setenv LD_LIBRARY_PATH .:${PREFIX}/build/lib

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${GCC}/lib64
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ROOT}/root/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${SCONS}/lib

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${BOOST}/lib
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CERN}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/lib64
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/lib

########################################################################
# PYTHONPATH
########################################################################

### python modules search path

setenv PYTHONPATH ${PREFIX}

########################################################################
# sources for build directories
########################################################################

setenv MYSQLINC /usr/include/mysql
setenv MYSQLLIB /usr/lib64/mysql

setenv BOOSTINC ${BOOST}
setenv BOOSTLIB ${BOOST}/lib

setenv CERNLIB ${CERN}/lib

setenv CLAS6INC ${PREFIX}/build/include
setenv CLAS6LIB ${PREFIX}/build/lib
setenv CLAS6BIN ${PREFIX}/build/bin

########################################################################
# misc
########################################################################

setenv CLAS_PARMS /group/clas/parms

rehash

