#!/bin/tcsh 

# The Basics 
source /group/clas12/packages/setup.csh 

# The JAVA/ROOT 
module load coatjava/6.3.1
#module load clas12/dev
module load root/6.14.04
module load groovy/2.5.6

# Important Env Variables 
setenv BASE_PATH /u/home/dmriser 
setenv JYPATH ${BASE_PATH}/analysis_code
setenv HIPO_TOOLS ${BASE_PATH}/hipo_tools
setenv PATH ${PATH}:${HIPO_TOOLS}/bin
setenv PKG_CONFIG_PATH ${HIPO_TOOLS}/share/pkgconfig
setenv LD_LIBRARY_PATH ${HIPO_TOOLS}/lib
setenv PYTHONPATH ${HIPO_TOOLS}/lib

#  Python3 and cython 
# setenv PATH /apps/anaconda/python3/anaconda201812/bin:${PATH}

# ROOT Newer Build 
# setenv ROOTSYS /work/clas12/tylern/software/root
 setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ROOTSYS}/lib
 setenv PATH ${ROOTSYS}/bin:${PATH}
 setenv PYTHONPATH ${PYTHONPATH}:${ROOTSYS}/lib
