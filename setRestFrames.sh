#!/bin/bash

if [ -d "${ROOTCOREBIN}" ]; then

    PATH=${ROOTCOREBIN}/../RestFrames/bin:$PATH; export PATH
    
    LD_LIBRARY_PATH=${ROOTCOREBIN}/../RestFrames/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
    
    DYLD_LIBRARY_PATH=${ROOTCOREBIN}/../RestFrames/lib:$DYLD_LIBRARY_PATH; export DYLD_LIBRARY_PATH
    
    SHLIB_PATH=${ROOTCOREBIN}/../RestFrames/lib:$SHLIB_PATH; export SHLIB_PATH
    
    LIBPATH=${ROOTCOREBIN}/../RestFrames/lib:$LIBPATH; export LIBPATH
else
    echo "ROOTCOREBIN is not defined. Please setup RootCore before calling this script."
fi
