#!/bin/csh

rm -rf build5
mkdir build5

if ( "$TORCH_SWITCH" == "ON" ) then
    echo "Building with Torch integration..."
    cd build5
    echo "TORCH_DIR = ${TORCH_INSTALLED_DIR}"
    setenv TORCH_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0"
    cmake -DCMAKE_PREFIX_PATH=${TORCH_INSTALLED_DIR} -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ..
    make VERBOSE=1
    make install
else
    echo "Building without Torch integration..."
    cmake -H. -Bbuild5
    cd build5
    make VERBOSE=1
    make install
endif

echo "Compilation complete. You can now run test files as '../build5/efpmd/efpmd input.in' from the tests directory."

