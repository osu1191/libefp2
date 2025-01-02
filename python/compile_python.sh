rm -rf build
mkdir build
cd build

cmake -DCMAKE_INSTALL_PREFIX=$LIBEFP_DIR/installed -DCMAKE_PREFIX_PATH="$CONDA_PREFIX;$LIBEFP_DIR/installed;$TORCH_INSTALLED_DIR" ..

make VERBOSE=1 
make install
