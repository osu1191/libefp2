cmake -S . -B build \
     -G Ninja \
    -D CMAKE_INSTALL_PREFIX=/scratch/gilbreth/paulsk/backup/pyMLMM/libefp2/installed \
    -D CMAKE_PREFIX_PATH=${TORCH_INSTALLED_DIR} 
    -D BUILD_SHARED_LIBS=ON \
    -D LIBEFP_ENABLE_OPENMP=ON

cmake --build build --target install
