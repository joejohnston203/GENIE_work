sudo cmake -Dmathmore:BOOL=ON \
    -DPYTHIA6_LIBRARY:FILEPATH=/opt/v6_424/lib/libPythia6.so \
    -Dpythia6:BOOL=ON \
    /home/joe/work/root

sudo cmake --build .

sudo cmake -DCMAKE_INSTALL_PREFIX=/opt/root \
    -P cmake_install.cmake \
    --build . --target install
