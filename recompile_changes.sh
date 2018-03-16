#! /bin/bash
clear
echo "Recompiling to integrate new changes..."
rm -r build
mkdir build
cd build/
echo "Changed to 'build' dir to start compilation..."
cmake ../src -DCMAKE_BUILD_TYPE=RELEASE
make -j4
echo "Compilation done, copying necessary executables to 'build'..."
cp ../ext_bin/sift_cli .
cp ../ext_bin/match_cli .
cp ../ext_bin/deepmatching .
echo "Everything complete, changing dir to 'scripts_python'" 
cd ../scripts_python
echo "All ready to execute next test!"
 
