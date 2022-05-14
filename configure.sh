# set "boost_root" and "yaml_root" to the your path of BOOST and YAML
boost_root=$HOME/BOOST
yaml_root=$HOME/YAML


# Start build
boost_include_dir=${boost_root}/include
boost_lib_dir=${boost_root}/lib
yaml_include_dir=${yaml_root}/include
yaml_lib_dir=${yaml_root}/lib

# Cleanup old cache before we configure
# Note:  This does not remove files produced by make.  Use "make clean" for this.
find . -name "CMakeFiles" -exec rm -rf {} \;
rm -f CMakeCache.txt


cmake \
    -DBOOST_INCLUDEDIR=${boost_include_dir} \
    -DBOOST_LIBRARY_DIR=${boost_lib_dir} \
    -DYAML_INCLUDE_DIR=${yaml_include_dir} \
    -DYAML_LIBRARY_DIR=${yaml_lib_dir} \
    ..