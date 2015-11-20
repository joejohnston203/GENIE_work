#!/bin/bash

# Set directory of repository cloned from github
root_download="/home/joe/work/ignore/root_download"
ver="v5-34-34" # Set version to use
build_dir="build_"$ver # Make this directory to store the build
# Install location
install_dir="/usr/local/root"
#directory of pythia6 library
pythia6lib="/opt/pythia/v6_424/lib/libPythia6.so"

# if the build or install directories exist stop
if [ -d $build_dir -o -d $install_dir ]; then
    exit
fi

cwd=$(pwd)

cd $root_download
# Check if a branch for this version exists already
git show-ref --verify --quiet refs/heads/$ver
# $? == 0 means local branch with <branch-name> exists.
if [ $? == 0 ]; then
    # The branch already exists. Print and continue.
    echo 'The branch '$ver' already exists. Continuing.'
else
    echo 'Creating branch '$ver'.'
    git checkout -b $ver $ver
fi

cd $cwd
mkdir $build_dir
cd $build_dir

# Set up build directory using mathmore and pythia6
cmake -Dmathmore:BOOL=ON \
    -DPYTHIA6_LIBRARY:FILEPATH=$pythia6lib -Dpythia6:BOOL=ON \
    $root_download 

# Create build directory. This command takes a while to run.
cmake --build .

# Pause and wait for user input, so the user knows that
# building has been completed. If installation is interrupted,
# only the install and source commands must be run again
read -p "Build directory completed. Press enter to install."

# Set install directory and install
sudo cmake -DCMAKE_INSTALL_PREFIX=$install_dir -P cmake_install.cmake \
    --build . --target install

# The website says to run source to set up the environment, but I
# don't know what this does.
#source $install_dir
