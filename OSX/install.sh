#!/bin/sh
 
##
# Install autoconf, automake and libtool smoothly on Mac OS X.
# Newer versions of these libraries are available and may work better on OS X
#
# This script is originally from http://jsdelfino.blogspot.com.au/2012/08/autoconf-and-automake-on-mac-os-x.html
#
 
export build=~/devtools # or wherever you'd like to build
mkdir -p $build
cd $build 
 
##
# Change version numbers here as needed
#

autoconf_version='2.69'
automake_version='1.15'
libtool_version='2.4'

##
# Autoconf
# http://ftpmirror.gnu.org/autoconf

curl -OL http://ftpmirror.gnu.org/autoconf/autoconf-$autoconf_version.tar.gz
tar xzf autoconf-$autoconf_version.tar.gz
cd autoconf-$autoconf_version
./configure --prefix=/usr/local
make
sudo make install
export PATH=/usr/local/bin:$PATH
cd ..
rm -r autoconf-*
 
##
# Automake
# http://ftpmirror.gnu.org/automake
 
curl -OL http://ftpmirror.gnu.org/automake/automake-$automake_version.tar.gz
tar xzf automake-$automake_version.tar.gz
cd automake-$automake_version
./configure --prefix=/usr/local
make
sudo make install
cd ..
rm -r automake-*
 
##
# Libtool
# http://ftpmirror.gnu.org/libtool
 
curl -OL http://ftpmirror.gnu.org/libtool/libtool-$libtool_version.tar.gz
tar xzf libtool-$libtool_version.tar.gz
cd libtool-$libtool_version
./configure --prefix=/usr/local
make
sudo make install
cd ..
rm -r libtool-*
 
echo "Installation complete."

