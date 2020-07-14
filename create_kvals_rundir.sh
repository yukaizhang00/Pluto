#!/bin/sh

#~ 05.07.2011 (last modified; TIM)

#------------------------------------------------------------------------------

#set -x;  # May be commented out if desired (if not debugging)

# Distance ourselves (a tiny bit) from any screen litter:
echo;

# Determine fully-qualified filesystem location where this script was
# executed from:
BEG_DIR=`pwd`;

# Determine fully-qualified filesystem location of this script:
#cd `dirname $0`; ANCHOR=`pwd`; cd "$ANCHOR";

#Determine where the script is called.
ANCHOR=`pwd`; cd "$ANCHOR";

# Store this script's name:
SCRIPT_NAME=`basename $0`;

# Function that displays usage information to the user:
#   Argument #1: Specific script name
usage() {
cat <<EOF

-----------------------------------------------------------
Usage: ./$1  <run_dir>  

**NOTE**: Must be executed in the main kdm_kvals_code directory, 
after the code is compiled/built.  The run directory <run_dir> must 
already exist.

run_dir: path to directory that holds the information required to run 
         the code, such as the executables, the ck directory, and an 
         out directory. This script will create some new directories and
         copy some files from the original build directory. 
-----------------------------------------------------------

EOF
exit 1;
}

# Process arguments:

if test "x$1" = x; then
   usage $SCRIPT_NAME;
fi
#if test "x$3" != x; then
#   usage $SCRIPT_NAME;
#fi
run_dir=$1;

# Check for the existence of the given WRF run directory path; if yes,
#  then proceed:

#if test ! -d "$run_dir"; then
#cat <<EOF
#
#$SCRIPT_NAME: ERROR: The given run directory path does not exist; create it before running this script again.
#
#EOF
#   usage $SCRIPT_NAME;
#fi

# Test to see if the run directory already exists. If so, stop the script.
#if test -d "$run_dir"; then
#   cat <<EOF
#   $SCRIPT_NAME: ERROR: That run directory already exists!
#   EOF
#   usage $SCRIPT_NAME;
#fi

cp -r $ANCHOR/*.exe $run_dir/.
cp -r $ANCHOR/ck $run_dir/.
cp -r $ANCHOR/kdm.nml $run_dir/.
if test ! -d "$run_dir/out"; then
   mkdir "$run_dir/out";
fi
echo "$SCRIPT_NAME: Finished."

echo; exit 0;


