#!/bin/bash
#echo "usage : ./test.sh ./[programme] [arg]"
rm callg*
valgrind --tool=callgrind ./hpc_NO.exe
kcachegrind &
