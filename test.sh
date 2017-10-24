#!/bin/bash
echo "usage : ./test.sh ./[programme] [arg]"
valgrind --tool=callgrind  $1 $2
kcachegrind &