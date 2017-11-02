#!/bin/bash
echo "usage : ./test.sh ./[programme] [arg]"
rm callg*
valgrind --tool=callgrind  $1 $2
kcachegrind &
