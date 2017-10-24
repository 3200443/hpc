#!/bin/bash

valgrind --tool=callgrind  $1 $2
kcachegrind &