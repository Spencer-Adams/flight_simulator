#!/bin/bash

printf "____________________________________________________________\n"
printf "____________________________________________________________\n"
printf "Building PID Controller Example ...\n"
SECONDS=0
g++ -std=c++11 \
../../src/LUDv.cpp \
../../src/interface/udp_port.cpp \
../../src/interface/serial_port.cpp \
../../src/interface/usb_port.cpp \
../../src/interface/connection.cpp \
../../src/interface/interface.cpp \
-I/usr/local/include/libusb-1.0 \
-L/usr/local/lib \
-lusb-1.0 \
main.cpp \
-o PID.out

duration=$SECONDS
printf "Time to build (min:sec) = $((duration / 60)):$((duration % 60))\n"
printf "Build Complete\n"
printf "____________________________________________________________\n"
printf "____________________________________________________________\n"

