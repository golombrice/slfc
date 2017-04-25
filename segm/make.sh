#!/bin/bash

g++ -c ms.cpp
g++ -c msImageProcessor.cpp
g++ -c msSys.cpp
g++ -c msSysPrompt.cpp
g++ -c RAList.cpp
g++ -c rlist.cpp
ar rcs libms.a *.o
