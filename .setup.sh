#!/bin/bash

SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SPH="$( cd $SRC/.. && pwd )"
PATH=$PATH:$SPH/bin

alias sph='$SPH/bin/sim.app/Contents/MacOS/sim'
