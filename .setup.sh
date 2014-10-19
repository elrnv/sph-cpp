#!/bin/bash

SPH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PATH=$PATH:$SPH/bin

alias sph='r $SPH/bin/sim.app'
