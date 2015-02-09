#!/bin/bash

SPH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SRC="$( cd $SPH/src && pwd )"
PATH=$PATH:$SPH/bin

alias sph_rel='$SPH/build/release/sph.app/Contents/MacOS/sph'
alias sph_dev='$SPH/build/debug/sph.app/Contents/MacOS/sph'
if [ -e $sph_rel ] ; then
  alias sph='sph_rel'
else
  alias sph='sph_dev'
fi

