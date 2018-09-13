#!/usr/bin/env bash

# R -e "shiny::runApp('~/shinyapp')"

# what goes in there? the project name?
# assuming the current working directory is proper for now
# ....could use $0 to ge the location of this script

DIRPATH=$(dirname $0)/../
echo "PATH? $DIRPATH"

R -e "shiny::runApp('$DIRPATH/peax-pmp')"
