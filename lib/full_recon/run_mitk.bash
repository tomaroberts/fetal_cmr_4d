#!/usr/bin/env bash

# WORK AROUND TO RUN MITK WORKBENCH FROM WITHIN MATLAB 2015a
# AS IT ONLY WORKS FROM 2016b AND I DON'T HAVE ACCESS RIGHTS

# THIS IS HORRIFIC I'M SORRY

matlab -nodisplay -nosplash -r "try [~,~]=system('/home/cmo19/MITK-2016.11.0-linux64/MitkWorkbench.sh'); catch; end; quit"







