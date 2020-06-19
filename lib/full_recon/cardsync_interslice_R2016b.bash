#!/usr/bin/env bash

# WORK AROUND TO RUN cardsync_interslice FROM WITHIN MATLAB 2015a
# AS IT ONLY WORKS FROM 2016b AND I DON'T HAVE ACCESS RIGHTS

# THIS IS HORRIFIC I'M SORRY


matlab -nodisplay -nosplash -r "try addpath( genpath( '/home/cmo19/MATLAB/fetal_cmr_4d-master' ) ); cardsync_interslice_terminal_wrapper(); catch; end; quit"
# matlab -nodisplay -r "addpath( genpath( '/home/cmo19/MATLAB/fetal_cmr_4d-master' ) ); cardsync_interslice_terminal_wrapper();"