#!/usr/bin/env bash


# RECONSTRUCT CINE VOL USING EXISTING TRANSFORMATION PARAMETERS

### DIDN'T USE THIS IN THE END AS LUCILIO RECON HAD ENOUGH SIGNAL, HENCE RECONSTRUCTED WITHOUT USING INFORMATION FROM STANDARD RECONSTRUCTION



# e.g., bash recon_cine_vol_harmonic.bash ~/path/to/top/level/recon/directory/ ~/path/to/top/level/first_recon/directory/ cine_vol_harmonic


# Input 

RECONDIR=$1
RECONTGTDIR=$2
VOLDESC=$3


# Check that Recon Directory Exists

if [ ! -d "$RECONDIR" ]; then
  echo directory $RECONDIR does not exist
  exit 1
else


# Manage Paths to Allow Queueing of Jobs using TaskSpooler

ORIG_PATH=$(pwd)
SCRIPT_PATH=$(dirname `which $0`)

RECONVOLDIR=$RECONDIR/$VOLDESC
mkdir -p $RECONVOLDIR
cd $RECONVOLDIR

echo RECON CINE VOLUME
echo $RECONVOLDIR


# Variables 

RECON=$VOLDESC.nii.gz
STACKS="../data/s*_rlt_ab_hrm.nii.gz"
MASKCINEVOL="$RECONTGTDIR/cine_vol/mask_cine_vol.nii.gz"
THICKNESSFILE="$RECONTGTDIR/data/slice_thickness.txt"
#EXCLUDESTACKFILE="../data/force_exclude_stack.txt"
EXCLUDESLICEFILE="../data/force_exclude_slice.txt"
EXCLUDEFRAMEFILE="../data/force_exclude_cine_vol.txt"
RESOLUTION=1.25
NMC=0  # 0 = no registration
NSR=10 # num super res iterations
NUMCARDPHASE=25
STACKDOFDIR="$RECONTGTDIR/dc_vol/stack_transformations"
SLICEDOFDIR="$RECONTGTDIR/dc_vol/slice_transformations"
CINETRANSDIR="$RECONTGTDIR/cine_vol/transformations"
#NUMSLICEENLARGE=20  # num. slices
#ROIBLURCINEVOL=10  # mm
#ROIBLURTHRESH=0.002699796063 # 3*sigma
#ROIBLURCINEVOLSIGMA=$(echo "$ROIBLURCINEVOL/3" | bc -l)
#NUMROICLOSINGITER=9
MEANRRFILE="$RECONTGTDIR/cardsync/mean_rrinterval.txt"
RRINTERVALSFILE="$RECONTGTDIR/cardsync/rrintervals.txt"
CARDPHASESFILE="$RECONTGTDIR/cardsync/cardphases_interslice_cardsync.txt"


# Setup

ITER=$(($NMC+1))
THICKNESS=$(cat $THICKNESSFILE)
MEANRR=$(cat $MEANRRFILE)
RRINTERVALS=$(cat $RRINTERVALSFILE)
CARDPHASES=$(cat $CARDPHASESFILE)
NUMSTACK=$(ls -1 ../data/s*_dc_ab_hrm.nii.gz | wc -l);
NUMSLICE=$(eval "wc -w $RRINTERVALSFILE | awk -F' ' '{print \$1}'" )
NUMFRAME=$(eval "wc -w $CARDPHASESFILE | awk -F' ' '{print \$1}'" )
STACKFILES=($STACKS)
#STACKMASKALL="../mask/s*_mask_heart.nii.gz"
#STACKMASKFILES=(../mask/s*_mask_heart.nii.gz)
#EXCLUDESTACK=$(cat $EXCLUDESTACKFILE)
#NUMEXCLUDESTACK=$(eval "wc -w $EXCLUDESTACKFILE | awk -F' ' '{print \$1}'" )
EXCLUDESLICE=$(cat $EXCLUDESLICEFILE)
NUMEXCLUDESLICE=$(eval "wc -w $EXCLUDESLICEFILE | awk -F' ' '{print \$1}'" )
EXCLUDEFRAME=$(cat $EXCLUDEFRAMEFILE)
NUMEXCLUDEFRAME=$(eval "wc -w $EXCLUDEFRAMEFILE | awk -F' ' '{print \$1}'" )


# # Identify Number of Slices in Each Stack

# declare -a ARRAYNUMSLICEINSTACK
# STACKINDEX=0
# for STACK in ${STACKFILES[@]}
# do
  # CMD="mirtk info $STACK | grep \"Image dimensions\" | awk -F' ' '{print \$6}'"
  # NUMSLICEINSTACK=$(eval $CMD)
  # ARRAYNUMSLICEINSTACK[$STACKINDEX]=$NUMSLICEINSTACK
  # ((STACKINDEX++))
# done


# Recon Cine Volume

echo reconstructing adaptive cine volume: $RECON
CMD="mirtk reconstructCardiac $RECON $NUMSTACK $STACKS -thickness $THICKNESS -dofin $STACKDOFDIR/stack-transformation*.dof -transformations $CINETRANSDIR -mask $MASKCINEVOL -iterations $ITER -rec_iterations $NSR -resolution $RESOLUTION -force_exclude_stack 0 -force_exclude_sliceloc $NUMEXCLUDESLICE $EXCLUDESLICE -force_exclude $NUMEXCLUDEFRAME $EXCLUDEFRAME -numcardphase $NUMCARDPHASE -rrinterval $MEANRR -rrintervals $NUMSLICE $RRINTERVALS -cardphase $NUMFRAME $CARDPHASES -debug > log-main.txt"

echo $CMD > recon.bash
eval $CMD


# > log-main.txt

# Clean Up

echo "" >> recon.bash
CMD="mkdir sr_iterations; mv averagebias* sr_iterations/; mv bias* sr_iterations/; mv correctedstack* sr_iterations/; mv errorstack* sr_iterations/; mv simstack* sr_iterations/; mv simweightstack* sr_iterations/; mv super_* sr_iterations/; mv weights* sr_iterations/;"
echo $CMD >> recon.bash
eval $CMD
echo "" >> recon.bash
CMD="mkdir transformations; mv transformation*.dof transformations/;"
echo $CMD >> recon.bash
eval $CMD


# Finish

echo "volume reconstruction complete"

cd $ORIG_PATH





fi

