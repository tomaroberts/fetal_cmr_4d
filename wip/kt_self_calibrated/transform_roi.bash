#!/usr/bin/env bash


## README:
#
# - Old script.
#
# - Transforms vessel ROIs drawn on standard 4dflow volumes into the space of self-calibrated k-t 4dflow volumes
#
# - Abandoned this for ISMRM 2020 abstract - instead just used the vessel ROIs without transforming them
#
# - Retaining script for future reference
#
#################################################################################################################


# e.g., bash transform_roi.bash roiname


# Input 

ROINAME=$1



# Manage Paths to Allow Queueing of Jobs using TaskSpooler

ORIG_PATH=$(pwd)
SCRIPT_PATH=$(dirname `which $0`)

echo "TRANSFORMING ROI ..."
echo $ROINAME


# extract time series
mirtk extract-image-volume -n 0 $ROINAME.nii.gz $ROINAME.nii.gz
ROINAMES="${ROINAME}_*.nii.gz"
ROIFILES=($ROINAMES)

# echo $ROINAMES
# echo ${ROIFILES[9]}
# echo ${#ROIFILES[@]}

### Perform transformation
	
for ((i=0;i<${#ROIFILES[@]};++i)); do

	# 0x File Numbering
	if [ $i -lt 10 ]; then
		ii=0$i
	else
		ii=$i
	fi

	CMD="mirtk transform-image ${ROIFILES[i]} ${ROINAME}_trans_$ii.nii.gz -dofin ../cine_vol_reg_$ii.dof; "
	echo $CMD >> transform.bash
	eval $CMD
	
done


CMD="mirtk combine-images ${ROINAME}_trans_*.nii.gz -output ${ROINAME}hhhTrans.nii.gz; " # temp name
echo $CMD >> transform.bash
eval $CMD


	
# Clean Up

CMD="rm ${ROINAME}_*.nii.gz;"
echo $CMD >> transform.bash
eval $CMD	

CMD="mv ${ROINAME}hhhTrans.nii.gz ${ROINAME}_hhh_reg.nii.gz"
echo $CMD >> transform.bash
eval $CMD

	
	
# Finish

echo "TRANSFORMED ROI."
echo

cd $ORIG_PATH





# fi