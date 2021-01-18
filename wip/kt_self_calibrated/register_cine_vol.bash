#!/usr/bin/env bash


# REGISTER CINE VOLUME

# e.g., bash register_cine_vol.bash ~/path/to/directory/containing/cine_vol/files/


# Input 

REGDIR=$1


# Check that Registration Directory Exists

if [ ! -d "$REGDIR" ]; then
  echo directory $REGDIR does not exist
  exit 1
else


	# Manage Paths to Allow Queueing of Jobs using TaskSpooler

	ORIG_PATH=$(pwd)
	SCRIPT_PATH=$(dirname `which $0`)

	echo REGISTERING CINE VOLUMES
	echo $REGDIR


	# Variables

	CINEVOLORIG="cine_vol_orig.nii.gz"
	CINEVOLHHH="cine_vol_hhh.nii.gz"
	VELVOL0HHH="vel_vol_0_hhh.nii.gz"
	VELVOL1HHH="vel_vol_1_hhh.nii.gz"
	VELVOL2HHH="vel_vol_2_hhh.nii.gz"

	# Convert 4D Volumes to 3D Volumes

	mirtk extract-image-volume -n 0 $CINEVOLORIG $CINEVOLORIG
	mirtk extract-image-volume -n 0 $CINEVOLHHH $CINEVOLHHH
	mirtk extract-image-volume -n 0 $VELVOL0HHH $VELVOL0HHH	
	mirtk extract-image-volume -n 0 $VELVOL1HHH $VELVOL1HHH	
	mirtk extract-image-volume -n 0 $VELVOL2HHH $VELVOL2HHH

	VOLS3DORIG="cine_vol_orig_*.nii.gz"
	VOLS3DHHH="cine_vol_hhh_*.nii.gz"
	VOLS3DHHHVEL0="vel_vol_0_hhh_*.nii.gz"
	VOLS3DHHHVEL1="vel_vol_1_hhh_*.nii.gz"
	VOLS3DHHHVEL2="vel_vol_2_hhh_*.nii.gz"

	VOLS3DORIGFILES=($VOLS3DORIG)
	VOLS3DHHHFILES=($VOLS3DHHH)
	VOLS3DHHHVEL0FILES=($VOLS3DHHHVEL0)
	VOLS3DHHHVEL1FILES=($VOLS3DHHHVEL1)
	VOLS3DHHHVEL2FILES=($VOLS3DHHHVEL2)



	### Perform Registration
	for ((i=0;i<${#VOLS3DORIGFILES[@]};++i)); do

		# 0x File Numbering
		if [ $i -lt 10 ]; then
			ii=0$i
		else
			ii=$i
		fi

		CMD="mirtk register ${VOLS3DORIGFILES[i]} ${VOLS3DHHHFILES[i]} -dofout cine_vol_reg_$ii.dof -output cine_vol_hhh_reg_$ii.nii.gz -bg -1 -interp BSpline -model Rigid; "
		echo $CMD >> register.bash
		eval $CMD
		
	done


	# Combine into 4D Volume

	CMD="mirtk combine-images cine_vol_hhh_reg_*.nii.gz -output cine_vol_hhhReg.nii.gz; " # temp name
	echo $CMD >> register.bash
	eval $CMD
	
	
	
	### Perform velocity volume transformation
	
	echo
	echo TRANSFORMING VELOCITY VOLUMES
	echo
	
	for ((i=0;i<${#VOLS3DORIGFILES[@]};++i)); do

		# 0x File Numbering
		if [ $i -lt 10 ]; then
			ii=0$i
		else
			ii=$i
		fi

		# vel component 0
		CMD="mirtk transform-image ${VOLS3DHHHVEL0FILES[i]} vel_vol_0_hhh_trans_$ii.nii.gz -target cine_vol_orig_00.nii.gz -interp BSpline -dofin cine_vol_reg_$ii.dof; "
		echo $CMD >> register.bash
		eval $CMD
					
		# vel component 1
		CMD="mirtk transform-image ${VOLS3DHHHVEL1FILES[i]} vel_vol_1_hhh_trans_$ii.nii.gz -target cine_vol_orig_00.nii.gz -interp BSpline -dofin cine_vol_reg_$ii.dof; "
		echo $CMD >> register.bash
		eval $CMD
				
		# vel component 2
		CMD="mirtk transform-image ${VOLS3DHHHVEL2FILES[i]} vel_vol_2_hhh_trans_$ii.nii.gz -target cine_vol_orig_00.nii.gz -interp BSpline -dofin cine_vol_reg_$ii.dof; "
		echo $CMD >> register.bash
		eval $CMD
		
		
	done


	# Combine into 4D Volume

	CMD="mirtk combine-images vel_vol_0_hhh_trans_*.nii.gz -output vel_vol_0_hhhTrans.nii.gz; " # temp name
	echo $CMD >> register.bash
	eval $CMD
		
	CMD="mirtk combine-images vel_vol_1_hhh_trans_*.nii.gz -output vel_vol_1_hhhTrans.nii.gz; " # temp name
	echo $CMD >> register.bash
	eval $CMD
		
	CMD="mirtk combine-images vel_vol_2_hhh_trans_*.nii.gz -output vel_vol_2_hhhTrans.nii.gz; " # temp name
	echo $CMD >> register.bash
	eval $CMD
		
		
	# Clean Up
	
	CMD="rm cine_vol_hhh_reg_*.nii.gz; rm cine_vol_orig_*.nii.gz; rm cine_vol_hhh_*.nii.gz; "
	echo $CMD >> register.bash
	eval $CMD		

	CMD="mv cine_vol_hhhReg.nii.gz cine_vol_hhh_reg.nii.gz; "
	echo $CMD >> register.bash
	eval $CMD		
	
	CMD="rm vel_vol_*_hhh_trans_*.nii.gz; rm vel_vol_*_hhh_*.nii.gz; "
	echo $CMD >> register.bash
	eval $CMD		

	CMD="mv vel_vol_0_hhhTrans.nii.gz vel_vol_0_hhh_reg.nii.gz; mv vel_vol_1_hhhTrans.nii.gz vel_vol_1_hhh_reg.nii.gz; mv vel_vol_2_hhhTrans.nii.gz vel_vol_2_hhh_reg.nii.gz; "
	echo $CMD >> register.bash
	eval $CMD	
	
	


	# Finish

	echo "cine_vol / vel_vol registration complete"

	cd $ORIG_PATH





fi

