#!/usr/bin/python

import os
import argparse
import subprocess


def parse_arguments():
    parser = argparse.ArgumentParser(description='Fetal CMR 4D flow render')
    parser.add_argument('-path', type=dir_path, metavar='reconDir', help='Path to reconDir (e.g.: /path/to/fcmr214)')
    parser.add_argument('-fcmrnum', type=int, metavar='fcmrnum', help='3-digit fcmr number (e.g.: 214)')

    return parser.parse_args()


def dir_path(path):
    if os.path.isdir(path):
        print('Directory found ... ')
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


def main():
    args = parse_arguments()

    # TODO: make pvpypath an optional function argument
    # TODO: make pvscriptpath relative to fcmr_make4dflow.py
    pvpypath = 'C:\\Program Files\\ParaView 5.4.1-Qt5-OpenGL2-Windows-64bit\\bin\\pvpython.exe'
    pvscriptpath = 'C:\\Users\\tr17\\Documents\\CODE_Projects\\fetal_cmr_4dflow\\paraview\\fcmr_4dflow_callpv.py' \
        + ' -path ' + args.path \
        + ' -fcmrnum ' + str(args.fcmrnum)

    print([pvpypath, pvscriptpath])

    subprocess.call([pvpypath, pvscriptpath])


if __name__ == "__main__":
    main()


# # positional arguments
# parser.add_argument('integers', metavar='reconDir', type=int, nargs='+',
#                     help='Path to reconDir (e.g.: c_fcmr214)')
#
# # optional arguments
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')



# args = parser.parse_args()
#
# # path/variable admin
# fcmrNum = 202
# foldExt = ''
# pvStateExt = ''
# velDir = '\\vel_vol_4d_4stk'
# # velDir = '\\vel_vol_trans_4d'
#
# # pvFold = '\\paraview' + '\\'
# pvFold = '\\paraview_polyCorr' + '\\'
# # pvFold = '\\paraview_polyCorr_aorta_v2_LV_RV_LOT_ROT_LA_RA_IVC_SVC_PA_DA_v2' + '\\'
#
# # fcmrDir = 'I:\\fcmr_4d_clinical_recons\\c_fcmr'
# # fcmrDir = 'I:\\fcmr_4d_chloe_recons\\c_fcmr'
# # fcmrDir = 'E:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\Chloe_Recons\\fcmr'
# fcmrDir = 'E:\\Users\\tr17\\Documents\\Projects\\PC_Fetal_CMR\\Data\\4D_Flow_Paper\\fcmr'
# path = fcmrDir + str(fcmrNum) + '' + foldExt + velDir + pvFold
#
# pvStateFilename = 'fcmr' + str(fcmrNum) + '_' + pvStateExt + 'paraview.pvsm'
# # pvStateFilename = 'c_fcmr' + str(fcmrNum) + '_' + pvStateExt + 'paraview.pvsm'
#
#
# filenameCine = 'cine_vol_masked_t'
# filenameVelVol = 'vel_vol_masked_VxVy-Vz_t'
# numFrames = 25
#
# print(args.accumulate(args.integers))
#
# print(path)
