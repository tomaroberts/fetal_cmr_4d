import os
import argparse

print('Entered test_parse.py ... ')

# options
# fcmrnum ??? --- parse filename instead?
# path
# pvStateFilename --- '*_paraview.pvsm'

# defaults
filenameCine = 'cine_vol_masked_t'
filenameVelVol = 'vel_vol_masked_VxVy-Vz_t'
numFrames = 25


# functions
def parse_arguments():
    parser = argparse.ArgumentParser(description='Fetal CMR 4D flow render')
    parser.add_argument('-path', type=dir_path, metavar='reconDir', help='Path to reconDir (e.g.: /path/to/fcmr214)')
    parser.add_argument('-fcmrnum', type=int, metavar='fcmrnum', help='3-digit fcmr number (e.g.: 214)')
    parser.add_argument('-outfilename', type=str, metavar='.pvsm Output Filename',
                        help='Paraview state file (e.g.: fcmr214_paraview.pvsm)')

    return parser.parse_args()


def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")


# Run
args = parse_arguments()
# print(makeArrayFilenames(args.path, filenameCine, numFrames))