import subprocess
import os
import shlex
import sys

def main():
    mode = sys.argv[1]
    if mode == "cluster":
        cwd = os.getcwd()
    elif mode == "local":
        cwd = "/home/yasmin/Desktop/run002/"   # local
    elif mode == "cluster" or "local":
        print("\nPlease supply mode in call - either *local* or *cluster*")
        sys.exit(1)
    initialize(mode)    # either "local" or "cluster"
    calculatecov(cwd, mode)


# load gromacs on local or cluster
def initialize(machine):
    if machine == "local":
        src = shlex.split("env -i bash -c 'source /usr/local/gromacs/bin/GMXRC'")
        process = subprocess.Popen(src, stdout=subprocess.PIPE)
        for line in process.stdout:
            print(line)
    elif machine == "cluster":
        # subprocess.call("module load cuda/10.2.89", shell=True)
        # subprocess.call("module load GROMACS/2020.2_GPU", shell=True)
        print("\nWriting initializaton for cluster to startPCA.sh.")
    else:
        print("\nInitialization didn't work.")

# small utility for subsetting list of files in directoy to trajectories only
def onlyxtc(elements):
    return [element for element in elements if element[-4:] == ".xtc"]


# paths to in and out files
def getpaths(cwd, mode):
    if mode == "local":
        trjs_path = os.path.join(cwd, "analyses_2_28", "trjs")          # trajectories local
        ndx_path = os.path.join(cwd, "analyses_2_28", "elements_ndx")   # index files local
    elif mode == "cluster":
        trjs_path = os.path.join(cwd, "results")                        # trajectories cluster
        ndx_path = os.path.join(cwd, "elements_ndx")                    # index files cluster
    trajectories = onlyxtc(os.listdir(trjs_path))
    subgroups = ["cAlpha", "arrcAlpha", "ntrcAlpha"]                    # supply appropriate index file names
    return trajectories, trjs_path, subgroups, ndx_path


# calculate covariance matrix and do eigenvector decomposition
def calculatecov(cwd, mode):
    trajectories, trjs_path, subgroups, ndx_path = getpaths(cwd, mode)
    if mode == "local":
        for traj in trajectories:
            trajn = traj.split(".")
            trj_path = os.path.join(trjs_path, traj)
            for element in subgroups:
                ndx_name = element + ".ndx"
                ndxel_path = os.path.join(ndx_path, ndx_name)
                outpath = os.path.join(cwd, "analyses_2_28", "pca", "gmx_covar", trajn[0], element)  # outpath local
                subprocess.call("mkdir -p " + outpath, shell=True)
                pcacommand = "gmx covar -f {trj} -s {c}/step7_production.tpr -n {ndx} -o {out}/{ele}_ev_spectrum.xvg -v {out}/{ele}_ev.trr -av {out}/{ele}_average.pdb -xpm {out}/{ele}_covar_matrix.xpm".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, ele=element)
                subprocess.call(pcacommand, shell=True)
                evdcommand = "gmx anaeig -v {out}/{ele}_ev.trr -f {trj} -eig {out}/{ele}_ev_spectrum.xvg -n {ndx} -s {c}/step7_production.tpr -first 1 -last 1 -nframes 100 -extr {out}/{ele}_eigv1.pdb".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, ele=element)
                subprocess.call(evdcommand, shell=True)
                ev12command = "gmx anaeig -v {out}/{ele}_ev.trr -f {trj} -eig {out}/{ele}_ev_spectrum.xvg -n {ndx} -s {c}/step7_production.tpr -first 1 -last 2 -2d {out}/{ele}_proj1v2.xvg".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, ele=element)
                subprocess.call(ev12command, shell=True)
                print("\nYay, covariance matrix and eigenvector decomposition done!")
    elif mode == "cluster":
        bashscript = []
        for traj in trajectories:
            trajn = traj.split(".")
            trj_path = os.path.join(trjs_path, traj)
            for element in subgroups:
                ndx_name = element + ".ndx"
                ndxel_path = os.path.join(ndx_path, ndx_name)
                outpath = os.path.join(cwd, "results", "pca", trajn[0], element)  # outpath cluster
                bashscript.append("mkdir -p " + outpath)
                bashscript.append("gmx covar -f {trj} -s {c}/step7_production.tpr -n {ndx} -o {out}/{ele}_ev_spectrum.xvg -v {out}/{ele}_ev.trr -av {out}/{ele}_average.pdb -xpm {out}/{ele}_covar_matrix.xpm".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, ele=element))
                bashscript.append("gmx anaeig -v {out}/{ele}_ev.trr -f {trj} -eig {out}/{ele}_ev_spectrum.xvg -n {ndx} -s {c}/step7_production.tpr -first 1 -last 1 -nframes 100 -extr {out}/{ele}_eigv1.pdb".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, ele=element))
                bashscript.append("gmx anaeig -v {out}/{ele}_ev.trr -f {trj} -eig {out}/{ele}_ev_spectrum.xvg -n {ndx} -s {c}/step7_production.tpr -first 1 -last 2 -2d {out}/{ele}_proj1v2.xvg".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, ele=element))
        with open("startPCA.sh", "w") as f:
            f.write("#!/bin/bash\n#SBATCH --gres=gpu:1\n#SBATCH --ntasks-per-node=9\n#SBATCH --cpus-per-task=2\n#SBATCH --mem=4096\n#SBATCH --mem-bind=local\n#SBATCH --nice=0\n#SBATCH --time=12:00:00\n#SBATCH --job-name=NTS1R_PCA\nmodule load cuda/10.2.89\nmodule load GROMACS/2020.2_GPU\n")
            f.write("\n".join(bashscript))
        print("\nYay, I wrote the startPCA.sh, let's go!")

if __name__ == "__main__":
    main()
    
    
