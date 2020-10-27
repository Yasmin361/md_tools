import subprocess
import os
import shlex
import sys


def main():
    mode = sys.argv[1]
    if mode == "cluster":
        cwd = os.getcwd()
    elif mode == "local":
        cwd = "/home/yasmin/Desktop/run002/"  # local
    elif mode == "cluster" or "local":
        print("\nPlease supply mode in call - either *local* or *cluster*")
        sys.exit(1)
    initialize(mode)  # either "local" or "cluster"
    calculatehelix(cwd, mode)


# small utility for subsetting list of files in directoy to trajectories only
def onlyxtc(elements):
    return [element for element in elements if element[-4:] == ".xtc"]


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
        print("\nWriting initializaton for cluster to starthelix.sh.")
    else:
        print("\nInitialization didn't work.")


# paths to in and out files
def getpaths(cwd, mode):
    if mode == "local":
        trjs_path = os.path.join(cwd, "analyses_2_28", "trjs")  # trajectories local
    elif mode == "cluster":
        trjs_path = os.path.join(cwd, "results")  # trajectories cluster
    trajectories = onlyxtc(os.listdir(trjs_path))
    #subgroups = os.listdir(ndx_path)
    return trajectories, trjs_path


# calculates helix properties of given groups over trajectories
def calculatehelix(cwd, mode):
    trajectories, trjs_path = getpaths(cwd, mode)
    if mode == "local":
        """
        for traj in trajectories:
            trajn = traj.split(".")[0]
            trj_path = os.path.join(trjs_path, traj)
            outpath = os.path.join(cwd, "analyses_2_28", "rmsd", trajn)  # outpath local
            subprocess.call("mkdir -p " + outpath, shell=True)
            for element in subgroups:
                ndx_name = element
                ndxel_path = os.path.join(ndx_path, ndx_name)
                rmscommand = "gmx rms -f {trj} -s {c}/step7_production.tpr -n {ndx} -o {out}/rmsd_{tn}_{ele}.xvg -tu ns".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, tn=trajn.split("_")[2], ele=element)
                subprocess.call(rmscommand, shell=True)
        """
    elif mode == "cluster":
        bashscript = []
        for traj in trajectories:
            trajn = traj.split(".")[0]
            trj_path = os.path.join(trjs_path, traj)
            outpath = os.path.join(cwd, "results", "helicity", trajn)  # outpath cluster
            bashscript.append("mkdir -p " + outpath)
            bashscript.append(
                "gmx rama -f {trj} -s {c}/step7_production.tpr -b 2500000 -e 3000000 "
                "-o {out}/ramachandran_{tn}.xvg".format(
                    trj=trj_path, c=cwd, out=outpath,
                    tn=trajn.split("_")[3]))  #check trajn is actually at position [3]
        with open("starthelicity.sh", "w") as f:
            f.write(
                "#!/bin/bash\n#SBATCH --gres=gpu:1\n#SBATCH --ntasks-per-node=9\n#SBATCH --cpus-per-task=2"
                "\n#SBATCH --mem=4096\n#SBATCH --mem-bind=local\n#SBATCH --nice=0\n#SBATCH --time=12:00:00"
                "\n#SBATCH --job-name=NTS1R_helicity\nmodule load cuda/10.2.89\nmodule load GROMACS/2020.2_GPU\n")
            f.write("\n".join(bashscript))
        print("\nYay, I wrote the starthelix.sh, let's go run it!")

if __name__ == "__main__":
    main()


