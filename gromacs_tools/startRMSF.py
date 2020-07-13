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
    calculatermsf(cwd, mode)


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
        print("\nWriting initializaton for cluster to startRMSF.sh.")
    else:
        print("\nInitialization didn't work.")


# paths to in and out files
def getpaths(cwd, mode):
    if mode == "local":
        trjs_path = os.path.join(cwd, "analyses_2_28", "trjs")  # trajectories local
        ndx_path = os.path.join(cwd, "analyses_2_28", "elements_ndx")  # index files local
    elif mode == "cluster":
        trjs_path = os.path.join(cwd, "results")  # trajectories cluster
        ndx_path = os.path.join(cwd, "elements_ndx")  # index files cluster
    trajectories = onlyxtc(os.listdir(trjs_path))
    subgroups = ["ntr.ndx", "arr.ndx"]  # supply appropriate index file names
    #subgroups = os.listdir(ndx_path)
    return trajectories, trjs_path, subgroups, ndx_path


# calculate rmsd over trajectories
def calculatermsf(cwd, mode):
    trajectories, trjs_path, subgroups, ndx_path = getpaths(cwd, mode)
    if mode == "local":
        for traj in trajectories:
            trajn = traj.split(".")[0]
            trj_path = os.path.join(trjs_path, traj)
            outpath = os.path.join(cwd, "analyses_2_28", "rmsf", trajn)  # outpath local
            subprocess.call("mkdir -p " + outpath, shell=True)
            for element in subgroups:
                ndx_name = element
                ndxel_path = os.path.join(ndx_path, ndx_name)
                rmscommand = "gmx rmsf -f {trj} -s {c}/step7_production.tpr -n {ndx} -res yes -o {out}/rmsf_{tn}_{ele}.xvg -od {out}/rmsf_dev_{tn}_{ele}.xvg -oc {out}/rmsf_correl_{tn}_{ele}.xvg -oq {out}/bfactor_{tn}_{ele}.pdb -ox {out}/xaver_{tn}_{ele}.pdb".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, tn=trajn.split("_")[2], ele=element.split(".")[0])
                subprocess.call(rmscommand, shell=True)
    elif mode == "cluster":
        bashscript = []
        for traj in trajectories:
            trajn = traj.split(".")[0]
            trj_path = os.path.join(trjs_path, traj)
            outpath = os.path.join(cwd, "results", "rmsf", trajn)  # outpath cluster
            bashscript.append("mkdir -p " + outpath)
            for element in subgroups:
                ndx_name = element
                ndxel_path = os.path.join(ndx_path, ndx_name)
                bashscript.append(
                    "gmx rmsf -f {trj} -s {c}/step7_production.tpr -n {ndx} -res yes -o {out}/rmsf_{tn}_{ele}.xvg -od {out}/rmsf_dev_{tn}_{ele}.xvg -oc {out}/rmsf_correl_{tn}_{ele}.xvg -oq {out}/bfactor_{tn}_{ele}.pdb -ox {out}/xaver_{tn}_{ele}.pdb".format(
                    trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath, tn=trajn.split("_")[2], ele=element.split(".")[0]))
        with open("startRMSF.sh", "w") as f:
            f.write(
                "#!/bin/bash\n#SBATCH --gres=gpu:1\n#SBATCH --ntasks-per-node=9\n#SBATCH --cpus-per-task=2\n#SBATCH --mem=4096\n#SBATCH --mem-bind=local\n#SBATCH --nice=0\n#SBATCH --time=12:00:00\n#SBATCH --job-name=NTS1R_RMSF\nmodule load cuda/10.2.89\nmodule load GROMACS/2020.2_GPU\n")
            f.write("\n".join(bashscript))
        print("\nYay, I wrote the startRMSF.sh, let's go run it!")

if __name__ == "__main__":
    main()


