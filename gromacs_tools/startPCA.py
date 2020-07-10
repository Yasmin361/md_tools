import subprocess
import os
import shlex


def main():
    #cwd = os.getcwd()
    cwd = "/home/yasmin/Desktop/run002/"

    initialize()
    calculatecov(cwd)

    
# load gromacs on local
def initialize():
    src = shlex.split("env -i bash -c 'source /usr/local/gromacs/bin/GMXRC'")
    process = subprocess.Popen(src, stdout=subprocess.PIPE)
    for line in process.stdout:
        print(line)

    # load gromacs on cluster:
    # subprocess.call("module load GROMACS/2020.2_GPU",shell=True)


# paths to in and out files
def getpaths(cwd):
    trjs_path = os.path.join(cwd, "analyses_2_28", "trjs")            #supply path to trajectories
    ndx_path = os.path.join(cwd, "analyses_2_28", "elements_ndx")     #supply path to index files
    trajectories = os.listdir(trjs_path)
    subgroups = ["cAlpha", "arrcAlpha", "ntrcAlpha"]                  #supply appropriate index file names
    return trajectories, trjs_path, subgroups, ndx_path


# invoke covariance matrix construction using gmx covar
def calculatecov(wd):
    trajectories, trjs_path, subgroups, ndx_path = getpaths(wd)
    for traj in trajectories:
        trajn = traj.split(sep=".")
        trj_path = os.path.join(trjs_path,traj)
        for element in subgroups:
            ndx_name = element + ".ndx"
            ndxel_path = os.path.join(ndx_path,ndx_name)
            outpath = os.path.join(wd, "analyses_2_28", "pca", "gmx_covar", trajn[0], element)
            subprocess.call("mkdir -p " + outpath, shell=True)
            pcacommand = "gmx covar -f {trj} -s {c}step7_production.tpr  -n {ndx} -o {out}.xvg -v {out}_ev.trr -av {out}_average.pdb -xpm {out}_covar_matrix.xpm".format(trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath)
            subprocess.call(pcacommand, shell=True)
            evdcommand = "gmx anaeig -v {out}_ev.trr -f {trj} -eig {out}_ev_spectrum.xvg -n {ndx} -s {c}step7_production.tpr -first 1  -last 1 -nframes 100 -extr {out}_eigv1.pdb".format(trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath)
            subprocess.call(evdcommand, shell=True)
            ev12command = "gmx anaeig -v {out}_ev.trr -f {trj} -eig {out}_ev_spectrum.xvg -n {ndx} -s {c}step7_production.tpr -first 1  -last 1 -nframes 100 -extr {out}_eigv1.pdb".format(trj=trj_path, c=cwd, ndx=ndxel_path, out=outpath)
            subprocess.call(ev12command, shell=True)


if __name__ == "__main__":
    main()
    
    
