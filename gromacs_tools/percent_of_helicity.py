import os, glob


def main():
    cwd = os.getcwd()
    residue_set = ("TYR-63", "GLY-64", "ARG-65", "GLU-66", "ASP-67", "LEU-68", "ASP-69", "VAL-70", "LEU-71", "GLY-72", "LEU-73", "THR-74", "PHE-75", "ARG-76", "LYS-77", "ASP-78")
    first_res = "TYR-63"
    calculate_percentage(cwd, first_res, residue_set)


def calculate_percentage(cwd, first_res, residue_set):
    frequencies = dict()
    subset = list()
    for res in residue_set:
        frequencies[res] = 0
    for plot in glob.glob(os.path.join(cwd, "**/ramachandran*.xvg"), recursive=True):
        with open(os.path.join(plot), "r") as rama:
            line = rama.readline()
            print(line)
            while line:
                if (len(line.split("#")) == 1) & (len(line.split("@")) == 1):  # skip comment lines
                    if line.split()[2] == first_res:  # starting from first res of sequence
                        for i in range(len(residue_set)):
                            subset.append(line)
                            if (-95 <= float(line.split()[0]) <= -35) & (-70 <= float(line.split()[1]) <= -10): # part of alpha helix def. by Cubellis et al.
                                frequencies[line.split()[2]] += 1
                            line = rama.readline()
                line = rama.readline()

    with open(os.path.join(cwd, "freq_fingerloop_helicity.list"), "w") as fo:
        for item in frequencies:
            fo.write("{k}\t{v}\t{f}\n".format(k=item, v=frequencies[item], f=frequencies[item]/(5*2500)))

    with open(os.path.join(cwd, "fingerloop_ramachandran.list"), "w") as fo:
        for line in subset:
            fo.write(line)


if __name__ == "__main__":
    main()
