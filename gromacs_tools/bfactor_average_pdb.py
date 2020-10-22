import os


def main():
    cwd = os.getcwd()
    averagebfactor(cwd)


def averagebfactor(cwd):
    allbfactors = list()
    for bffile in os.listdir(cwd):
        bfactors = list()
        with open(os.path.join(cwd, bffile), "r") as pdb:
            print("Reading " + bffile)
            line = pdb.readline()
            while line:
                if line.split()[0] == "ATOM":              # ATOM record
                    bfactor = line[60:65]
                    bfactors.append(bfactor)
                    line = pdb.readline()
                    continue
                line = pdb.readline()
        allbfactors.append(bfactors)
    averages = [(float(i) + float(j) + float(k) + float(l) + float(m))/5 for i, j, k, l, m in zip(allbfactors[0], allbfactors[1],allbfactors[2],allbfactors[3],allbfactors[4])]
    output = list()
    with open(os.path.join(cwd, os.listdir(cwd)[1]), "r") as temp:
        index = 0
        line = temp.readline()
        while line:
            if line.split()[0] == "ATOM":
                bf = str(round(averages[index],2)).rjust(6)
                line = line[:59] + bf + line[66:]
                output.append(line)
                index += 1
                line = temp.readline()
                continue
            output.append(line)
            line = temp.readline()

    with open(os.path.join(cwd, ((os.listdir(cwd)[1]).split("_")[0] + "_" + (os.listdir(cwd)[1]).split("_")[2] + "_avg.pdb")), "w") as outf:
        for line in output:
            outf.write(line)


if __name__ == "__main__":
    main()
