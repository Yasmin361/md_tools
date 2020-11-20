import os, glob


def main():
    cwd = os.getcwd()
    limits = (1, 5559, 5684, 10807)     # still same, since PIO is reported later..
    filter_inter_intra(cwd, limits)


def sum_by_element(lista, listb):
    sum_list = []
    if lista is None and listb is not None:
        lista = [0] * len(listb)
    if listb is None and lista is not None:
        listb = [0] * len(lista)
    for (item1, item2) in zip(lista,listb):
        sum_list.append(max(item1,item2))   # contacts are either 0 or 1, cannot be bigger than 1
    return sum_list


def binarize(inplist):
    inplist[:] = ["1" if int(x) > 0 else x for x in inplist]


def filter_inter_intra(cwd, limits):
    arrmin = limits[0]
    arrmax = limits[1]
    ntrmin = limits[2]
    ntrmax = limits[3]
    for rep in glob.glob(os.path.join(cwd, "*.tml"), recursive=True):
        interarr = dict()   # intermolecular contacts arrestin
        interntr = dict()   # intermolecular contacts receptor
        with open(os.path.join(rep), "r") as r:
            line = r.readline()
            while line:
                while "#" in line:
                    line = r.readline()
                pair = "{r}-{n}".format(r=line.split()[1], n=line.split()[2].split("==")[0])  # name of residue
                line = r.readline()
                parta = int(line.split()[5].rstrip(")"))  # atomindex of first residue
                pair_description = parta
                ctcs = list()
                line = r.readline()
                while line and "freeSelLabel" not in line:
                    a = line.split()[0]
                    b = line.split()[1]
                    ctcs.insert(int(a), b)
                    line = r.readline()
                binarize(ctcs)
                if (arrmin <= parta) & (parta <= arrmax):   # inter arr
                    if pair not in interarr:
                        interarr[pair] = [pair_description, ctcs]
                    else:
                        interarr[pair] = [pair_description, sum_by_element(interarr[pair][1], ctcs)]
                elif (ntrmin <= parta) & (parta <= ntrmax):   # inter ntr
                    if pair not in interntr:
                        interntr[pair] = [pair_description, ctcs]
                    else:
                        interntr[pair] = [pair_description, sum_by_element(interntr[pair][1], ctcs)]
                else:
                    pass

        with open(os.path.join(cwd, "{rep}_PIO_contacts_inter_arr.list".format(rep=rep)), "w") as fo:
            for k,v in interarr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_PIO_contacts_inter_ntr.list".format(rep=rep)), "w") as fo:
            for k,v in interntr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))


if __name__ == "__main__":
    main()
