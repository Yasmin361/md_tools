import os, glob


def main():
    cwd = os.getcwd()
    limits = (1, 5559, 5684, 10807)
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


def correct_pair_order(dictionary):
    for key in dictionary:
        indices = dictionary[key][0].split()    # list with two elements
        if int(indices[0]) > int(indices[1]):   # wrong order
            residues = key.split("--")
            newkey = "{}--{}".format(residues[1], residues[0])
            new_description = "{}\t{}".format(indices[1], indices[0])
            dictionary[key][0] = new_description
            dictionary[newkey] = dictionary.pop(key)
        else:
            pass
    return dictionary


def filter_inter_intra(cwd, limits):
    arrmin = limits[0]
    arrmax = limits[1]
    ntrmin = limits[2]
    ntrmax = limits[3]

    for rep in glob.glob(os.path.join(cwd, "*.tml"), recursive=True):
        inter = dict()      # intermolecular contacts
        intraarr = dict()   # intramolecular contacts arrestin
        intrantr = dict()   # intramolecular contacts receptor
        ligand = dict()     # contacts that fit in neither of above, because they contain ligand residues
        with open(os.path.join(rep), "r") as r:
            line = r.readline()
            while line:
                while "#" in line:
                    line = r.readline()
                pair = line.split()[1]  # name of residues
                line = r.readline()
                parta = int(line.split()[5]) # atomindex of first residue
                partb = int(line.split()[6].rstrip(")")) # atomindex of second residue
                pair_description = "{}\t{}".format(parta, partb)
                ctcs = list()
                line = r.readline()
                while line and "freeSelLabel" not in line:
                    a = line.split()[0]
                    b = line.split()[1]
                    ctcs.insert(int(a), b)
                    line = r.readline()
                if (arrmin <= parta) & (parta <= arrmax) & (arrmin <= partb) & (partb <= arrmax):   # intra arr
                    # intraarr[pair] = [pair_description, ctcs]
                    if pair not in intraarr:
                        intraarr[pair] = [pair_description, ctcs]
                    else:
                        intraarr[pair] = [pair_description, sum_by_element(intraarr[pair][1], ctcs)]
                elif (ntrmin <= parta) & (parta <= ntrmax) & (ntrmin <= partb) & (partb <= ntrmax):   # intra ntr
                    # intrantr[pair] = [pair_description, ctcs]
                    if pair not in intrantr:
                        intrantr[pair] = [pair_description, ctcs]
                    else:
                        intrantr[pair] = [pair_description, sum_by_element(intrantr[pair][1], ctcs)]
                elif ((arrmin <= parta) & (parta <= arrmax) & (ntrmin <= partb) & (partb <= ntrmax)) | \
                        ((ntrmin <= parta) & (parta <= ntrmax) & (arrmin <= partb) & (partb <= arrmax)):  # inter
                    # inter[pair] = [pair_description, ctcs]
                    if pair not in inter:
                        inter[pair] = [pair_description, ctcs]
                    else:
                        inter[pair] = [pair_description, sum_by_element(inter[pair][1], ctcs)]
                else:
                    # ligand[pair] = [pair_description, ctcs]
                    if pair not in ligand:
                        ligand[pair] = [pair_description, ctcs]
                    else:
                        ligand[pair] = [pair_description, sum_by_element(ligand[pair][1], ctcs)]
        correct_pair_order(inter)
        correct_pair_order(intraarr)
        correct_pair_order(intrantr)
        correct_pair_order(ligand)

        with open(os.path.join(cwd, "{rep}_saltbridges_intra_arr.list".format(rep=rep)), "w") as fo:
            for k,v in intraarr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_saltbridges_intra_ntr.list".format(rep=rep)), "w") as fo:
            for k,v in intrantr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_saltbridges_inter.list".format(rep=rep)), "w") as fo:
            for k,v in inter.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_saltbridges_ligand.list".format(rep=rep)), "w") as fo:
            for k,v in ligand.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))


if __name__ == "__main__":
    main()
