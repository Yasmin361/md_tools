import glob
import os


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
    for (item1, item2) in zip(lista, listb):
        sum_list.append(max(item1, item2))   # contacts are either 0 or 1, cannot be bigger than 1
    return sum_list


def correct_pair_order(dictionary):
    outdict = dict()
    for key in dictionary:
        indices = dictionary[key][0]   # list with two elements
        if indices[0] > indices[1]:   # wrong order
            newkey = (key[1], key[0])
            new_description = (indices[1], indices[0])
            dictionary[key][0] = new_description
            outdict[newkey] = dictionary[key]
        else:
            pass
    return outdict


def onlysubstantial(inpdict,cutoff):
    outdict = dict()
    for key in inpdict:
        sum = 0
        ctclist = inpdict[key][1]
        for entry in ctclist:
            if int(entry) == 1:
                sum += 1
        if sum >= (cutoff/100)*len(ctclist):
            outdict[key] = inpdict[key]
    return outdict


def extract_pair(dictionary, partnera, partnerb):
    outdict = dict()
    for key in dictionary:
        partners = key
        if partnera in partners and partnerb in partners:
            outdict[key] = dictionary[key]
    return outdict


def get_residue_indices(cwd, file):
    print("Starting to read residue names...")
    indexdict = dict()
    with open(os.path.join(cwd, file), "r") as r:
        line = r.readline()
        while line:
            if "ATOM" not in line:
                line = r.readline()
                if not line:
                    break
                continue
            index = int(line[5:11])
            indexdict[index] = "".join((line[17:20], str(int(line[22:26]))))  # looks weird but gets rid of whitespaces
            line = r.readline()
    print("...finished reading residue names!")
    return indexdict


def assign_residue_names(dictin, namedict):
    dictout = dict()
    for key in dictin:
        nametupel = (namedict[key[0]], namedict[key[1]])
        dictout[nametupel] = dictin[key]
    return dictout


def filter_inter_intra(cwd, limits):
    arrmin = limits[0]
    arrmax = limits[1]
    ntrmin = limits[2]
    ntrmax = limits[3]
    resnames = get_residue_indices(cwd, "step6.6_equilibration.part0005_prot_nopbc.pdb")

    for rep in glob.glob(os.path.join(cwd, "*.tml"), recursive=True):
        inter = dict()      # intermolecular contacts
        intraarr = dict()   # intramolecular contacts arrestin
        intrantr = dict()   # intramolecular contacts receptor
        ligand = dict()     # contacts that fit in neither of above, because they contain ligand residues
        with open(os.path.join(rep), "r") as r:
            print("Stating to read from {} ...".format(rep))
            line = r.readline()
            while line:
                while "#" in line:
                    line = r.readline()
                pair = (int(line.split()[1])+1, int(line.split()[2])+1)  # indices of residues
                line = r.readline()  # jump freeSelString
                ctcs = list()
                line = r.readline()
                while line and "freeSelLabel" not in line:
                    a = line.split()[0]
                    b = line.split()[1]
                    ctcs.insert(int(a), b)
                    line = r.readline()
                if (arrmin <= pair[0]) & (pair[0] <= arrmax) & (arrmin <= pair[1]) & (pair[1] <= arrmax):   # intra arr
                    # intraarr[pair] = [pair_description, ctcs]
                    if pair not in intraarr:
                        intraarr[pair] = [pair, ctcs]
                    else:
                        intraarr[pair] = [pair, sum_by_element(intraarr[pair][1], ctcs)]
                elif (ntrmin <= pair[0]) & (pair[0] <= ntrmax) & (ntrmin <= pair[1]) & (pair[1] <= ntrmax):   # intra ntr
                    # intrantr[pair] = [pair_description, ctcs]
                    if pair not in intrantr:
                        intrantr[pair] = [pair, ctcs]
                    else:
                        intrantr[pair] = [pair, sum_by_element(intrantr[pair][1], ctcs)]
                elif ((arrmin <= pair[0]) & (pair[0] <= arrmax) & (ntrmin <= pair[1]) & (pair[1] <= ntrmax)) | \
                        ((ntrmin <= pair[0]) & (pair[0] <= ntrmax) & (arrmin <= pair[1]) & (pair[1] <= arrmax)):  # inter
                    # inter[pair] = [pair_description, ctcs]
                    if pair not in inter:
                        inter[pair] = [pair, ctcs]
                    else:
                        inter[pair] = [pair, sum_by_element(inter[pair][1], ctcs)]
                else:
                    # ligand[pair] = [pair_description, ctcs]
                    if pair not in ligand:
                        ligand[pair] = [pair, ctcs]
                    else:
                        ligand[pair] = [pair, sum_by_element(ligand[pair][1], ctcs)]

        #  change indextupel to resnames in dict keys
        print("Assigning residue names ...")
        inter = assign_residue_names(inter, resnames)
        intraarr = assign_residue_names(intraarr, resnames)
        intrantr = assign_residue_names(intrantr, resnames)
        ligand = assign_residue_names(ligand, resnames)

        #  rearrange resnames and indextupel to ascending order
        print("Bringing order to residue pairs ...")
        inter = correct_pair_order(inter)
        intraarr = correct_pair_order(intraarr)
        intrantr = correct_pair_order(intrantr)
        ligand = correct_pair_order(ligand)

        with open(os.path.join(cwd, "{rep}_hbonds_intra_arr.list".format(rep=rep)), "w") as fo:
            for k,v in intraarr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_hbonds_intra_ntr.list".format(rep=rep)), "w") as fo:
            for k,v in intrantr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_hbonds_inter.list".format(rep=rep)), "w") as fo:
            for k,v in inter.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_hbonds_ligand.list".format(rep=rep)), "w") as fo:
            for k,v in ligand.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))

        cutoff = 5        # cutoff in percent
        inter = onlysubstantial(inter, cutoff)
        intraarr = onlysubstantial(intraarr, cutoff)
        intrantr = onlysubstantial(intraarr, cutoff)
        with open(os.path.join(cwd, "{rep}_hbonds_inter_min_cutoff_{c}percent.list".format(rep=rep,c=5)), "w") as fo:
            for k,v in inter.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_hbonds_intra_arr_min_cutoff_{c}percent.list".format(rep=rep,c=5)), "w") as fo:
            for k,v in intraarr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))
        with open(os.path.join(cwd, "{rep}_hbonds_intra_ntr_min_cutoff_{c}percent.list".format(rep=rep, c=5)), "w") as fo:
            for k, v in intrantr.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))

        extract = ["PHE277", "ARG282"]
        glyleu = extract_pair(intraarr, extract[0], extract[1])
        with open(os.path.join(cwd, "{}_hbonds_{}_{}.list".format(rep, extract[0], extract[1])), "w") as fo:
            for k, v in glyleu.items():
                fo.write("{k}\t{v1}\t{v2}\n".format(k=k, v1=v[0], v2="\t".join(v[1])))


if __name__ == "__main__":
    main()
