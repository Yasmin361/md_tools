import os, glob, pandas


def main():
    cwd = os.getcwd()
    sum_mdc_interface(cwd)
    sum_mdc_pairs(cwd)

def sum_mdc_pairs(cwd):
    interface = list()
    for mdcres in glob.glob(os.path.join(cwd, "**/*.xlsx"), recursive=True):
        contacts = pandas.read_excel(mdcres, sheet_name=0, header=1)
        interface.append(contacts)
    pairlist = dict()
    for replica in interface:
        lines = replica.to_numpy()
        for line in lines:
            line = str(line).lstrip('[').rstrip(']')
            line = line.split()
            if len(line[0].split("@")) == 2:                #pair-contact entry
                pair = "{p1} - {p2}".format(p1=line[0].lstrip("'"), p2=line[2].rstrip("'")) # only write residue name
                if pair not in pairlist:
                    pairlist[pair] = float(line[3].lstrip("'").rstrip("'"))
                    continue
                pairlist[pair] += float(line[3].lstrip("'").rstrip("'"))
    for k in pairlist:
        pairlist[k] = str(round(pairlist[k], 2))

    with open(os.path.join(cwd, "sum_of_replicas_pairs.list"), "w") as fo:
        for item in pairlist:
            fo.write("{k}\t{v}\n".format(k=item, v=pairlist[item]))


def sum_mdc_interface(cwd):
    interface = list()
    for mdcres in glob.glob(os.path.join(cwd, "**/*.xlsx"), recursive=True):
        contacts = pandas.read_excel(mdcres, sheet_name=1, header=1)
        interface.append(contacts)
    contactlist_arr = dict()
    contactlist_ntr = dict()
    for replica in interface:
        lines = replica.to_numpy()
        for line in lines:
            line = str(line).lstrip('[').rstrip(']')
            line = line.split()
            if len(line[0].split("@")) == 2:                #frag0 = arrestin
                frag0 = line[0].lstrip("'").rstrip("'").split("@")[0] # only write residue name
                if frag0 not in contactlist_arr:
                    contactlist_arr[frag0] = float(line[1].lstrip("'").rstrip("'"))
                    continue
                contactlist_arr[frag0] += float(line[1].lstrip("'").rstrip("'"))
            if len(line[3].split("@")) == 2:               # frag1 = receptor
                frag1 = line[3].lstrip("'").rstrip("'").split("@")[0]
                if frag1 not in contactlist_ntr:
                    contactlist_ntr[frag1] = float(line[4].lstrip("'").rstrip("'"))
                    continue
                contactlist_ntr[frag1] += float(line[4].lstrip("'").rstrip("'"))
    for k in contactlist_arr:
        contactlist_arr[k] = str(round(contactlist_arr[k], 2))
    for k in contactlist_ntr:
        contactlist_ntr[k] = str(round(contactlist_ntr[k], 2))

    with open(os.path.join(cwd, "sum_of_replicas_arr_contacts.list"), "w") as fo:
        for item in contactlist_arr:
            fo.write("{k}\t{v}\n".format(k=item, v=contactlist_arr[item]))

    with open(os.path.join(cwd, "sum_of_replicas_ntr_contacts.list"), "w") as fo:
        for item in contactlist_ntr:
            fo.write("{k}\t{v}\n".format(k=item, v=contactlist_ntr[item]))


if __name__ == "__main__":
    main()
