import os, glob, pandas, re


def main():
    cwd = os.getcwd()
    sum_mdc_interface(cwd)


def sum_mdc_interface(cwd):
    interface = list()
    for mdcres in glob.glob(os.path.join(cwd, "**/*.xlsx"), recursive=True):
        contacts = pandas.read_excel(mdcres, sheet_name=1, header=1)
        interface.append(contacts)
    contactlist = dict()
    for replica in interface:
        lines = replica.to_numpy()
        for line in lines:
            line = str(line).lstrip('[').rstrip(']')
            line = line.split()
            if len(line[0].split("@")) == 2:                #frag0
                frag0 = line[0].lstrip("'").rstrip("'")
                if frag0 not in contactlist:
                    contactlist[frag0] = float(line[1].lstrip("'").rstrip("'"))
                    continue
                contactlist[frag0] += float(line[1].lstrip("'").rstrip("'"))
            if len(line[3].split("@")) == 2:               # frag1
                frag1 = line[3].lstrip("'").rstrip("'")
                if frag1 not in contactlist:
                    contactlist[frag1] = float(line[4].lstrip("'").rstrip("'"))
                    continue
                contactlist[frag1] += float(line[4].lstrip("'").rstrip("'"))
    for k in contactlist:
        contactlist[k] = str(round(contactlist[k], 2))
    with open(os.path.join(cwd, "sum_of_replicas.list"), "w") as fo:
        for item in contactlist:
            fo.write("{k}\t{v}\n".format(k=item, v=contactlist[item]))


if __name__ == "__main__":
    main()
