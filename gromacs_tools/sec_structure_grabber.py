def getsecstructure(x):
    name = x.split(".")                                 # get pdb code
    with open(x, "r") as pdb:
        print("Reading " + pdb.name)
        line = pdb.readline()
        helices = list()
        sheets = list()
        while line:
            if line.split()[0] == "HELIX":              # helix record
                hrecord = line.split()
                chain = hrecord[4]
                startres = hrecord[5]
                stopres = hrecord[8]
                helices.append([chain, startres, stopres])
                line = pdb.readline()
                continue
            if line.split()[0] == "SHEET":              # sheet record
                srecord = line.split()
                chain = srecord[5]
                startres = srecord[6]
                stopres = srecord[9]
                sheets.append([chain, startres, stopres])
                line = pdb.readline()
                continue
            line = pdb.readline()
    with open(name[0]+"_helices.csv", "w") as outf:     # write helices
        for helix in helices:
            outf.write(", ".join(helix) + "\n")
    with open(name[0]+"_sheets.csv", "w") as outf:      # write sheets
        for sheet in sheets:
            outf.write(", ".join(sheet) + "\n")


getsecstructure("6up7.pdb")


