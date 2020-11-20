#this generates the gmx trjcat command for a span of parts generated with sub.sh
def trjcatspan(x,y):
    first_part = x
    last_part = y
    select_parts = range(first_part, last_part+1)
    part_list = []
    for i in select_parts:
        partnumber = "step7_production.part"+str(i).zfill(4)+".xtc"
        part_list.append(partnumber)
    #part_list = ["step7_production.part"+str(i).zfill(4)+".xtc" for i in select_parts if i != None]
    return "gmx trjcat -f "+" ".join(part_list)+" -o step7_production.cat"+str(first_part)+"_"+str(last_part)+".xtc"

print(trjcatspan(75, 80))



#T187 - H198 + T221-K230 + T261-T267+ V328-D342

#this generates index group sequence for merging operation in gmx make_ndx
def ndx_groups(x, y):
    sequence = range(x, y+1)
    index_list = []
    for i in sequence:
        index_list.append(str(i))

    return " | ".join(index_list)


print(ndx_groups(1, 9))






def ndx_multigroups(x, z, o):
    index_list = []
    for i in range(0,z):
        sequence = x+(i*o)
        index_list.append(str(sequence))

    return " | ".join(index_list)


#print(ndx_multigroups(20, 311, 134))


