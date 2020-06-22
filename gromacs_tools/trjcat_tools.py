#this generates the gmx trjcat command for a span of parts generated with sub.sh
def trjcatspan(x,y):
    first_part = x
    last_part = y
    select_parts = range(first_part, last_part+1)
    part_list = []
    for i in select_parts:
        partnumber = "step7_production.part"+str(i).zfill(4)+".xtc"
        part_list.append(partnumber)
    return "gmx trjcat -f "+" ".join(part_list)+" -o step7_production.cat"+str(first_part)+"_"+str(last_part)+".xtc"


print(trjcatspan(75,83))
