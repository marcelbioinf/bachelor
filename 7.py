import re

in_f = open("./ZZWYNIKI/t_d_loopsv3.txt")

#text = in_f.readlines()
#print(text)

Dlengths, Tlenghts = {}, {}

for line in in_f:
    if line.startswith("D-loop: "):
        seq = line.lstrip("D-loop: ").rstrip("\n")
        if len(seq) in Dlengths.keys():
            Dlengths[len(seq)] += 1
        else:
            Dlengths[len(seq)] = 1
    elif line.startswith("T-loop: "):
        #seq = line.rstrip("\n").lstrip("T-loop: ")  #UCINA ZA DUZO BO WCINA TEZ 2 T Z BRZEGU
        seq = re.sub(r"T-loop: ", "", line.rstrip("\n"))
        if len(seq) == 6:
            print(seq)
        if len(seq) in Tlenghts.keys():
            Tlenghts[len(seq)] += 1
        else:
            Tlenghts[len(seq)] = 1


print(Dlengths)
print(Tlenghts)