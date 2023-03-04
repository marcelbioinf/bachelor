import pandas as pd
import re
import math

a_boxes = []   #A            #T            #G             #C
frequencies = {0: 0.31224, 1: 0.31226, 2: 0.18772, 3: 0.18778}
with open("./ZZWYNIKI/a_box_11.txt") as file:
    for line in file:
        if line.startswith('>'):
            continue
        else:
            a_boxes.append(line.rstrip("\n"))
file.close()

print(len(a_boxes))   #1387 abox11  1106 abox12  2493 bbox
unique_a = set(a_boxes)
print(len(unique_a))      #102   129     266

profile_abox = [[0 for i in range(4)] for j in range(11)]

for seq in a_boxes:
    for i, nuc in enumerate(seq):
        if nuc == 'A':
            profile_abox[i][0] += 1
        elif nuc == 'T':
            profile_abox[i][1] += 1
        elif nuc == 'G':
            profile_abox[i][2] += 1
        elif nuc == 'C':
            profile_abox[i][3] += 1

profile_abox = [[round(i / len(a_boxes), 3) for i in profile_abox[j]] for j in range(11)]

for i in range(len(profile_abox)):
    for j in range(len(profile_abox[i])):
        if profile_abox[i][j] == 0.0:
            profile_abox[i][j] += 0.0001
        profile_abox[i][j] = round(math.log(profile_abox[i][j]/frequencies[j],10),2)


for row in profile_abox:
    print(f"{row}\n")



pd.DataFrame(profile_abox).to_csv("./A_box11_scoring_matrix.csv")