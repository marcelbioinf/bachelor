import pandas as pd
import statistics
matrix = pd.read_csv('A_box12_scoring_matrix.csv',  header = None, skiprows=1)
scoring_matrix = matrix.iloc[:, 1:]
print(scoring_matrix)   #najpierw kolumna potem wiersz  wiersze od 0 do 10 a kolumny od 1 do 4

boxes_confirmed = []
with open("./ZZWYNIKI/a_box_12.txt") as file:
    for line in file:
        if line.startswith('>'):
            continue
        else:
            boxes_confirmed.append(line.upper().lstrip('SEQ: ').rstrip("\n"))

print(len(boxes_confirmed))

# seqs = []
# with open("./ZZWYNIKI/a_box_11.txt") as file:
#     for line in file:
#         if line.startswith('>'):
#             #header = line.split()[1]
#             continue
#         else:
#             seq = line.upper().rstrip("\n")
#             seqs.append(seq)
# file.close()
# global_score = 0
# scores = []
# for seq in boxes_confirmed:
#     score = 0
#     for i, nuc in enumerate(seq):
#         if nuc == 'A':
#             score += scoring_matrix[1][i]
#         elif nuc == 'T':
#             score += scoring_matrix[2][i]
#         elif nuc == 'G':
#             score += scoring_matrix[3][i]
#         elif nuc == 'C':
#             score += scoring_matrix[4][i]
#     global_score += score
#     scores.append(score)
#
# print(global_score/1106)                         #srednia punktów to 5.2502                              #UNIQUE 4.49
# print(min([i for i in scores if i > 0]))       #min to -9.3 najmniejsza dodatnia to 0.12                -3.15  #0.23
# print(min(scores))
# print(statistics.median(scores))                 #mediana to 6.14                                       5.01
# print(max(scores))                               #max to 6.56                                           5.2


seqs = []
with open("./ZZWYNIKI/unique_trnasv2_nonrev.txt") as file:
    for line in file:
        if line.startswith('>'):
            #header = line.split()[1]
            continue
        elif line.upper().startswith('SEQ: '):
            seq = line.upper().lstrip('SEQ: ').rstrip("\n")
            seqs.append(seq)
        else:
            continue
file.close()


boxes = []

cutoff = -1     ##BAWIC SIE ZEBY ZNALEZC DOBRY
treshold = 4.1

for seq in seqs:
    best_score = -1
    box = ''
    for i in range(len(seq)):
        potential_box = seq[i:i+12]
        if len(potential_box) < 12:
            break
        #print(potential_box)
        score = 0
        for j, nuc in enumerate(potential_box):
            if nuc == 'A':
                score += scoring_matrix[1][j]
            elif nuc == 'T':
                score += scoring_matrix[2][j]
            elif nuc == 'G':
                score += scoring_matrix[3][j]
            elif nuc == 'C':
                score += scoring_matrix[4][j]
            if score < cutoff:
                break
        #print(score)
        if score > best_score:
            best_score = score
            box = potential_box
        if best_score > treshold:
            boxes.append(box)
            break
    # print(box)
    #print(best_score)
    #boxes.append(box)

#print(len(seqs))                                                                #1106 a boxow12,   2546 tRNA
print(f'Znaleziono: {len(boxes)} aboxow')                            #921 aboxow11 znaleziono we wszystkich tRNA
if boxes == boxes_confirmed:
    print("TUTTO BENE")

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3



print(f'Ilosc prawidlowych aboxow: {len(intersection(boxes, boxes_confirmed))}')             #916 aboxw ktore znalazlem znajduje sie w abox confirmed