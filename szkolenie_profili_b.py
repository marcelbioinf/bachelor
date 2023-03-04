import pandas as pd
import statistics
import collections
import math
matrix = pd.read_csv('B_box_scoring_matrix.csv',  header = None, skiprows=1)
scoring_matrix = matrix.iloc[:, 1:]
#print(scoring_matrix)   #najpierw kolumna potem wiersz  wiersze od 0 do 10 a kolumny od 1 do 4

boxes_confirmed = []
with open("./ZZWYNIKI/b_box.txt") as file:
    for line in file:
        if line.startswith('>'):
            continue
        else:
            boxes_confirmed.append(line.upper().lstrip('SEQ: ').rstrip("\n"))

def estimate_shannon_entropy(dna_sequence):
    m = len(dna_sequence)
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])

    shannon_entropy_value = 0
    for base in bases:
        # number of residues
        n_i = bases[base]
        # n_i (# residues type i) / M (# residues in column)
        p_i = n_i / float(m)
        entropy_i = p_i * (math.log(p_i, 2))
        shannon_entropy_value += entropy_i

    return shannon_entropy_value * -1
shan = []
seqs = []
with open("./ZZWYNIKI/a_box_11.txt") as file:
    for line in file:
        if line.startswith('>'):
            #header = line.split()[1]
            continue
        else:
            seq = line.upper().rstrip("\n")
            seqs.append(seq)
file.close()
print(len(boxes_confirmed))
global_score = 0
scores = []
for seq in boxes_confirmed:
    score = 0
    shan.append(estimate_shannon_entropy(seq))
    print(estimate_shannon_entropy(seq))
    for i, nuc in enumerate(seq):
        if nuc == 'A':
            score += scoring_matrix[1][i]
        elif nuc == 'T':
            score += scoring_matrix[2][i]
        elif nuc == 'G':
            score += scoring_matrix[3][i]
        elif nuc == 'C':
            score += scoring_matrix[4][i]
    global_score += score
    scores.append(score)

print(global_score/2493)                         #srednia punktów to 4.49                             #UNIQUE 3.99
print(min([i for i in scores if i > 0]))         #min to -12 najmniejsza dodatnia to 0.07             0.04
print(min(scores))                                                                                     #-6.06
print(statistics.median(scores))                 #mediana to 4.98                                     4.28
print(max(scores))                               #max to 5.72                                          4.9
print(f'Shan: {sum(shan)/len(shan)}')

# seqs = []
# with open("./ZZWYNIKI/unique_trnasv2_nonrev.txt") as file:
#     for line in file:
#         if line.startswith('>'):
#             #header = line.split()[1]
#             continue
#         elif line.upper().startswith('SEQ: '):
#             seq = line.upper().lstrip('SEQ: ').rstrip("\n")
#             seqs.append(seq)
#         else:
#             continue
# file.close()
#
#
# boxes = []
#
# cutoff = -1     ##BAWIC SIE ZEBY ZNALEZC DOBRY
# treshold = 2.5
#
# for seq in seqs:
#     best_score = -1
#     box = ''
#     for i in range(len(seq)):
#         potential_box = seq[i:i+11]
#         print(estimate_shannon_entropy(potential_box))
#         if len(potential_box) < 11:
#             break
#         #print(potential_box)
#         score = 0
#         for j, nuc in enumerate(potential_box):
#             if nuc == 'A':
#                 score += scoring_matrix[1][j]
#             elif nuc == 'T':
#                 score += scoring_matrix[2][j]
#             elif nuc == 'G':
#                 score += scoring_matrix[3][j]
#             elif nuc == 'C':
#                 score += scoring_matrix[4][j]
#             if score < cutoff:
#                 break
#         #print(score)
#         if score > best_score:
#             best_score = score
#             box = potential_box
#         if best_score > treshold:
#             boxes.append(box)
#             break
#     # print(box)
#     #print(best_score)
#     #boxes.append(box)
#
#
# print(len(seqs))                                                                #2005 b boxow,   2546 tRNA, 2493 wykrytych b boxow(a tez w sumie 11 i 12)
# print(f'Znaleziono: {len(boxes)} aboxow')                            #1481 aboxow11 znaleziono we wszystkich tRNA
# if boxes == boxes_confirmed:
#     print("TUTTO BENE")
#
# def intersection(lst1, lst2):
#     lst3 = [value for value in lst1 if value in lst2]
#     return lst3
#
#
#
# print(f'Ilosc prawidlowych aboxow: {len(intersection(boxes, boxes_confirmed))}')             #2000 bboxow ktore znalazlem znajduje sie w abox confirmed, czyli 199 znalazlem nieprawidlowo