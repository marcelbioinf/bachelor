import pandas as pd
import statistics
import collections
import math

matrix_a11 = pd.read_csv('A_box11_scoring_matrix.csv',  header = None, skiprows=1)
scoring_a11_matrix = matrix_a11.iloc[:, 1:]

matrix_a12 = pd.read_csv('A_box12_scoring_matrix.csv',  header = None, skiprows=1)
scoring_a12_matrix = matrix_a12.iloc[:, 1:]

matrix_b = pd.read_csv('B_box_scoring_matrix.csv',  header = None, skiprows=1)
scoring_b_matrix = matrix_b.iloc[:, 1:]

# boxes_confirmed_b = []
# with open("./ZZWYNIKI/b_box.txt") as file:
#     for line in file:
#         if line.startswith('>'):
#             continue
#         else:
#             boxes_confirmed_b.append(line.upper().lstrip('SEQ: ').rstrip("\n"))

# boxes_confirmed_a11 = []
# with open("./ZZWYNIKI/a_box_11.txt") as file:
#     for line in file:
#         if line.startswith('>'):
#             continue
#         else:
#             boxes_confirmed_a11.append(line.upper().lstrip('SEQ: ').rstrip("\n"))

# boxes_confirmed_a12 = []
# with open("./ZZWYNIKI/a_box_12.txt") as file:
#     for line in file:
#         if line.startswith('>'):
#             continue
#         else:
#             boxes_confirmed_a12.append(line.upper().lstrip('SEQ: ').rstrip("\n"))

# boxes_arabidopsis_a12 = []
# with open("./atabidobsis_a12_confirmed.txt") as file:
#     for line in file:
#         boxes_arabidopsis_a12.append(line.upper().rstrip("\n"))
#
# boxes_arabidopsis_a11 = []
# with open("./atabidobsis_a11_confirmed.txt") as file:
#     for line in file:
#         boxes_arabidopsis_a11.append(line.upper().rstrip("\n"))
#
# boxes_arabidopsis_b = []
# with open("./atabidobsis_b_confirmed.txt") as file:
#     for line in file:
#         boxes_arabidopsis_b.append(line.upper().rstrip("\n"))

#to na później żeby wszystkie pliki ogarnąć
# files = []
# for directory in os.listdir():
#     if os.path.isdir(f'./{directory}'):
#         if directory.endswith('idea'):
#             continue
#         for file in os.listdir(f'./{directory}'):
#             if file.endswith('.fna'):
#                 files.append(f'{directory}/{file}')
trnas = 0
trnas_list = []
headers_list = []
boxes_a11 = []
boxes_a12 = []
boxes_b = []
cutoff = -1
treshold_a = 3.4
treshold_b = 2.5
print(f'Abox treshold: {treshold_a}\nBbox treshold: {treshold_b}')


with open('./arabidopsis_thaliana/arabidopsis_thaliana_g.fna') as file:
    for chromosome in file.read().split('>'):
        if chromosome:
            header, *content = chromosome.splitlines()
            seq = ''.join(content).upper()
            seq = seq.replace("\n", "")
            print(seq)
            a11_box = ''
            a12_box = ''
            b_box = ''
            is_boxb = False
            is_boxa11 = False
            is_boxa12 = False
            a_box_ending_pos = 0
            b_box_ending_pos = 0
            seq_len, i = len(seq), 0
            print(f'Długosc sekwencji {seq_len}\n')
            while i < seq_len:
                #print(f'Pozycja {i}')
                potential_b_box = seq[i:i + 11]
                #print(f'Entr: {estimate_shannon_entropy(potential_b_box)}')
                if len(potential_b_box) < 11:
                    break
                score_b = 0
                for j, nuc in enumerate(potential_b_box):
                    if nuc == 'A':
                        score_b += scoring_b_matrix[1][j]
                    elif nuc == 'T':
                        score_b += scoring_b_matrix[2][j]
                    elif nuc == 'G':
                        score_b += scoring_b_matrix[3][j]
                    elif nuc == 'C':
                        score_b += scoring_b_matrix[4][j]
                    if score_b < cutoff:
                        break
                if score_b > treshold_b:
                    b_box = potential_b_box
                    boxes_b.append(b_box)
                    b_box_ending_pos = i + 11
                    is_boxb = True
                    #print(f'Znaleziono boxa b: {b_box} punkty: {score_b} pozycja: {i} pozycja konca: {b_box_ending_pos}')

                if is_boxb:
                    a_box_searching_seq = seq[b_box_ending_pos-64:b_box_ending_pos-30]
                    #print(f'    Abox searching seq : {a_box_searching_seq}')
                    for k in range(len(a_box_searching_seq)):
                        potential_a11_box = a_box_searching_seq[k:k + 11]
                        potential_a12_box = a_box_searching_seq[k:k + 12]
                        #print(f'     Potencjalny abox11: {potential_a11_box}')
                        #print(f'     Potencjalny abox12: {potential_a12_box}')
                        if len(potential_a11_box) < 11:
                            break
                        score_a11 = 0
                        score_a12 = 0
                        for j, nuc in enumerate(potential_a12_box):
                            if j == 11:
                                if nuc == 'A':
                                    score_a12 += scoring_a12_matrix[1][j]
                                elif nuc == 'T':
                                    score_a12 += scoring_a12_matrix[2][j]
                                elif nuc == 'G':
                                    score_a12 += scoring_a12_matrix[3][j]
                                elif nuc == 'C':
                                    score_a12 += scoring_a12_matrix[4][j]
                                break
                            if nuc == 'A':
                                score_a11 += scoring_a11_matrix[1][j]
                                score_a12 += scoring_a12_matrix[1][j]
                            elif nuc == 'T':
                                score_a11 += scoring_a11_matrix[2][j]
                                score_a12 += scoring_a12_matrix[2][j]
                            elif nuc == 'G':
                                score_a11 += scoring_a11_matrix[3][j]
                                score_a12 += scoring_a12_matrix[3][j]
                            elif nuc == 'C':
                                score_a11 += scoring_a11_matrix[4][j]
                                score_a12 += scoring_a12_matrix[4][j]
                            if score_a11 < cutoff and score_a12 < cutoff:
                                break
                        if score_a11 > treshold_a and score_a11 > score_a12:
                            a11_box = potential_a11_box
                            boxes_a11.append(a11_box)
                            a_box_ending_pos = b_box_ending_pos - 64 + k + 11
                            is_boxa11 = True
                            #print(f'Znaleziono boxa a11: {a11_box} punkty: {score_a11} pozycja: {b_box_ending_pos - 64 + k} pozycja konca: {a_box_ending_pos}')
                        elif score_a12 > treshold_a and score_a12 > score_a11:
                            a12_box = potential_a12_box
                            boxes_a12.append(a12_box)
                            a_box_ending_pos = b_box_ending_pos - 64 + k + 12
                            is_boxa12 = True
                            #print(f'Znaleziono boxa a12: {a12_box} punkty: {score_a12} pozycja: {b_box_ending_pos - 64 + k} pozycja konca: {a_box_ending_pos}')
                        #else:
                            #print(f'        score potemcjalnego a-boxa11: {score_a11}')
                            #print(f'        score potemcjalnego a-boxa12: {score_a12}')
                    if is_boxa11 == False and is_boxa12 == False:
                        i += 1
                        is_boxb = False
                else:
                    i += 1

                if is_boxb and (is_boxa11 or is_boxa12):
                    trnas += 1
                    h = header.split()
                    if is_boxa11:
                        trna = seq[a_box_ending_pos - 18:b_box_ending_pos + 11]
                        headers_list.append(f">{' '.join(h[0:3])} tRNA {trnas} {a_box_ending_pos - 18}:{b_box_ending_pos + 11} {(b_box_ending_pos+11)-(a_box_ending_pos-18)}bp")
                        print(f'!!!Znaleziony tRNA - {trna} \n!!!Pozycje: {a_box_ending_pos - 18, b_box_ending_pos + 11} {(b_box_ending_pos+11)-(a_box_ending_pos-18)}bp')
                    if is_boxa12:
                        trna = seq[a_box_ending_pos - 19:b_box_ending_pos + 11]
                        headers_list.append(f">{' '.join(h[0:3])} tRNA {trnas} {a_box_ending_pos - 19}:{b_box_ending_pos + 11} {(b_box_ending_pos+11)-(a_box_ending_pos-19)}bp")
                        print(f'!!!Znaleziony tRNA - {trna} \n!!!Pozycje: {a_box_ending_pos - 19, b_box_ending_pos + 11} {(b_box_ending_pos+11)-(a_box_ending_pos-19)}bp')
                    trnas_list.append(trna)
                    is_boxa11, is_boxa12, is_boxb = False, False, False
                    a11_box = ''
                    a12_box = ''
                    b_box = ''
                    i = b_box_ending_pos + 1
                    a_box_ending_pos = 0
                    b_box_ending_pos = 0
                    print(f"Znaleziono: {trnas} tRNA dotychczas")
            print('Kolejny chromosom!!!!!\n')


print(f'Liczba tRNA zakwalifikowanych jako tRNA wg profilu: {trnas}\n')

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

print(f'Ilosc aboxów11 prawdziwych: {len(boxes_confirmed_a11)}')
print(f'Ilosc aboxów11 prawdziwych arabidopsis: {len(boxes_arabidopsis_a11)}')
print(f'Ilosc boxow a11 znalezionyc wg profilu {len(boxes_a11)}')
print(f'Ilosc prawidlowych aboxow_11: {len(intersection(boxes_a11, boxes_confirmed_a11))}')
print(f'Ilosc prawidlowych aboxow_11 z arabidopsis: {len(intersection(boxes_a11, boxes_arabidopsis_a11))}\n')

print(f'Ilosc aboxów12 prawdziwych: {len(boxes_confirmed_a12)}')
print(f'Ilosc aboxów12 prawdziwych arabidopsis: {len(boxes_arabidopsis_a12)}')
print(f'Ilosc boxow a12 znalezionyc wg profilu {len(boxes_a12)}')
print(f'Ilosc prawidlowych aboxow_12: {len(intersection(boxes_a12, boxes_confirmed_a12))}')
print(f'Ilosc prawidlowych aboxow_12 z arabidopsis: {len(intersection(boxes_a12, boxes_arabidopsis_a12))}\n')

print(f'Ilosc bboxow prawdziwych: {len(boxes_confirmed_b)}')
print(f'Ilosc b boxow prawdziwych arabidopsis: {len(boxes_arabidopsis_b)}')
print(f'Ilosc boxow b znalezionyc wg profilu {len(boxes_b)}')
print(f'Ilosc prawidlowych bboxow: {len(intersection(boxes_b, boxes_confirmed_b))}')
print(f'Ilosc prawidlowych bboxow z arabidopsis: {len(intersection(boxes_b, boxes_arabidopsis_b))}\n')

# out1 = open('./profile/abxoes11_arabidopsis_legit_fasterv1.txt', 'w')
# for item in boxes_a11:
#     out1.write(f'{item}\n')
#
# out2 = open('./profile/abxoes12_arabidopsis_legit_fasterv1.txt', 'w')
# for item in boxes_a12:
#     out2.write(f'{item}\n')
#
# out3 = open('./profile/bboxes_arabidopsis_legit_fasterv1.txt', 'w')
# for item in boxes_b:
#     out3.write(f'{item}\n')

out4 = open('./profile/trnas_arabidopsisch1_legit_fasterULTIMATEV3.txt', 'w')
for i, item in enumerate(trnas_list):
    out4.write(f'{headers_list[i]}\n')
    out4.write(f'{item}\n')
out4.close()