import re
from itertools import groupby
numms = {}
matches = {}
frequencies = {}
false = 0
caas = []
caas_f = []
both = 0
tata_n = 0
with open('ZZWYNIKI/POPRAWKI/50_upstreamv1.txt', 'r') as f:
    for line in f:
        if not line.startswith(">"):
            flag = False
            upstream = line.rstrip('\n').lower()
            #print(list(map(seq.count, "ATGC"))) #1 sposób
            nucs = dict(zip(['A', 'T', 'G', 'C'], list(map(upstream.count, "atgc"))))
            for key, value in nucs.items():
                frequencies[key] = frequencies.get(key, 0) + (value / len(upstream))
            match = re.search('tataaaa', upstream)
            if match:
                matches[match.group(0)] = matches.get(match.group(0), 0) + 1
            elif re.search('tatataa', upstream):
                match = re.search('tatataa', upstream)
                matches[match.group(0)] = matches.get(match.group(0), 0) + 1
            elif re.search('tataaat', upstream):
                match = re.search('tataaat', upstream)
                matches[match.group(0)] = matches.get(match.group(0), 0) + 1
            elif re.search('tatatat', upstream):
                match = re.search('tatatat', upstream)
                matches[match.group(0)] = matches.get(match.group(0), 0) + 1
            elif re.search("tata", upstream):
                match = re.search("tata", upstream)
                matches[match.group(0)] = matches.get(match.group(0), 0) + 1
                tata_n += 1
            else:
                false += 1
            # match = re.search("tata[tcga][tcga][tcga]", upstream)
            # if match:
            #     if match.group(0) not in matches.keys():
            #         matches[match.group(0)] = matches.get(match.group(0), 0) + 1
            #         false -= 1
            pos = upstream.rfind('caa') + 1
            if pos > 0:
                caas.append(pos)
                flag = True
            pos_f = upstream.find('caa') + 1
            if pos_f > 0:
                print(pos_f)
                caas_f.append(pos_f)
                if flag and pos != pos_f:
                    both += 1
            numm = upstream.count('caa')
            numms[numm] = numms.get(numm, 0) + 1

print(false)  #1533   5399
print(matches)  #{'tata': 666, 'tatatat': 64, 'tataaat': 64, 'tatataa': 87, 'tataaaa': 132}  {'tata': 2195, 'tatatat': 192, 'tataaat': 207, 'tatataa': 254, 'tataaaa': 463}
for key in frequencies.keys():
    frequencies[key] = round(frequencies[key] / 8720, 4)
print(frequencies)  #{'A': 0.422, 'T': 0.2949, 'G': 0.1504, 'C': 0.1322}
print(len(caas_f)) #1975 - ma jeden lub wiecej motyw CAA, 571 nie ma wcale, 1094 ma 2 lub wiecej  6820 ma jeden lub wiecej,  2195 nie ma
print(tata_n) #649 tata...

freqs = {}
for elem in caas:
    if not elem in freqs.keys():
        freqs[elem] = caas.count(elem)

list_a = []
list_b = []
list_c = []
list_d = []
for key, value in freqs.items():
    list_a.append(key)
    list_b.append(value)

print(f'Ostatnie CAA pozycje: {list_a}\n')
print(f'Ostatnie CAA liczebnosc: {list_b}\n')

for i in range(len(list_b)):
    if list_a[i] >= 35:
        list_c.append(list_a[i])
        list_d.append(list_b[i])


freqs_f = {}
for elem in caas_f:
    if not elem in freqs_f.keys():
        freqs_f[elem] = caas_f.count(elem)

list_af = []
list_bf = []

for key, value in freqs_f.items():
    list_af.append(key)
    list_bf.append(value)

print(f'Pierwsze CAA pozycje: {list_af}\n')
print(f'Pierwsze CAA liczebnosc: {list_bf}\n')

print(f'Liczba sekwencji z 2 lub wiecej motywami CAA: {both}')

print(f'Rozklad motywów CAA: {numms}')

print(list_c)
print(list_d)