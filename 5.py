import re

out = open('termination_distance_100down.txt', 'w')
false = 0
disrupted = {'a': 0, 'g': 0, 'c': 0}
four_ts = []
four_ts_a = []
dis_five_23 = []
dis_five_32 = []
dis_five_322 = []
dis_five_14 = []
dis_four_13 = []
dis_four_31 = []
dis_four_22 = []
dis_four_222 = []

with open('./uniqe_trnas_genomic_seqs_v3_100down_nonrev.txt', 'r') as f:
    for line in f:
        if line.startswith('>'):
            out.write(f"{line.split()[0]}\n")
            pass
        else:
            downstream = line.rstrip('\n')[-100:].lower()
            match = re.search('tttt+', downstream)
            if match:
                #out.write(f"{downstream.find('tttt') + 1}\n")
                four_ts.append(downstream.find('tttt') + 1)
            elif re.search(r'tt[agc]ttt', downstream):
                x = re.search(r'tt[agc]ttt', downstream).group(0)
                dis_five_23.append(downstream.find(str(x))+1)
                disrupted[x[2]] += 1
            elif re.search(r'ttt[agc]tt', downstream):
                x = re.search(r'ttt[agc]tt', downstream).group(0)
                dis_five_32.append(downstream.find(str(x))+1)
                disrupted[x[3]] += 1
            elif re.search(r'ttt[agc][agc]tt', downstream):
                x = re.search(r'ttt[agc][agc]tt', downstream).group(0)
                dis_five_322.append(downstream.find(str(x))+1)
                disrupted[x[3]] += 1
                disrupted[x[4]] += 1
            elif re.search(r't[agc]tttt', downstream):
                x = re.search(r't[agc]tttt', downstream).group(0)
                dis_five_14.append(downstream.find(str(x))+1)
                disrupted[x[1]] += 1
            elif re.search(r't[agc]ttt', downstream):
                x = re.search(r't[agc]ttt', downstream).group(0)
                dis_four_13.append(downstream.find(str(x))+1)
                disrupted[x[1]] += 1
            elif re.search(r'ttt[agc]t', downstream):
                x = re.search(r'ttt[agc]t', downstream).group(0)
                dis_four_31.append(downstream.find(str(x))+1)
                disrupted[x[3]] += 1
            elif re.search(r'tt[agc]tt', downstream):
                x = re.search(r'tt[agc]tt', downstream).group(0)
                dis_four_22.append(downstream.find(str(x))+1)
                disrupted[x[2]] += 1
            elif re.search(r'tt[agc][agc]tt', downstream):
                x = re.search(r'tt[agc][agc]tt', downstream).group(0)
                dis_four_222.append(downstream.find(str(x))+1)
                disrupted[x[2]] += 1
                disrupted[x[3]] += 1
            else:
                false += 1

print(false)
print(len(four_ts))     #100 downstream - 513 bez 4 lub więcej T, 102 bez 3 lub wiecej T, 371 bez 4 lub wiecej T lub TTATT+, 80 bez 3 lub wiecej T lub TTATT+
print(len(dis_five_23))
print(len(dis_five_32))
print(len(dis_five_322))
print(len(dis_five_14))
print(len(dis_four_13))
print(len(dis_four_31))
print(len(dis_four_22))
print(len(dis_four_222))
print(disrupted)