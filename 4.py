save = open('ZZWYNIKI/POPRAWKI/50_downstreamv1.txt', 'w')

with open('ZZWYNIKI/POPRAWKI/all_trnas_genomic_seqs_v1_50downup_nonrev_bothstrands.txt', 'r') as f:
    for line in f:
        if line.startswith('>'):
            header = line.split()
            save.write(f'{header[0]}\n')
        else:
            if not len(line.rstrip('\n')) < 171:
                save.write(f'{line[-51:].upper()}')

