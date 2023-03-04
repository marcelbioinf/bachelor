import os

for file in os.listdir():
    if file.endswith('fna'):
         os.system(f"tRNAscan-SE -E -o {file.rstrip('_g.fna')}_result.txt -f {file.rstrip('_g.fna')}result_ss.txt --breakdown --detail {file}")


