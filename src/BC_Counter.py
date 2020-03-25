import os
from sys import argv
script, input_file = argv

#open file

print('Opening...')
filename = input('Enter file name with fq extension to analyze >: ')
handle = open(os.path.expanduser(f"~/Desktop/{filename}"))
print('File opened')

barcodes = open(input_file)

BC_list = dict()

for line in barcodes:
    line.rstrip()
    y = line.split()
    if len(line) > 12:
        BC_list[y[1]] = y[2]

print(BC_list)

inputBC = input('Press ENTER to confirm Barcode list: > ')

trimmed = list()

trimming_count = 0

#change number for analyzing up to this # of reads

readstocount = int(input('Enter a max number of reads to count >: '))

#this loop trims the 72 first nucleotides of the sequence

for line in handle:
    line.rstrip()
    if line.startswith('GCTGG'):
        trimmed.append(line[72:78])
        trimming_count = trimming_count + 1
        print('Sequences trimmed:', round(100*(trimming_count/readstocount),2),'%')
        if trimming_count >= readstocount:
            break
    else:
        continue

print('Total sequences trimmed:', trimming_count)

print('Trimming Completed, now analyzing barcodes')
print()
print('**********')
print()

counts = dict()

#this loop generates an histogram of the frequency of each barcode

for barcode in trimmed:
    counts[barcode] = counts.get(barcode, 0) + 1

final = counts.items()

count_reads = 0

final_reads = list()

BC_final = BC_list.items()

#this loop matches AAV with Barcode Counts

for bc,aav in BC_final:
    for k,v in final:
        if k == bc:
            count_reads = count_reads + v
            final_reads.append([aav, v])
        else:
            continue


print(f'Total valid reads: {count_reads} of {trimming_count} -->', round(100*count_reads/trimming_count,2),'%')
print()


final_readss = sorted(final_reads)

#this loop prints AAV Barcode Counts and Percentage


print(filename)
print('*****')
for sd, fd in final_readss:
    print(sd, 'Total reads:', fd)
    print(fd)


print('Thanks for your patience, mate')
