import argparse
import pysam

def classifier(regions,save_path):
    for region in regions:
        chrom, start2, end1, start1, end2, pos = region  
        class1 = 0
        class2 = 0
        class3 = 0
        total = 0
        threshot = 30

        for read in bamfile.fetch(chrom):
            if read.query_alignment_sequence is None:
                continue
            tag = 0
            
            if int(start1) - int(threshot) <= read.reference_start <= int(start1) + int(threshot) and int(end1) - int(threshot) <= read.reference_end <= int(end1) + int(threshot):
                class1 += 1
                tag += 1
                # print(1)
            elif int(start1) - int(threshot) <= read.reference_start <= int(start1) + int(threshot) and int(end2) - int(threshot) <= read.reference_end <= int(end2) + int(threshot):
                class1 += 1
                tag += 1
                # print(2)
            elif int(start2) - int(threshot) <= read.reference_start <= int(start2) + int(threshot) and int(end1) - int(threshot) <= read.reference_end <= int(end1) + int(threshot):
                class1 += 1
                tag += 1
            elif int(start2) - int(threshot) <= read.reference_start <= int(start2) + int(threshot) and int(end2) - int(threshot) <= read.reference_end <= int(end2) + int(threshot):
                class1 += 1
                tag += 1
            elif int(start2) + int(threshot) < read.reference_start < int(start1) - int(threshot) or int(end1) + int(threshot) < read.reference_end < int(end2) - int(threshot):
                class1 += 1
                tag += 1

            if int(start1) - int(threshot) <= read.reference_start <= int(start1) + int(threshot) and read.reference_end > int(end2) + int(threshot):
                class2 += 1
                tag += 1
            elif int(start1) - int(threshot) <= read.reference_start <= int(start1) + int(threshot) and read.reference_end > int(end1) + int(threshot) and read.reference_end < int(end2) - int(threshot):
                class2 += 1
                tag += 1
            elif int(start2) - int(threshot) <= read.reference_start <= int(start2) + int(threshot) and read.reference_end > int(end2) + int(threshot):
                class2 += 1
                tag += 1
            elif int(start2) - int(threshot) <= read.reference_start <= int(start2) + int(threshot) and read.reference_end > int(end1) + int(threshot) and read.reference_end < int(end2) - int(threshot):
                class2 += 1
                tag += 1
            elif read.reference_start < int(start2) - int(threshot) and int(end1) - int(threshot) <= read.reference_end <= int(end1) + int(threshot):
                class2 += 1
                tag += 1
            elif int(start2) + int(threshot) < read.reference_start < int(start1) - int(threshot) and int(end1) - int(threshot) <= read.reference_end <= int(end1) + int(threshot):
                class2 += 1
                tag += 1
            elif read.reference_start < int(start2) - int(threshot) and int(end2) - int(threshot) <= read.reference_end <= int(end2) + int(threshot):
                class2 += 1
                tag += 1
            elif int(start2) + int(threshot) < read.reference_start < int(start1) - int(threshot) and int(end2) - int(threshot) <= read.reference_end <= int(end2) + int(threshot):
                class2 += 1
                tag += 1

            
            if int(start1) + int(threshot) < read.reference_start < int(end1) and int(end1) - int(threshot) < read.reference_end < int(end1) + int(threshot):
                class3 += 1
                tag += 1
            elif int(start1) - int(threshot) < read.reference_start < int(start1) + int(threshot) and read.reference_end < int(end1) - int(threshot):
                class3 += 1
                tag += 1
            elif int(start1) + int(threshot) < read.reference_start and read.reference_end < int(end1) - int(threshot):
                class3 += 1
                tag += 1
            elif int(start2) - int(threshot) < read.reference_start < int(start2) + int(threshot) and read.reference_end < int(end1) - int(threshot):
                class3 += 1
                tag += 1
            elif int(start1) + int(threshot) < read.reference_start and int(end2) - int(threshot) < read.reference_end < int(end2) + int(threshot):
                class3 += 1
                tag += 1
            total += 1

        with open(save_path,'a') as f:
            f.write(f'{pos}:  class1:{class1}   class2:{class2}   class3:{class3}   mapped:{class1 + class2 + class3}   total:{total}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--bam','-b',type=str,help='bam path',required=True)
    parser.add_argument('--bed',type=str,help='bed path',required=True)
    parser.add_argument('--output','-o',type=str,help='output path',default='result.txt')
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.bam, "rb")

    with open(args.bed, "r") as bedfile:
       regions = [line.strip().split() for line in bedfile]

    classifier(regions,args.output)

    bamfile.close()

