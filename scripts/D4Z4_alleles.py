# D4Z4End2End: processing spanning reads and consensus calling for D4Z4 alleles
# usage: python3 D4Z4_alleles.py -i [haplotypes.txt file] -b [bam file] -o [output directory] -4q [number of 4q alleles] -10q [number of 10q alleles]

import os
import sys
import subprocess
import argparse
import numpy as np
from jenkspy import JenksNaturalBreaks

parser = argparse.ArgumentParser(description='Process some inputs')

parser.add_argument('-i', '--input', type=str, default='D4Z4_haplotyping/haplotypes.txt', help='haplotypes file (default: D4Z4_haplotyping/haplotypes.txt')
parser.add_argument('-b', '--bam', type=str, help='BAM file containing 4q and 10q reads')
parser.add_argument('-o', '--output', type=str, default='D4Z4_alleles', help='Output directory (default: D4Z4_alleles)')
parser.add_argument('-4q', '--clusters4q', type=int, help='number of 4q alleles')
parser.add_argument('-10q', '--clusters10q', type=int, help='number of 10q alleles')
parser.add_argument('--lists_only', action='store_true', help='only create lists of spanning read ids')
parser.add_argument('--no_consensus', action='store_true', help='do not create consensus sequences')
args = parser.parse_args()

haplotypes_file = os.path.join("output", args.input)
modBAM = args.bam
out_dir = os.path.join("output", args.output)
jnb_4q = args.clusters4q
jnb_10q = args.clusters10q

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

def count_reads(file_path):
    with open(file_path, 'r') as fq_file:
        total_lines = sum(1 for _ in fq_file)
    return total_lines // 4
    
def spanningreads_byhaplotype(haplotypesfile, chr, jnb_clusters):

    ## Assigning spanning reads to haplotypes based on number of D4Z4 units

    haplotypes_file = open(haplotypesfile, 'r')
    haplotypes_list = haplotypes_file.readlines()
    haplotypes_file.close()

    spanning = []

    for record in haplotypes_list[1:]:
        readinfo_list = record.rstrip().split('\t')
        if readinfo_list[3] == 'yes' and readinfo_list[4] == chr:
            spanning.append(readinfo_list)
        else:
            pass

    if jnb_clusters > 1:
        D4Z4count = np.asarray(spanning)[:,2]
        D4Z4count = np.asarray(D4Z4count, dtype=int)

        jnb = JenksNaturalBreaks(jnb_clusters) # number of 4q/10q alleles
        jnb.fit(D4Z4count) # create the clusters according to values in 'D4Z4count'

        labelled_spanning = list(zip(np.asarray(spanning),jnb.labels_))
    else:
        labelled_spanning = list(zip(np.asarray(spanning), np.zeros(len(spanning), dtype=int)))

    haplotypes = [[] for _ in range(max(jnb_clusters, 1))]

    for read in labelled_spanning:
        haplotypes[read[-1]].append(read[0].tolist())

    ## Processing spanning reads for each haplotype

    ref_files = []
    bam_files = []

    for i in range(jnb_clusters):
        
        if chr == 'chr10':
            chrom = '10q'
        if chr == 'chr4':
            chrom = '4q'

        hap_dir = os.path.join(out_dir, f'spanning_{chrom}_hap{i+1}')
        os.makedirs(hap_dir, exist_ok=True)

        # save lists of read ids for each haplotype

        with open(os.path.join(hap_dir, f'spanning_{chrom}_hap{i+1}.lst'), 'w') as output:
            for readid in haplotypes[i]:
                output.write(f'{readid[0]}\n')

        with open(os.path.join(hap_dir, f'spanning_{chrom}_hap{i+1}.txt'), 'w') as output:
            for readid in haplotypes[i]:
                output.write('\t'.join(readid) + '\n')
        
        # extract longest spanning read for each allele to use as the reference sequence
        
        if args.lists_only == False:
            os.system(f'samtools view -h -N {os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.lst")} {modBAM} | samtools fastq -T "*" - > {os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.fq")}')
            input_file = os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.fq")
            rawref = os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.rawref.fa")
            rawrest = os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.rawrest.fa")

            command = (
                f"seqtk seq -a {input_file} | "
                f"awk '/^>/ {{print $1}} !/^>/ {{print}}' | "
                f"seqkit sort -l -r | seqkit head -n 1 > {rawref} && "
                f"seqtk seq -a {input_file} | "
                f"awk '/^>/ {{print $1}} !/^>/ {{print}}' | "
                f"seqkit sort -l -r | seqkit range -r 2:-1 > {rawrest} "
            )
            try:
                subprocess.run(command, shell=True, check=True, executable="/bin/bash")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while processing haplotype {i+1}: {e}")
            
            with open(rawref, 'r') as file:
                read_id = file.readline().strip()[1:]
            with open(os.path.join(hap_dir, f'spanning_{chrom}_hap{i+1}.txt'), 'r') as file:
                lines = file.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith(read_id):
                    strand = line.split('\t')[6]
                    break
            # reverse complement if strand == '-'
            if strand == '-':
                os.system(f'seqtk seq -r {rawref} > {rawref}.revcomp')
                os.system(f'mv {rawref}.revcomp {rawref}')

            with open(rawref, 'r') as file:
                lines = file.readlines()
            lines[0] = f">{chrom}_hap{i+1} {lines[0][1:]}"
            with open(rawref, 'w') as file:
                file.writelines(lines)

            # create consensus sequence using Racon if there are >= 3 spanning reads

            if args.no_consensus == False:
                if count_reads(os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.fq")) >= 3:
                    os.system(f'minimap2 -t 8 -ax map-ont {rawref} {rawrest} > {os.path.join(hap_dir, f"overlaps.sam")}')
                    os.system(f'racon -t 8 --no-trimming {rawrest} {os.path.join(hap_dir, f"overlaps.sam")} {rawref} > {os.path.join(hap_dir, f"{chrom}_hap{i+1}.consensus.fa")}')
                    with open(os.path.join(hap_dir, f"{chrom}_hap{i+1}.consensus.fa"), 'r') as file:
                        lines = file.readlines()
                    lines[0] = f">{chrom}_hap{i+1}\n"
                    with open(os.path.join(hap_dir, f"{chrom}_hap{i+1}.consensus.fa"), 'w') as file:
                        file.writelines(lines)
                    ref_files.append(os.path.join(hap_dir, f"{chrom}_hap{i+1}.consensus.fa"))
                    os.system(f'rm {os.path.join(hap_dir, f"overlaps.sam")}')
                else:
                    ref_files.append(rawref)
            else:
                ref_files.append(rawref)

            os.system(f'rm {rawrest}')

            # re-map spanning reads for each allele against the corresponding reference sequence (either longest raw spanning read or consensus sequence)

            if os.path.exists(os.path.join(hap_dir, f"{chrom}_hap{i+1}.consensus.fa")):
                os.system(f'minimap2 -y -t 8 -ax map-ont {os.path.join(hap_dir, f"{chrom}_hap{i+1}.consensus.fa")} {os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.fq")} | samtools view -b -F 0x900 | samtools sort > {os.path.join(hap_dir, f"{chrom}_hap{i+1}.bam")}')
                bam_files.append(os.path.join(hap_dir, f"{chrom}_hap{i+1}.bam"))
            else:
                os.system(f'minimap2 -y -t 8 -ax map-ont {rawref} {os.path.join(hap_dir, f"spanning_{chrom}_hap{i+1}.fq")} | samtools view -b -F 0x900 | samtools sort > {os.path.join(hap_dir, f"{chrom}_hap{i+1}.bam")}')
                bam_files.append(os.path.join(hap_dir, f"{chrom}_hap{i+1}.bam"))  
    
    return ref_files, bam_files

files_4q = spanningreads_byhaplotype(haplotypes_file, 'chr4', jnb_4q)
files_10q = spanningreads_byhaplotype(haplotypes_file, 'chr10', jnb_10q)

## Write reference sequences for each allele to D4Z4_alleles.fa
## Merge bam files for each allele into D4Z4_alleles.bam

if args.lists_only == False:
    ref_files = files_4q[0] + files_10q[0]
    with open(os.path.join(out_dir, 'D4Z4_alleles.fa'), 'w') as outfile:
        for ref_file in ref_files:
            with open(ref_file, 'r') as infile:
                outfile.write(infile.read())
    bam_files = files_4q[1] + files_10q[1]
    bam_files_str = " ".join(bam_files)
    os.system(f'samtools merge {os.path.join(out_dir, "D4Z4_alleles.bam")} {bam_files_str}')
    os.system(f'samtools index {os.path.join(out_dir, "D4Z4_alleles.bam")}')
