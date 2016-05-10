#!/usr/bin/env python

from __future__ import print_function

import pysam
import argparse
import sys
import re

chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, 
          '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, 
          '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, 
          '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25, 
          'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 
          'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 
          'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 
          'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24, 'chrM':25}

def parse_info(info):
    res = {}
    for tag in info.split(';'):
        if '=' in tag:
            k, v = tag.split('=')
            if ',' in v:
                v = v.split(',')
            res[k] = v
        else:
            res[tag] = None
    return res

def process_args():
    parser = argparse.ArgumentParser(description="Validate detected structural variants in PacBio data.")
    parser.add_argument("--sv_vcf", help="A VCF file containing detected structural variants [stdin]")
    parser.add_argument("--use_chr", help="Prefix chromosome names with 'chr'", action="store_true")
    parser.add_argument("--outfile", help="The output VCF file")
    parser.add_argument("pacbio_bam", help="A BAM file of PacBio alignments for the same individual")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()
    
    if args.sv_vcf:
        args.sv_vcf = open(args.sv_vcf)
    else:
        args.sv_vcf = sys.stdin
    if args.outfile:
        args.outfile = open(args.outfile, 'w')
    else:
        args.outfile = sys.stdout

    args.pacbio_bam = pysam.AlignmentFile(args.pacbio_bam)

    for line in args.sv_vcf:
        if line[0] == '#':
            print(line.rstrip())
            continue
        line = line.strip().split('\t')
        f_chrom, f_pos = line[:2]
        if args.use_chr:
            f_chrom = 'chr' + f_chrom
        f_pos = int(f_pos)
        alt = []
        if not f_chrom in chroms:
            continue
        if ',' in line[4]:
            _alt = line[4].split(',')
            for a in _alt:
                alt.append(re.split(r'[\[\]]', a))
        else:
            alt.append(re.split(r'[\[\]]', line[4]))
        val_res = []
        cnt_res = []
        info_dict = parse_info(line[7])
        for i, a in enumerate(alt):
            t_chrom, t_pos = a[1].split(':')
            if not t_chrom in chroms:
                val_res.append('*')
                cnt_res.append('*')
                continue
            if args.use_chr:
                t_chrom = "chr" + t_chrom
            t_pos = int(t_pos)
            uncertainty = 10
            qgap = 0
            if 'QGAP' in info_dict:
                if not type(info_dict['QGAP']) is str:
                    qgap = int(info_dict['QGAP'][i])
                else:
                    qgap = int(info_dict['QGAP'])
            if qgap < 0:
                uncertainty += abs(qgap)
            read_dict = {}
            reads = args.pacbio_bam.fetch(f_chrom, f_pos - uncertainty, f_pos + uncertainty)
            n_reads = 0
            supporting_reads = 0
            for read in reads:
                r_pos, q_pos, clip, n_op = 0, 0, 0, len(read.cigartuples)
                if not read.query_name in read_dict:
                    read_dict[read.query_name] = []
                    n_reads += 1
                for j, tup in enumerate(read.cigartuples):
                    op, op_len = tup
                    if op == 0:
                        r_pos += op_len
                        q_pos += op_len
                    elif op == 1:
                        if qgap * 0.8 < op_len and qgap * 1.2 > op_len:
                            supporting_reads += 1
                        q_pos += op_len
                    elif op == 2:
                        if f_pos - uncertainty < int(read.reference_start) + r_pos and f_pos + uncertainty > int(read.reference_start) + r_pos:
                            if t_pos - uncertainty < int(read.reference_start) + r_pos + op_len and t_pos + uncertainty > int(read.reference_start) + r_pos + op_len:
                                supporting_reads += 1
                        elif t_pos - uncertainty < int(read.reference_start) + r_pos and t_pos + uncertainty > int(read.reference_start) + r_pos:
                            if f_pos - uncertainty < int(read.reference_start) + r_pos + op_len and f_pos + uncertainty > int(read.reference_start) + r_pos + op_len:
                                supporting_reads += 1
                        r_pos += op_len
                    elif op == 4 or op == 5:
                        if j == 0:
                            clip = 1
                        else:
                            clip = 2
                        if read.query_name in read_dict:
                            read_dict[read.query_name].append({"clip":clip, "ref":int(read.reference_start) + r_pos, "chrom":read.reference_name, "qpos":q_pos})
                        else:
                            read_dict[read.query_name] = [{"clip":clip, "ref":int(read.reference_start) + r_pos, "chrom":read.reference_name, "qpos":q_pos}]
                        q_pos += op_len
                    else:
                        print("Error: Unrecognized operation {} in read {}".format(op, read.query_name), file=sys.stderr)
                        sys.exit()
            if not supporting_reads: # no ins/del tags, look for split alignment
                reads = args.pacbio_bam.fetch(t_chrom, t_pos - uncertainty, t_pos + uncertainty)
                for read in reads:
                    r_pos, q_pos, clip, n_op = 0, 0, 0, len(read.cigartuples)
                    if not read.query_name in read_dict:
                        continue
                    for j, tup in enumerate(read.cigartuples):
                        op, op_len = tup
                        if op == 0:
                            q_pos += op_len
                            r_pos += op_len
                        elif op == 1:
                            q_pos += op_len
                        elif op == 2:
                            r_pos += op_len
                        elif op == 4 or op == 5:
                            if j == 0:
                                clip = 1
                            else:
                                clip = 2
                            if read.query_name in read_dict:
                                read_dict[read.query_name].append({"clip":clip, "ref":int(read.reference_start) + r_pos, "chrom":read.reference_name, "qpos":q_pos})
                            else:
                                read_dict[read.query_name] = [{"clip":clip, "ref":int(read.reference_start) + r_pos, "chrom":read.reference_name, "qpos":q_pos}]
                            q_pos += op_len
                for read, val in read_dict.items():
                    if len(val) > 2:
                        continue
                    found = False
                    for j, v in enumerate(val):
                        if f_chrom == v["chrom"] and f_pos - 50 < v["ref"] and f_pos + 50 > v["ref"]:
                            for v2 in val:
                                if t_chrom == v2["chrom"] and t_pos - 50 < v2["ref"] and t_pos + 50 > v2["ref"]:
                                    found = True
                                    break
                            if found:
                                break
                        elif t_chrom == v["chrom"] and t_pos - 50 < v["ref"] and t_pos + 50 > v["ref"]:
                            for v2 in val:
                                if f_chrom == v2["chrom"] and f_pos - 50 < v2["ref"] and f_pos + 50 > v2["ref"]:
                                    found = True
                                    break
                            if found:
                                break
                    if found:
                        supporting_reads += 1
            cnt_res.append("{}-{}".format(supporting_reads, n_reads))
            if supporting_reads > n_reads * 0.1:
                val_res.append('T')
            elif n_reads > 5:
                val_res.append('F')
            else:
                val_res.append('*')
        line[7] = line[7] + ";PACBIO_VALID={};PACBIO_CNT={}".format(','.join(val_res), ','.join(cnt_res))
        print('\t'.join(line), file=args.outfile)
    args.outfile.close()
    args.sv_vcf.close()

if __name__ == "__main__":
    main(None)
