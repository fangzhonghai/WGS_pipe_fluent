# -*- coding:utf-8 -*-
from functools import reduce
import pandas as pd
import subprocess
import optparse
import logging
import yaml
import glob
import time
import sys
import os
import re


class CreateJob:
    def __init__(self, config, path, sample, gender, fq_info_df):
        self.config = config
        self.path = path
        self.sample = sample
        self.gender = gender
        self.fq_info_df = fq_info_df
        fastp_split = self.config['parameters']['fastp_split_num']
        self.bam_split = ["{:0>4d}".format(i) for i in range(1, int(fastp_split) + 1)]

    def fq_filter(self):
        job_dic = {}
        filter_path = os.path.join(self.path, 'filter')
        if not os.path.exists(filter_path):
            os.makedirs(filter_path)
        filter_script_path = os.path.join(filter_path, 'script')
        if not os.path.exists(filter_script_path):
            os.makedirs(filter_script_path)
        fastp = self.config['software']['fastp']
        param = self.config['parameters']['fastp']
        resource = self.config['resource']['fastp']
        for i in range(self.fq_info_df.shape[0]):
            chip_lane_index = self.fq_info_df.loc[i, 'chip_lane_index']
            filter_fq_path = os.path.join(filter_path, chip_lane_index)
            if not os.path.exists(filter_fq_path):
                os.makedirs(filter_fq_path)
            fq1_path = self.fq_info_df.loc[i, 'fq1_path']
            fq2_path = self.fq_info_df.loc[i, 'fq2_path']
            fq1_name = self.fq_info_df.loc[i, 'fq1_name']
            fq2_name = self.fq_info_df.loc[i, 'fq2_name']
            script = os.path.join(filter_script_path, 'filter.' + chip_lane_index + '.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'filter'
            job_dic[script]['parent'] = []
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chip_lane_index} fq filter and split start `date`
{fastp} {param} \
-i {fq1_path} -I {fq2_path} \
-o {filter_fq_path}/fastp.{fq1_name} -O {filter_fq_path}/fastp.{fq2_name} \
-j {filter_fq_path}/{chip_lane_index}.fastp.json -h {filter_fq_path}/{chip_lane_index}.fastp.html
if [ $? -ne 0 ]; then
    echo {self.sample} {chip_lane_index} fq filter and split failed.
    exit 1
else
    echo {self.sample} {chip_lane_index} fq filter and split end `date`
    touch {script}.check
fi
'''.format(**locals())
                f.write(shell)
        return job_dic

    def fq_align(self):
        job_dic = {}
        align_path = os.path.join(self.path, 'align')
        if not os.path.exists(align_path):
            os.makedirs(align_path)
        align_script_path = os.path.join(align_path, 'script/alignment')
        if not os.path.exists(align_script_path):
            os.makedirs(align_script_path)
        bwa = self.config['software']['bwa']
        param = self.config['parameters']['bwa']
        resource = self.config['resource']['bwa']
        reference = self.config['database']['hg19']
        samtools = self.config['software']['samtools']
        for i in range(self.fq_info_df.shape[0]):
            chip_lane_index = self.fq_info_df.loc[i, 'chip_lane_index']
            align_bam_path = os.path.join(align_path, chip_lane_index)
            if not os.path.exists(align_bam_path):
                os.makedirs(align_bam_path)
            fq1_name = self.fq_info_df.loc[i, 'fq1_name']
            fq2_name = self.fq_info_df.loc[i, 'fq2_name']
            for n in self.bam_split:
                align_bam_split_path = os.path.join(align_bam_path, n)
                if not os.path.exists(align_bam_split_path):
                    os.makedirs(align_bam_split_path)
                script = os.path.join(align_script_path, 'align.' + n + "." + chip_lane_index + '.sh')
                script_flag = script + ".check"
                job_dic[script] = {}
                job_dic[script]['flow'] = 'align'
                job_dic[script]['parent'] = [self.path + '/filter/script/filter.' + chip_lane_index + '.sh']
                job_dic[script]['child'] = []
                job_dic[script]['resource'] = resource
                job_dic[script]['jobid'] = ''
                if os.path.exists(script_flag):
                    job_dic[script]['status'] = 'complete'
                else:
                    job_dic[script]['status'] = 'incomplete'
                with open(script, 'w') as f:
                    shell = r'''#!/bin/bash
echo {self.sample} {chip_lane_index} fq {n} part align start `date`
{bwa} {param} -R "@RG\tID:{self.sample}\tSM:{self.sample}\tPL:illumina\tLB:{self.sample}" {reference} \
{self.path}/filter/{chip_lane_index}/{n}.fastp.{fq1_name} \
{self.path}/filter/{chip_lane_index}/{n}.fastp.{fq2_name} | {samtools} view -S -b -o \
{align_bam_split_path}/{chip_lane_index}.{n}.align.bam -
if [ $? -ne 0]; then
    echo {self.sample} {chip_lane_index} fq {n} part align failed
    exit 1
else
    echo {self.sample} {chip_lane_index} fq {n} part align end `date`
    touch {script}.check
fi
'''.format(**locals())
                    f.write(shell)
        return job_dic

    def bam_sort(self):
        job_dic = {}
        align_path = os.path.join(self.path, 'align')
        sort_script_path = os.path.join(align_path, 'script/sort')
        if not os.path.exists(sort_script_path):
            os.makedirs(sort_script_path)
        sambamba = self.config['software']['sambamba']
        param = self.config['parameters']['sambamba_sort']
        resource = self.config['resource']['sambamba_sort']
        samtools = self.config['software']['samtools']
        for i in range(self.fq_info_df.shape[0]):
            chip_lane_index = self.fq_info_df.loc[i, 'chip_lane_index']
            align_bam_path = os.path.join(align_path, chip_lane_index)
            for n in self.bam_split:
                align_bam_split_path = os.path.join(align_bam_path, n)
                script = os.path.join(sort_script_path, 'sort.' + n + "." + chip_lane_index + '.sh')
                script_flag = script + ".check"
                job_dic[script] = {}
                job_dic[script]['flow'] = 'sort'
                job_dic[script]['parent'] = [self.path + '/align/script/alignment/align.' + n + '.' + chip_lane_index + '.sh']
                job_dic[script]['child'] = []
                job_dic[script]['resource'] = resource
                job_dic[script]['jobid'] = ''
                if os.path.exists(script_flag):
                    job_dic[script]['status'] = 'complete'
                else:
                    job_dic[script]['status'] = 'incomplete'
                with open(script, 'w') as f:
                    shell = r'''#!/bin/bash
echo {self.sample} {chip_lane_index} bam {n} part sort start `date`
{sambamba} {param} --tmpdir=TMPDIR -o {align_bam_split_path}/{chip_lane_index}.{n}.align.sort.bam \
{align_bam_split_path}/{chip_lane_index}.{n}.align.bam
if [ $? -ne 0]; then
    echo {self.sample} {chip_lane_index} bam {n} part sort failed
    exit 1
else
    {samtools} index {align_bam_split_path}/{chip_lane_index}.{n}.align.sort.bam
    echo {self.sample} {chip_lane_index} bam {n} part sort end `date`
    touch {script}.check
fi
'''.format(**locals())
                    f.write(shell)
        return job_dic

    def bam_merge(self):
        job_dic = {}
        align_path = os.path.join(self.path, 'align')
        sort_script_path = os.path.join(align_path, 'script/sort')
        bam_path = os.path.join(self.path, 'bam_chr')
        merge_script_path = os.path.join(bam_path, 'script/merge')
        if not os.path.exists(merge_script_path):
            os.makedirs(merge_script_path)
        samtools = self.config['software']['samtools']
        param = self.config['parameters']['samtools_merge']
        resource = self.config['resource']['samtools_merge']
        chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM_NC_012920.1']
        sort_scripts = []
        sort_bams = os.path.join(merge_script_path, 'all.sort.bam.list')
        with open(sort_bams, 'w') as f:
            for i in range(self.fq_info_df.shape[0]):
                chip_lane_index = self.fq_info_df.loc[i, 'chip_lane_index']
                align_bam_path = os.path.join(align_path, chip_lane_index)
                for n in self.bam_split:
                    sort_scripts.append(os.path.join(sort_script_path, 'sort.' + n + "." + chip_lane_index + '.sh'))
                    align_bam_split_path = os.path.join(align_bam_path, n)
                    shell = '''{align_bam_split_path}/{chip_lane_index}.{n}.align.sort.bam
'''.format(**locals())
                    f.write(shell)
        for chrom in chroms:
            script = os.path.join(merge_script_path, chrom + '.merge.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'merge'
            job_dic[script]['parent'] = sort_scripts
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chrom} bam sort merge start `date`
{samtools} {param} -b {sort_bams} -R {chrom} {bam_path}/{chrom}.sort.bam
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} bam sort merge failed
    exit 1
else
    echo {self.sample} {chrom} bam sort merge end `date`
    touch {script}.check
fi
'''.format(**locals())
                f.write(shell)
        return job_dic

    def dupmark(self):
        job_dic = {}
        bam_path = os.path.join(self.path, 'bam_chr')
        merge_script_path = os.path.join(bam_path, 'script/merge')
        dup_script_path = os.path.join(bam_path, 'script/dupmark')
        if not os.path.exists(dup_script_path):
            os.makedirs(dup_script_path)
        gatk = self.config['software']['gatk']
        param = self.config['parameters']['dupmark']
        resource = self.config['resource']['dupmark']
        chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM_NC_012920.1']
        for chrom in chroms:
            tmp_dir = os.path.join(bam_path, 'script/tmp/' + chrom)
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)
            script = os.path.join(dup_script_path, chrom + '.dupmark.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'dupmark'
            job_dic[script]['parent'] = [os.path.join(merge_script_path, chrom + '.merge.sh')]
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chrom} bam dupmark start `date`
{gatk} MarkDuplicates {param} --TMP_DIR {tmp_dir} \
-I {bam_path}/{chrom}.sort.bam -O {bam_path}/{chrom}.sort.dup.bam -M {bam_path}/{chrom}.sort.dup.bam.metrics
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} bam dupmark failed
    exit 1
else
    echo {self.sample} {chrom} bam dupmark end `date`
    touch {script}.check
fi
'''.format(**locals())
                f.write(shell)
        return job_dic

    def fixmate(self):
        job_dic = {}
        bam_path = os.path.join(self.path, 'bam_chr')
        dup_script_path = os.path.join(bam_path, 'script/dupmark')
        fix_script_path = os.path.join(bam_path, 'script/fixmate')
        if not os.path.exists(fix_script_path):
            os.makedirs(fix_script_path)
        gatk = self.config['software']['gatk']
        samtools = self.config['software']['samtools']
        param = self.config['parameters']['fixmate']
        resource = self.config['resource']['fixmate']
        chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM_NC_012920.1']
        for chrom in chroms:
            tmp_dir = os.path.join(bam_path, 'script/tmp/' + chrom)
            script = os.path.join(fix_script_path, chrom + '.fixmate.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'fixmate'
            job_dic[script]['parent'] = [os.path.join(dup_script_path, chrom + '.dupmark.sh')]
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chrom} bam fixmate start `date`
{gatk} FixMateInformation {param} --TMP_DIR {tmp_dir} \
-I {bam_path}/{chrom}.sort.dup.bam -O {bam_path}/{chrom}.sort.dup.fix.bam
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} bam fixmate failed
    exit 1
else
    {samtools} index {bam_path}/{chrom}.sort.dup.fix.bam
    echo {self.sample} {chrom} bam fixmate end `date`
    touch {script}.check
fi
'''.format(**locals())
                f.write(shell)
        return job_dic

    def bqsr(self):
        job_dic = {}
        bam_path = os.path.join(self.path, 'bam_chr')
        fix_script_path = os.path.join(bam_path, 'script/fixmate')
        bqsr_script_path = os.path.join(bam_path, 'script/bqsr')
        if not os.path.exists(bqsr_script_path):
            os.makedirs(bqsr_script_path)
        gatk = self.config['software']['gatk']
        reference = self.config['database']['hg19']
        samtools = self.config['software']['samtools']
        param = self.config['parameters']['bqsr']
        resource = self.config['resource']['bqsr']
        chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM_NC_012920.1']
        for chrom in chroms:
            tmp_dir = os.path.join(bam_path, 'script/tmp/' + chrom)
            script = os.path.join(bqsr_script_path, chrom + '.bqsr.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'bqsr'
            job_dic[script]['parent'] = [os.path.join(fix_script_path, chrom + '.fixmate.sh')]
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chrom} bam bqsr start `date`
{gatk} BaseRecalibrator --tmp-dir {tmp_dir} \
-I {bam_path}/{chrom}.sort.dup.fix.bam -O {bam_path}/{chrom}.bqsr.table -R {reference} \
{param}
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} bam BaseRecalibrator failed
    exit 1
fi
{gatk} ApplyBQSR --tmp-dir {tmp_dir} \
-I {bam_path}/{chrom}.sort.dup.fix.bam -bqsr {bam_path}/{chrom}.bqsr.table -O {bam_path}/{chrom}.bqsr.bam 
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} bam ApplyBQSR failed
    exit 1
else
    {samtools} index {bam_path}/{chrom}.bqsr.bam
    echo {self.sample} {chrom} bam bqsr end `date`
    touch {script}.check
fi
'''.format(**locals())
                f.write(shell)
        return job_dic

    def bam_part_split(self):
        job_dic = {}
        bam_path = os.path.join(self.path, 'bam_chr')
        bam_split_path = os.path.join(self.path, 'QC/bam_chr_part_split')
        bam_split_script_path = os.path.join(bam_split_path, 'script')
        if not os.path.exists(bam_split_script_path):
            os.makedirs(bam_split_script_path)
        sambamba = self.config['software']['sambamba']
        resource = self.config['resource']['sambamba_view']
        beds = glob.glob(self.config['database']['genome_bed'])
        for b in beds:
            bed = os.path.basename(b)
            bed_name_split = bed.split('.')
            if not bed.startswith('chrM'):
                chrom = str(bed_name_split[0])
                part = str(bed_name_split[1])
            else:
                chrom = str(bed_name_split[0]) + '.' + str(bed_name_split[1])
                part = str(bed_name_split[2])
            bam_split_part_path = os.path.join(bam_split_path, chrom + '/' + part)
            if not os.path.exists(bam_split_part_path):
                os.makedirs(bam_split_part_path)
            script = os.path.join(bam_split_script_path, chrom + '.' + part + '.split.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'split'
            job_dic[script]['parent'] = [os.path.join(bam_path, 'script/bqsr/' + chrom + '.bqsr.sh')]
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chrom} {part} bam split start `date`
{sambamba} view -f bam -h -o {bam_split_part_path}/{chrom}_{part}.bam -L {b} {bam_path}/{chrom}.bqsr.bam
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} {part} bam split failed
    exit 1
else
    {sambamba} index {bam_split_part_path}/{chrom}_{part}.bam
    echo {self.sample} {chrom} {part} bam split end `date`
    touch {script}.check
fi
'''.format(**locals())
                f.write(shell)
        return job_dic

    def qc(self):
        job_dic = {}
        bam_split_path = os.path.join(self.path, 'QC/bam_chr_part_split')
        bam_split_qc_path = os.path.join(self.path, 'QC/QC_BED_SPLIT')
        bam_split_qc_script_path = os.path.join(bam_split_qc_path, 'script')
        if not os.path.exists(bam_split_qc_script_path):
            os.makedirs(bam_split_qc_script_path)
        bamdst = self.config['software']['bamdst']
        resource = self.config['resource']['bamdst']
        beds = glob.glob(self.config['database']['genome_bed'])
        for b in beds:
            bed = os.path.basename(b)
            bed_name_split = bed.split('.')
            if not bed.startswith('chrM'):
                chrom = str(bed_name_split[0])
                part = str(bed_name_split[1])
            else:
                chrom = str(bed_name_split[0]) + '.' + str(bed_name_split[1])
                part = str(bed_name_split[2])
            bam_split_part_path = os.path.join(bam_split_path, chrom + '/' + part)
            bam_split_qc_part_path = os.path.join(bam_split_qc_path, chrom + '/' + part)
            if not os.path.exists(bam_split_qc_part_path):
                os.makedirs(bam_split_qc_part_path)
            script = os.path.join(bam_split_qc_script_path, chrom + '.' + part + '.qc.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'qc'
            job_dic[script]['parent'] = [os.path.join(bam_split_path, 'script/' + chrom + '.' + part + '.split.sh')]
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chrom} {part} bam qc start `date`
{bamdst} -p {b} -o {bam_split_qc_part_path} {bam_split_part_path}/{chrom}_{part}.bam
echo {self.sample} {chrom} {part} bam qc end `date`
touch {script}.check
'''.format(**locals())
                f.write(shell)
        return job_dic

    def haplotypecaller(self):
        job_dic = {}
        bam_path = os.path.join(self.path, 'bam_chr')
        variant_path = os.path.join(self.path, 'variant')
        variant_script_path = os.path.join(variant_path, 'script')
        tmp_dir = os.path.join(variant_script_path, 'tmp')
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        gatk = self.config['software']['gatk']
        resource = self.config['resource']['HaplotypeCaller']
        param = self.config['parameters']['HaplotypeCaller']
        beds = glob.glob(self.config['database']['genotypeGVCFs_bed'])
        reference = self.config['database']['hg19']
        for b in beds:
            bed = os.path.basename(b)
            bed_name_split = bed.split('.')
            if not bed.startswith('chrM'):
                chrom = str(bed_name_split[0])
                part = str(bed_name_split[1])
            else:
                chrom = str(bed_name_split[0]) + '.' + str(bed_name_split[1])
                part = str(bed_name_split[2])
            variant_part_path = os.path.join(variant_path, chrom + '/' + part)
            if not os.path.exists(variant_part_path):
                os.makedirs(variant_part_path)
            script = os.path.join(variant_script_path, chrom + '.' + part + '.hc.sh')
            script_flag = script + ".check"
            job_dic[script] = {}
            job_dic[script]['flow'] = 'haplotypecaller'
            job_dic[script]['parent'] = [os.path.join(bam_path, 'script/bqsr/' + chrom + '.bqsr.sh')]
            job_dic[script]['child'] = []
            job_dic[script]['resource'] = resource
            job_dic[script]['jobid'] = ''
            if os.path.exists(script_flag):
                job_dic[script]['status'] = 'complete'
            else:
                job_dic[script]['status'] = 'incomplete'
            with open(script, 'w') as f:
                shell = r'''#!/bin/bash
echo {self.sample} {chrom} {part} haplotypecaller start `date`
{gatk} HaplotypeCaller --tmp-dir {tmp_dir} -ERC GVCF \
-L {b} \
-R {reference} \
-I {bam_path}/{chrom}.bqsr.bam \
-O {variant_part_path}/{chrom}_{part}.gvcf.gz \
{param}
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} {part} haplotypecaller failed
    exit 1
fi
{gatk} GenotypeGVCFs --tmp-dir {tmp_dir} \
-R {reference} \
-V {variant_part_path}/{chrom}_{part}.gvcf.gz \
-L {b} \
-O {variant_part_path}/{chrom}_{part}.vcf.gz
if [ $? -ne 0 ]; then
    echo {self.sample} {chrom} {part} GenotypeGVCFs failed
    exit 1
else
    echo {self.sample} {chrom} {part} GenotypeGVCFs end `date`
    touch {script}.check
fi
'''.format(**locals())
                f.write(shell)
        return job_dic

    def concat_gvcf(self):
        job_dic = {}
        variant_path = os.path.join(self.path, 'variant')
        variant_script_path = os.path.join(variant_path, 'script')
        bcftools = self.config['software']['bcftools']
        resource = self.config['resource']['bcftools']
        param = self.config['parameters']['bcftools_concat']
        tabix = self.config['software']['tabix']
        beds = glob.glob(self.config['database']['genotypeGVCFs_bed'])
        variant_hc_scripts = []
        gvcfs = []
        vcfs = []
        for b in beds:
            bed = os.path.basename(b)
            bed_name_split = bed.split('.')
            if not bed.startswith('chrM'):
                chrom = str(bed_name_split[0])
                part = str(bed_name_split[1])
            else:
                chrom = str(bed_name_split[0]) + '.' + str(bed_name_split[1])
                part = str(bed_name_split[2])
            variant_part_path = os.path.join(variant_path, chrom + '/' + part)
            variant_hc_scripts.append(os.path.join(variant_script_path, chrom + '.' + part + '.hc.sh'))
            gvcfs.append(os.path.join(variant_part_path, chrom + '_' + part + '.gvcf.gz'))
            vcfs.append(os.path.join(variant_part_path, chrom + '_' + part + '.vcf.gz'))
        gvcfs_file = pd.DataFrame({'gvcf': gvcfs}).to_csv(variant_path + '/part_gvcf.list', sep='\t', index=False, header=None)
        vcfs_file = pd.DataFrame({'vcf': vcfs}).to_csv(variant_path + '/part_vcf.list', sep='\t', index=False, header=None)
        script = os.path.join(variant_script_path, 'concat.part_vcf.' + self.sample + '.sh')
        script_flag = script + ".check"
        job_dic[script] = {}
        job_dic[script]['flow'] = 'concat_vcfs'
        job_dic[script]['parent'] = variant_hc_scripts
        job_dic[script]['child'] = []
        job_dic[script]['resource'] = resource
        job_dic[script]['jobid'] = ''
        if os.path.exists(script_flag):
            job_dic[script]['status'] = 'complete'
        else:
            job_dic[script]['status'] = 'incomplete'
        with open(script, 'w') as f:
            shell = r'''#!/bin/bash
echo {self.sample} concat vcf start `date`
{bcftools} {param} -o {variant_path}/{self.sample}.wgs.vcf.gz -f {variant_path}/part_vcf.list
if [ $? -ne 0 ]; then
    echo {self.sample} concat vcf failed
    exit 1
fi
{tabix} -p vcf {variant_path}/{self.sample}.wgs.vcf.gz
if [ $? -ne 0 ]; then
    echo {self.sample} tabix vcf failed
    exit 1
fi
{bcftools} {param} -o {variant_path}/{self.sample}.wgs.gvcf.gz -f {variant_path}/part_gvcf.list
if [ $? -ne 0 ]; then
    echo {self.sample} concat gvcf failed
    exit 1
fi
{tabix} -p vcf {variant_path}/{self.sample}.wgs.gvcf.gz
if [ $? -ne 0 ]; then
    echo {self.sample} tabix gvcf failed
    exit 1
fi
echo {self.sample} concat vcf end `date`
touch {script}.check
'''.format(**locals())
            f.write(shell)
        return job_dic

    def vqsr_snp(self):
        job_dic = {}
        variant_path = os.path.join(self.path, 'variant')
        variant_script_path = os.path.join(variant_path, 'script')
        tmp_dir = os.path.join(variant_path, 'script/snp_tmp')
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        gatk = self.config['software']['gatk']
        resource = self.config['resource']['vqsr']
        param = self.config['parameters']['vqsr_snp']
        reference = self.config['database']['hg19']
        script = os.path.join(variant_script_path, 'vqsr.snp.' + self.sample + '.sh')
        script_flag = script + ".check"
        job_dic[script] = {}
        job_dic[script]['flow'] = 'vqsr_snp'
        job_dic[script]['parent'] = [os.path.join(variant_script_path, 'concat.part_vcf.' + self.sample + '.sh')]
        job_dic[script]['child'] = []
        job_dic[script]['resource'] = resource
        job_dic[script]['jobid'] = ''
        if os.path.exists(script_flag):
            job_dic[script]['status'] = 'complete'
        else:
            job_dic[script]['status'] = 'incomplete'
        with open(script, 'w') as f:
            shell = r'''#!/bin/bash
echo {self.sample} vqsr snp start `date`
{gatk} SelectVariants -select-type SNP --tmp-dir {tmp_dir} \
-R {reference} \
--variant {variant_path}/{self.sample}.wgs.vcf.gz \
-O {variant_path}/{self.sample}.wgs.snp.vcf
if [ $? -ne 0 ]; then
    echo {self.sample} select snp failed
    exit 1
fi
{gatk} VariantRecalibrator -mode SNP --max-gaussians 4 \
--tmp-dir {tmp_dir} \
-R {reference} \
-V {variant_path}/{self.sample}.wgs.snp.vcf \
-O {variant_path}/{self.sample}.wgs.snp_recal.vcf \
--tranches-file {variant_path}/{self.sample}.wgs.snp_tranches.vcf \
{param}
if [ $? -ne 0 ]; then
    echo {self.sample} snp VariantRecalibrator failed
    exit 1
fi
{gatk} ApplyVQSR --tmp-dir {tmp_dir} \
-R {reference} \
-V {variant_path}/{self.sample}.wgs.snp.vcf \
-O {variant_path}/{self.sample}.wgs.snp.filter.vcf \
--truth-sensitivity-filter-level 99.0 -mode SNP \
--tranches-file {variant_path}/{self.sample}.wgs.snp_tranches.vcf \
--recal-file {variant_path}/{self.sample}.wgs.snp_recal.vcf
if [ $? -ne 0 ]; then
    echo {self.sample} snp ApplyVQSR failed
    exit 1
fi
{gatk} SortVcf --TMP_DIR {tmp_dir} \
-I {variant_path}/{self.sample}.wgs.snp.filter.vcf \
-O {variant_path}/{self.sample}.wgs.snp.filter.sort.vcf
if [ $? -ne 0 ]; then
    echo {self.sample} snp sort vcf failed
    exit 1
fi
echo {self.sample} vqsr snp end `date`
touch {script}.check
'''.format(**locals())
            f.write(shell)
        return job_dic

    def vqsr_indel(self):
        job_dic = {}
        variant_path = os.path.join(self.path, 'variant')
        variant_script_path = os.path.join(variant_path, 'script')
        tmp_dir = os.path.join(variant_path, 'script/indel_tmp')
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        gatk = self.config['software']['gatk']
        resource = self.config['resource']['vqsr']
        param = self.config['parameters']['vqsr_indel']
        reference = self.config['database']['hg19']
        script = os.path.join(variant_script_path, 'vqsr.indel.' + self.sample + '.sh')
        script_flag = script + ".check"
        job_dic[script] = {}
        job_dic[script]['flow'] = 'vqsr_indel'
        job_dic[script]['parent'] = [os.path.join(variant_script_path, 'concat.part_vcf.' + self.sample + '.sh')]
        job_dic[script]['child'] = []
        job_dic[script]['resource'] = resource
        job_dic[script]['jobid'] = ''
        if os.path.exists(script_flag):
            job_dic[script]['status'] = 'complete'
        else:
            job_dic[script]['status'] = 'incomplete'
        with open(script, 'w') as f:
            shell = r'''#!/bin/bash
echo {self.sample} vqsr indel start `date`
{gatk} SelectVariants -select-type INDEL --tmp-dir {tmp_dir} \
-R {reference} \
--variant {variant_path}/{self.sample}.wgs.vcf.gz \
-O {variant_path}/{self.sample}.wgs.indel.vcf
if [ $? -ne 0 ]; then
    echo {self.sample} select indel failed
    exit 1
fi
{gatk} VariantRecalibrator -mode INDEL --max-gaussians 4 \
--tmp-dir {tmp_dir} \
-R {reference} \
-V {variant_path}/{self.sample}.wgs.indel.vcf \
-O {variant_path}/{self.sample}.wgs.indel_recal.vcf \
--tranches-file {variant_path}/{self.sample}.wgs.indel_tranches.vcf \
{param}
if [ $? -ne 0 ]; then
    echo {self.sample} indel VariantRecalibrator failed
    exit 1
fi
{gatk} ApplyVQSR --tmp-dir {tmp_dir} \
-R {reference} \
-V {variant_path}/{self.sample}.wgs.indel.vcf \
-O {variant_path}/{self.sample}.wgs.indel.filter.vcf \
--truth-sensitivity-filter-level 99.0 -mode INDEL \
--tranches-file {variant_path}/{self.sample}.wgs.indel_tranches.vcf \
--recal-file {variant_path}/{self.sample}.wgs.indel_recal.vcf
if [ $? -ne 0 ]; then
    echo {self.sample} indel ApplyVQSR failed
    exit 1
fi
{gatk} SortVcf --TMP_DIR {tmp_dir} \
-I {variant_path}/{self.sample}.wgs.indel.filter.vcf \
-O {variant_path}/{self.sample}.wgs.indel.filter.sort.vcf
if [ $? -ne 0 ]; then
    echo {self.sample} indel sort vcf failed
    exit 1
fi
echo {self.sample} vqsr indel end `date`
touch {script}.check
'''.format(**locals())
            f.write(shell)
        return job_dic

    def merge_vqsr(self):
        job_dic = {}
        variant_path = os.path.join(self.path, 'variant')
        variant_script_path = os.path.join(variant_path, 'script')
        tmp_dir = os.path.join(variant_script_path, 'tmp')
        gatk = self.config['software']['gatk']
        resource = self.config['resource']['merge_vcf']
        tabix = self.config['software']['tabix']
        script = os.path.join(variant_script_path, 'vqsr.merge.' + self.sample + '.sh')
        script_flag = script + ".check"
        job_dic[script] = {}
        job_dic[script]['flow'] = 'merge_vqsr'
        job_dic[script]['parent'] = [os.path.join(variant_script_path, 'vqsr.snp.' + self.sample + '.sh'),
                                     os.path.join(variant_script_path, 'vqsr.indel.' + self.sample + '.sh')]
        job_dic[script]['child'] = []
        job_dic[script]['resource'] = resource
        job_dic[script]['jobid'] = ''
        if os.path.exists(script_flag):
            job_dic[script]['status'] = 'complete'
        else:
            job_dic[script]['status'] = 'incomplete'
        with open(script, 'w') as f:
            shell = r'''#!/bin/bash
echo {self.sample} merge vqsr vcf start `date`
{gatk} MergeVcfs --TMP_DIR {tmp_dir}
-I {variant_path}/{self.sample}.wgs.snp.filter.sort.vcf \
-I {variant_path}/{self.sample}.wgs.indel.filter.sort.vcf \
-O {variant_path}/{self.sample}.wgs.filter.vcf.gz
if [ $? -ne 0 ]; then
    echo {self.sample} merge vqsr vcf failed
    exit 1
fi
{tabix} -p vcf {variant_path}/{self.sample}.wgs.filter.vcf.gz
if [ $? -ne 0 ]; then
    echo {self.sample} tabix vcf failed
    exit 1
fi
echo {self.sample} merge vqsr vcf end `date`
touch {script}.check
'''.format(**locals())
            f.write(shell)
        return job_dic


class RunJob:
    def __init__(self, job, job_dic, pro):
        self.job = job
        self.job_dic = job_dic
        self.pro = pro

    @property
    def parent_status(self):
        if len(self.job_dic['parent']) == 0:
            return 'finished'
        status_list = [parent['status'] for parent in self.job_dic['parent']]
        if 'incomplete' not in status_list:
            return 'finished'
        else:
            return 'unfinished'

    @property
    def status_in_sge(self):
        command = "qstat | grep " + "\"" + self.job_dic['jobid'] + " " + "\""
        status, output = subprocess.getstatusoutput(command)
        return status, output

    @staticmethod
    def job_num_in_sge():
        command = "qstat | grep `whoami` |wc -l"
        status, output = subprocess.getstatusoutput(command)
        return status, int(output)

    @staticmethod
    def job_id_in_sge(command):
        status, output = subprocess.getstatusoutput(command)
        jobid = re.findall(r"Your job (\d+) ", output)[0]
        return status, jobid

    @property
    def submit(self):
        job_num_status, job_num = RunJob.job_num_in_sge()
        if job_num_status != 0:
            sys.exit(1)
        while job_num >= 1999:
            time.sleep(600)
            job_num_status, job_num = RunJob.job_num_in_sge()
            if job_num_status != 0:
                sys.exit(1)
        job_path = os.path.dirname(self.job)
        command = "qsub -wd " + job_path + " -P " + self.pro + " " + self.job_dic['resource'] + " " + self.job
        status, jobid = RunJob.job_id_in_sge(command)
        return status, jobid


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def get_sample_info(fq_ls):
    df = pd.read_csv(fq_ls, sep='\t', header=None)
    df.columns = ['sample', 'fq1_path', 'gender']
    if len(df['sample'].unique()) > 1 or len(df['gender'].unique()) > 1 or len(df['fq1_path'].unique()) < df.shape[0]:
        print("Fatal: Sample Info ERROR!")
        sys.exit(1)
    df['fq2_path'] = df['fq1_path'].str.replace('_1.fq.gz', '_2.fq.gz')
    df['chip_lane_index'] = df['fq1_path'].map(os.path.basename).str.replace('_1.fq.gz', '')
    df['fq1_name'] = df['fq1_path'].map(os.path.basename)
    df['fq2_name'] = df['fq2_path'].map(os.path.basename)
    sample = df['sample'].values[0]
    gender = df['gender'].values[0]
    return sample, gender, df[['fq1_path', 'fq2_path', 'fq1_name', 'fq2_name', 'chip_lane_index']]


def jobs_map(jobs_dic) -> dict:
    jobs = list(filter(lambda x: jobs_dic[x]['flow'] != 'filter', jobs_dic.keys()))
    for i in jobs:
        for j in jobs_dic[i]['parent']:
            jobs_dic[j]['child'].append(i)
    return jobs_dic


def work_flow(jobs_dic, pro):
    # jobs_status = ['r', 'hr', 'qw', 'Eqw', 'hqw', 't', 's']
    queue = list(filter(lambda x: jobs_dic[x]['flow'] == 'filter', jobs_dic.keys()))
    while len(queue) > 0:
        queue_add = []
        queue_remove = []
        for q in queue:
            if os.path.exists(q + '.check'):
                jobs_dic[q]['status'] = 'complete'
            if jobs_dic[q]['status'] == 'incomplete':
                if jobs_dic[q]['jobid'] == '':
                    if RunJob(q, jobs_dic[q], pro).parent_status == 'finished':
                        command_status, command_jobid = RunJob(q, jobs_dic[q], pro).submit
                        if command_status != 0:
                            logger_main.error(q + ' Submit Failed!')
                            sys.exit(1)
                        jobs_dic[q]['jobid'] = command_jobid
                        logger_main.info(q + ' Submit Success!')
                else:
                    command_status, command_output = RunJob(q, jobs_dic[q], pro).status_in_sge
                    if command_status != 0 and command_output == '' and not os.path.exists(q + '.check'):
                        logger_main.error(q + ' Run Failed!')
                        sys.exit(1)
            else:
                queue_remove.append(q)
                queue_add += jobs_dic[q]['child']
                logger_main.info(q + ' Finished!')
        if len(queue_add) == 0:
            time.sleep(60)
        queue = list(set(list(set(queue) - set(queue_remove)) + queue_add))
    logger_main.info('All Jobs Finished!')


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-c', '--config', dest='config', help='config file', default=None, metavar='file')
    parser.add_option('-s', '--sample', dest='sample_info', help='sample fq info', default=None, metavar='file')
    parser.add_option('-p', '--project', dest='project', help='project queue in sge', default=None, metavar='string')
    parser.add_option('-o', '--outdir', dest='outdir', help='wgs project dir', default=None, metavar='string')
    parser.add_option('--onlyshell', dest='onlyshell', help='create scripts only', action='store_true')
    (opts, args) = parser.parse_args()
    config_yaml = opts.config
    sample_info = opts.sample_info
    project = opts.project
    outdir = opts.outdir
    onlyshell = opts.onlyshell
    config_dic = yaml_read(config_yaml)
    specimen, sex, fq_info = get_sample_info(sample_info)
    wkdir = os.path.join(outdir, specimen)
    if not os.path.exists(wkdir):
        os.makedirs(wkdir)
    os.system('cp ' + config_yaml + " " + wkdir)
    fh = logging.FileHandler(wkdir + '/main.log', 'w', encoding='utf-8')
    logger_main = logging.getLogger()
    logger_main.setLevel(logging.INFO)
    fm = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger_main.addHandler(fh)
    fh.setFormatter(fm)
    logger_main.info('WGS Pipeline Create Jobs Start!')
    filter_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).fq_filter()
    align_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).fq_align()
    sort_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).bam_sort()
    merge_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).bam_merge()
    dupmark_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).dupmark()
    fixmate_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).fixmate()
    bqsr_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).bqsr()
    split_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).bam_part_split()
    qc_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).qc()
    hc_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).haplotypecaller()
    concat_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).concat_gvcf()
    vqsrsnp_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).vqsr_snp()
    vqsrindel_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).vqsr_indel()
    vqsrmerge_job_dic = CreateJob(config_dic, wkdir, specimen, sex, fq_info).merge_vqsr()
    if onlyshell:
        logger_main.info('WGS Pipeline Create Jobs Finish and Exit!')
        sys.exit(0)
    nodes = [filter_job_dic, align_job_dic, sort_job_dic, merge_job_dic, dupmark_job_dic, fixmate_job_dic, bqsr_job_dic, split_job_dic,
             qc_job_dic, hc_job_dic, concat_job_dic, vqsrsnp_job_dic, vqsrindel_job_dic, vqsrmerge_job_dic]
    nodes_dic = reduce(lambda x, y: dict(x, **y), nodes)
    nodes_dic = jobs_map(nodes_dic)
    with open(wkdir + '/jobs_map.yaml', 'w') as fp:
        yaml.dump(nodes_dic, fp, default_flow_style=False)
        logger_main.info('WGS Pipeline Create Jobs & Map Finish!')
    work_flow(nodes_dic, project)
