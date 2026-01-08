# -*- coding: utf-8 -*-
"""
Fix on Thu Dec 18 14:41:21 CST 2025
telomere expand pipeline V9

@author: zhengshang@frasergen.com zhengshang-zn@qq.com

"""
import re
import sys
import os
data_from=sys.argv[1].strip()
work_dir=os.getcwd()
##conf reading
f=open(data_from,'r')
data=f.readlines()
f.close()
data=[i[:-1] for i in data if i[0]!="#" and i!="\n"]
s_a=[i.split('=')[0].strip() for i in data]
s_b=[i.split('=')[1].strip() for i in data]
data_dict=dict(zip(s_a,s_b))

script_path=os.path.dirname(os.path.realpath(__file__)) + "/software.txt"
f=open(script_path,'r')
data=f.readlines()
f.close()
data=[i[:-1] for i in data if i[0]!="#" and i!="\n"]
s_a=[i.split('=')[0].strip() for i in data]
s_b=[i.split('=')[1].strip() for i in data]
script_dict=dict(zip(s_a,s_b))

os.chdir(work_dir)
os.system('mkdir 00.gapcloser_merge  01.check_without_telomere  02.minimap  03.expand')
os.chdir('./00.gapcloser_merge')
f=open('00.gapcloser_merge.sh','w')
content="""#!/bin/bash
echo begin  at `date`
seqkit={0}
ln -s {1} gap_free.fa &&
ln -s {2} all.asm.fa &&
$seqkit seq all.asm.fa -w 10000 | grep -A 2 ">" | grep -i {3}{3}{3}{3}""".format(script_dict['seqkit'],data_dict['genome'],data_dict['asm'],data_dict['telomere_monomer'])+""" -B 1 |grep ">" | sed 's/>//g' |awk '{print $1}' > head_cis.id"""+"""
$seqkit seq all.asm.fa -r -p -w 10000 | grep -A 2 ">" | grep -i {0}{0}{0}{0}""".format(data_dict['telomere_monomer'])+""" -B 1 |grep ">" | sed 's/>//g' |awk '{print $1}' > tail_trans.id
$seqkit grep -f head_cis.id all.asm.fa   | $seqkit seq -w 0 | sed 's/>/>cis_/g'>candidate_contig.fa &&
$seqkit grep -f tail_trans.id all.asm.fa | $seqkit seq -r -p -w 0 | sed 's/>/>trans_/g'>> candidate_contig.fa &&
$seqkit seq -w0 candidate_contig.fa  |awk '{print $1}' | ${seqkit} seq > candidate_contig.fa.tmp &&
mv candidate_contig.fa.tmp candidate_contig.fa &&
$seqkit faidx candidate_contig.fa &&
echo end at `date` &&
touch 00.gapcloser_merge.sh.finish
"""
f.write(content)
f.close()
f=open('01.candidate_filter.sh','w')
if data_dict['category']!="f":
    content="""#!/bin/bash
echo begin  at `date`
seqkit={}
VGP_PIPELINE={}
sdust={}
Repetitive_monomer='{}'""".format(script_dict['seqkit'],script_dict['vgp'],script_dict['sdust'],data_dict['telomere_monomer'])+"""
file=candidate_contig.fa
prefix=telomere_candidate
$seqkit faidx candidate_contig.fa
$VGP_PIPELINE/telomere/find_telomere $file $Repetitive_monomer awk '{print $1"\\t"$(NF-4)"\\t"$(NF-3)"\\t"$(NF-2)"\\t"$(NF-1)"\\t"$NF}' - > ${prefix}.telomere
$sdust $file > ${prefix}.sdust
java -cp $VGP_PIPELINE/telomere/telomere.jar SizeFasta $file > $prefix.lens

java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereBreaks $prefix.lens $prefix.sdust $prefix.telomere > $prefix.breaks

awk '$4<1000{print $11}' $prefix.breaks | sort -u > tel_filter.id
$seqkit grep -f tel_filter.id candidate_contig.fa > candidate_contig.filter.fa
$seqkit faidx candidate_contig.filter.fa
touch 01.candidate_filter.sh.finish
echo end at `date`"""
    f.write(content)
else:
    content="""#!/bin/bash
echo begin  at `date`
seqkit={}
seqkt={}
Repetitive_monomer='{}'""".format(script_dict['seqkit'],script_dict['seqkt']data_dict['telomere_monomer'])+"""
file=candidate_contig.fa
$seqtk telo -s 20 -m ${Repetitive_monomer} candidate_contig.fa > candidate_contig.telo.bed 2> candidate_contig.telo.count
awk '$2<100{print $1}' candidate_contig.telo.bed | sort -u > tel_filter.id
$seqkit grep -f tel_filter.id candidate_contig.fa > candidate_contig.filter.fa
$seqkit faidx candidate_contig.filter.fa
touch 01.candidate_filter.sh.finish
echo end at `date`"""
	f.write(content)

f.close()

os.chdir('../01.check_without_telomere')
f=open('02.check_without_telomere.sh','w')
content="""#!/bin/bash
echo begin  at `date`
seqkit={}
VGP_PIPELINE={}
sdust={}
Repetitive_monomer='{}'
file=../00.gapcloser_merge/gap_free.fa
prefix={}
threshold=0.4
""".format(script_dict['seqkit'],script_dict['vgp'],script_dict['sdust'],data_dict['telomere_monomer'],data_dict['sp'])+"""
$VGP_PIPELINE/telomere/find_telomere $file $Repetitive_monomer awk '{print $1"\\t"$(NF-4)"\\t"$(NF-3)"\\t"$(NF-2)"\\t"$(NF-1)"\\t"$NF}' - > ${prefix}.telomere
$sdust $file > ${prefix}.sdust
java -cp $VGP_PIPELINE/telomere/telomere.jar SizeFasta $file > $prefix.lens

java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereBreaks $prefix.lens $prefix.sdust $prefix.telomere > $prefix.breaks

awk '$4<1000{print $11}' $prefix.breaks | sort -u > tel_head_chr.id
awk '($NF-$6)<1000{print $11}' $prefix.breaks | sort -u > tel_tail_chr.id
awk 'NR==FNR{a[$1]=$1}NR>FNR&&a[$1]!=$1{print $1":1-20000"}' tel_head_chr.id $prefix.lens | xargs $seqkit faidx ../00.gapcloser_merge/gap_free.fa | $seqkit seq -w 0 > telomere_without.fa
awk 'NR==FNR{a[$1]=$1}NR>FNR&&a[$1]!=$1{print $1":"$2-19999"-"$2}' tel_tail_chr.id $prefix.lens | xargs $seqkit faidx ../00.gapcloser_merge/gap_free.fa | $seqkit seq -w 0 -r -p >> telomere_without.fa

echo end at `date` &&
touch 02.check_without_telomere.sh.finish
"""
f.write(content)
f.close()

os.chdir('../02.minimap')
f=open('03.minimap.sh','w')
content="""#!/bin/bash
echo begin  at `date`
threshold=100
{} -t 20 -x {} ../01.check_without_telomere/telomere_without.fa ../00.gapcloser_merge/candidate_contig.filter.fa > candidate.paf""".format(script_dict['minimap'],data_dict['map_para'])+"""
awk -v t=$threshold '$3>t&&$5=="+"&&$8<t{print $10/$11"\\t"$0}' candidate.paf |sort -k 7,7 -k 11n -k 1n |awk '{a[$7]=$2"\\t"$11"\\t"$4}END{for(i in a)print i"\\t"a[i]}'|sort -k1,1 |awk '{print $2,$1,$3,$4}' > expand.info
echo end at `date` &&
touch 03.minimap.sh.finish
"""
f.write(content)
f.close()

os.chdir('../03.expand')
f=open('04.expand.sh','w')
content="""#!/bin/bash
echo begin  at `date`
python {2}/telomere_expand.py --info {0}/02.minimap/expand.info --genome {0}/00.gapcloser_merge/gap_free.fa --candidate {0}/00.gapcloser_merge/candidate_contig.filter.fa --map_para {1} --category {3}""".format(work_dir,data_dict['map_para'],os.path.dirname(os.path.realpath(__file__)),data_dict['category'])+"""
ls -d */ |awk '{print "cd "$1";sh work.sh;cd .."}' > run.sh
sh run.sh
echo end at `date` &&
touch 04.expand.sh.finish
"""
f.write(content)
f.close()

os.chdir('..')
f=open('shell.list','w')
f.write("""{0}/00.gapcloser_merge/00.gapcloser_merge.sh:4g:1
{0}/00.gapcloser_merge/00.gapcloser_merge.sh:4g:1\t{0}/00.gapcloser_merge/01.candidate_filter.sh:4g:1
{0}/00.gapcloser_merge/00.gapcloser_merge.sh:4g:1\t{0}/01.check_without_telomere/02.check_without_telomere.sh:4g:1
{0}/01.check_without_telomere/02.check_without_telomere.sh:4g:1\t{0}/02.minimap/03.minimap.sh:80g:20
{0}/02.minimap/03.minimap.sh:80g:20\t{0}/03.expand/04.expand.sh:4g:1
""".format(work_dir))
