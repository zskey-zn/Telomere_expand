# -*- coding: utf-8 -*-
"""
Fix on Mon Dec 30 11:49:35 CST 2024
telomere expand

@author: zhengshang@frasergen.com   2398848440@qq.com

"""
import os
import sys
import re
import argparse
parser = argparse.ArgumentParser(description="telomere expand")
parser.add_argument('--info',help='expand infomation path',required=True)
parser.add_argument('--genome',help='gap-free genome path',required=True)
parser.add_argument('--candidate',help='candidate contig path',required=True)
parser.add_argument('--map_para',help='minimap2 parameter',default="map-hifi")
argv=vars(parser.parse_args())
info_file=argv['info'].strip()
genome_file=argv['genome'].strip()
candidate_file=argv['candidate'].strip()
map_para=argv['map_para'].strip()
script_path=os.path.dirname(os.path.realpath(__file__)) + "/software.txt"
f=open(script_path,'r')
data=f.readlines()
f.close()
data=[i[:-1] for i in data if i[0]!="#" and i!="\n"]
s_a=[i.split('=')[0].strip() for i in data]
s_b=[i.split('=')[1].strip() for i in data]
script_dict=dict(zip(s_a,s_b))

info=[]
with open(info_file) as obj:
    for line in obj:
        tmp=line.split()
        contig_id=tmp[0]
        if "cis" in contig_id:
            phase="+"
        else:
            phase="-"
        chr_id,chr_1,_=re.split("[:-]",tmp[1])
        if int(chr_1)>10000:
            pos_type='tail'
        else:
            pos_type='head'
        match_pos=int(tmp[3])
        info.append([chr_id,pos_type,contig_id,phase,match_pos])
        
expand_chr_list=list(set(i[0] for i in info))
for cell in expand_chr_list:
    cell_content=[i for i in info if i[0]==cell]
    os.system("mkdir {}".format(cell))
    os.chdir("./{}".format(cell))
    head_content=[i for i in cell_content if i[1]=="head"]
    tail_content=[i for i in cell_content if i[1]=="tail"]
    f=open("work.sh",'w')
    content="""#!/bin/bash
echo begin at `date`
seqkit={}
echo ">{}" > telomere_without.fa
touch candicate.fa;rm candicate.fa;
""".format(script_dict['seqkit'],cell)
    #head
    if len(head_content)!=0:
        content+="""$seqkit faidx {0} {1}:1-{2} | $seqkit seq -w {3} | head -2 | tail -1 >>  telomere_without.fa
echo {4} >> telomere_without.fa
$seqkit faidx {0} {1}:1-{5} >> candicate.fa
""".format(candidate_file,head_content[0][2],head_content[0][4],int(head_content[0][4]/100-1)*100,"N"*100,int(head_content[0][4])+20000)
    #gap-free chr
    content+="""echo {} > chr.id
$seqkit grep -f chr.id {} | $seqkit seq -w 0 | tail -1 >> telomere_without.fa
""".format(cell,genome_file)
    #tail
    if len(tail_content)!=0:
        content+="""echo {0} >> telomere_without.fa
$seqkit faidx {1} {2}:1-{3} | $seqkit seq -w {4} | head -2 |$seqkit seq -r -p -w 0 | tail -1 >>  telomere_without.fa
$seqkit faidx {1} {2}:1-{5} >> candicate.fa
""".format("N"*100,candidate_file,tail_content[0][2],tail_content[0][4],int(tail_content[0][4]/100-1)*100,int(tail_content[0][4])+20000)
    
    #gapcloser
    content+="""{} --scaff telomere_without.fa --reads candicate.fa --minmap_arg '-x {}' --min_idy 0.2 --min_match 200 --thread 20  --output gapclose --ne >pipe.log 2>pipe.err
echo end at `date`
touch work.sh.finish
""".format(script_path['GapCloser'],map_para)
    f.write(content)
    f.close()
    os.chdir("..")

