#!/usr/bin/env bash

#set -x # show each running command
# @File    :   DNA-Seq.sh
# @Time    :   2021/05/13 16:02:09
# @Author  :   Wang Xiangeng 
# @Version :   1.0
# @Contact :   meinetbosekampf@gmail.com
# @License :   (C)Copyright 2017-2018, Liugroup-NLPR-CASIA
# should be run inside /local_data/WGS
# 一些心得体会：
#   删除一定谨慎明确

########################################
#  
# GLOBAL CONTROL
# 
########################################

# shopt -s extglob
# touch /data/Manager/DNA-Seq_Manage/Finished_WGS_PREP
# touch /data/Manager/DNA-Seq_Manage/Finished_GERMLINE
touch /data/Manager/DNA-Seq_Manage/CHECK
#touch /data/Manager/DNA-Seq_Manage/FQ
#touch /data/Manager/DNA-Seq_Manage/BQSR
#ALLOC=/data/Manager/DNA-Seq_Manage/WGS$(printf "%02d" $(hostname|cut -f2 -d"-")) # 自动获取对应的文件安排
ALLOC=/data/Manager/DNA-Seq_Manage/WGS
FINISHED=/data/Manager/DNA-Seq_Manage/CHECK # 记录分析完成的样本名
#GERMLINE=/data/Manager/DNA-Seq_Manage/Finished_GERMLINE
RESQ=/data/Manager/DNA-Seq_Manage/RESQ
RESULT=/data/Results/DNA-Seq
#FQ_UPLOADED=/data/Manager/DNA-Seq_Manage/FQ
#BQSRED=/data/Manager/DNA-Seq_Manage/BQSR


