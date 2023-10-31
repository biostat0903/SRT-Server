#!/usr/bin/bash
#PBS -N ref_train
#PBS -l nodes=node3:ppn=1,mem=50g
#PBS -l walltime=360:00:00
#PBS -p +1023
#PBS -t 1
#PBS -j oe
#PBS -o /public/home/biostat03/project/stwebProject/00_temp/ref_train_test.out
#PBS -q batch

#
bash
let k=0

project_path=/public/home/biostat03/project/stwebProject/
ref_file=${project_path}02_data/reference_data/SRT-Server/scRNA-seq/ref_df.txt
covertR=${project_path}01_code/srt_server/covert_seurat.R
trainpy=${project_path}01_code/srt_server/ref_train.py

for nref in 26
do
let k=${k}+1
if [ ${k} -eq ${PBS_ARRAYID} ]
then

refx=`head -n ${nref} ${ref_file} | tail -n 1`
spec=`echo ${refx} | awk '{print $1}'`
ref=`echo ${refx} | awk '{print $2}'`
#
Rscript ${covertR} --spec ${spec} --ref ${ref}
#
source activate /public/home/biostat04/anaconda3/envs/cell2loc_env
path=/public/home/biostat03/project/stwebProject/01_code/case_study
ref_train_file=${path}/case_study1_ref_train.txt
/usr/bin/time -v -o ${ref_train_file} python ${trainpy} --spec ${spec} --ref ${ref}
conda deactivate

fi
done

