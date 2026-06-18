#!/bin/bash
# Train DeepTopic CNN using CREsted with Keras/TF backend
date
echo -e "Job ID: $SLURM_JOB_ID\n"

conda deactivate 2>/dev/null
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate

NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH
export KERAS_BACKEND=tensorflow

cmd="python /cellar/users/aklie/opt/deeptopic/scripts/crested/train_topic_model.py \
    --beds_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine/binarized/beds \
    --regions_file /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine/binarized/consensus_regions.bed \
    --genome /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --output_dir /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/endocrine/crested_model \
    --project_name sc-islet-differentiation_10X-Multiome_endocrine \
    --backend tensorflow"
echo -e "Running:\n$cmd\n"
eval $cmd

date
