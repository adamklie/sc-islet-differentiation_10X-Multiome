#!/bin/bash
#####
# Topic-level (DeepTopic / CREsted topic classification) motif marginalization.
# Uses the n100 51-motif catalog (combined_mapped.meme) and the best
# n100 checkpoint (epoch 60). Per-topic background sequences come from each
# topic's binarized BED file (binarized/beds/Topic{N}.bed).
#
# Submit as a 100-topic array (concurrency 4 to keep GPU busy without hogging):
#   /cellar/users/aklie/opt/SLURM/gpu_array.sh \
#     -s /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_topic_models/scripts/marginalization/sc-islet-differentiation_10X-Multiome_topic_marginalization.sh \
#     -j topic_marg -g a30:1 -c 4 -m 32G -t 14-00:00:00 -n 100 -x 4 \
#     -o /cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/slurm_logs/4_topic_models/marginalization
#####

date
echo -e "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID\n"

# CREsted Keras venv (same env used for contributions)
source /cellar/users/aklie/opt/deeptopic/.venv-keras/bin/activate
NVIDIA_LIBS=$VIRTUAL_ENV/lib/python3.11/site-packages/nvidia
export LD_LIBRARY_PATH=$NVIDIA_LIBS/cublas/lib:$NVIDIA_LIBS/cuda_runtime/lib:$NVIDIA_LIBS/cudnn/lib:$NVIDIA_LIBS/cufft/lib:$NVIDIA_LIBS/curand/lib:$NVIDIA_LIBS/cusolver/lib:$NVIDIA_LIBS/cusparse/lib:$NVIDIA_LIBS/nccl/lib:$NVIDIA_LIBS/nvjitlink/lib:$LD_LIBRARY_PATH
export KERAS_BACKEND=tensorflow

TOPIC="Topic${SLURM_ARRAY_TASK_ID}"
echo "Topic: $TOPIC"

SCRIPT=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/bin/4_topic_models/scripts/marginalization/marginalize_topic.py
BASE=/cellar/users/aklie/data/datasets/sc-islet-differentiation_10X-Multiome/results/4_topic_models/all/n100

n_seqs=100

python ${SCRIPT} \
    --model ${BASE}/crested_model/checkpoints/60.keras \
    --motif_file ${BASE}/crested_model/motifs/combined_mapped.meme \
    --genome_fasta /cellar/users/aklie/data/ref/genomes/hg38/hg38.fa \
    --peaks ${BASE}/binarized/beds/${TOPIC}.bed \
    --adata ${BASE}/crested_model/contributions_specific/adata_with_predictions.h5ad \
    --topic ${TOPIC} \
    --output_dir ${BASE}/crested_model/marginalization \
    --num_seqs ${n_seqs} \
    --batch_size 64

date
