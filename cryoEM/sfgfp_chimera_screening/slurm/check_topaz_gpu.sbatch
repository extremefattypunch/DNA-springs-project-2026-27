#!/bin/bash -l
#SBATCH -J sfgfp-topaz-gpu
#SBATCH -p gpu_test
#SBATCH --gres=gpu:1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0-00:10
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
nvidia-smi -L
env -u PYTHONPATH -u LD_LIBRARY_PATH /n/holylabs/hekstra_lab/Everyone/ianpoon/software/topaz-0.3/bin/python -c "
import torch; print('cuda available', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else '')
x=torch.randn(2048,2048,device='cuda'); print('matmul ok', float((x@x).sum())!=0)"
