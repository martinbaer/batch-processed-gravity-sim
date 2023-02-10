#!/bin/bash -l
#
#SBATCH --job-name=proj-cuda-m
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#The two line below are need to run code on the GPU 
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=00:00:01

#You could add these to your bashrc if you wanted
module load cuda
module load gnu

./bh_cuda input/constants.txt &> stdouterr_cuda.txt

#This what you could run if you wanted and interactive gpu session
#srun --cpus-per-task=1 --mem=2GB  --partition=gpu --gres=gpu:1 --pty bash 