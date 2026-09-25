#! /bin/bash
#SBATCH --job-name=epoch_decodetype
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
#SBATCH --account=torch_pr_467_general

# load matlab
module load matlab/2025b 

#where output will be pla ed
savedir=/scratch/dh148/dynamics/results/maggie/20260602_decode_projection/
rundir=/home/dh148/projects/constantinoplelab/Analysis/david/BiasModel/opto_hetero/

epoch=reward
decodertype=psth

echo "starting matlab job, MLB"
cd $rundir
usemlb=true
savename=\'$savedir/mlb_${epoch}_${decodertype}.mat\'
echo "run_decoding_mlb($savename, $usemlb, '$epoch', '$decodertype','true')"
matlab -batch "run_decoding_mlb($savename, $usemlb, '$epoch', '$decodertype','true')"

usemlb=false
savename=\'$savedir/mlb2_${epoch}_${decodertype}.mat\'
echo "run_decoding_mlb($savename, $usemlb, '$epoch', '$decodertype','true')"
matlab -batch "run_decoding_mlb($savename, $usemlb, '$epoch', '$decodertype','true')"



echo "done"

