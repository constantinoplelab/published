#! /bin/bash
#SBATCH --job-name=SHUFF_allofc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%j.out
#SBATCH --account=torch_pr_467_general

# load matlab
module load matlab/2025b 

#where output will be pla ed
savedir=/scratch/dh148/dynamics/results/maggie/
datadir=/scratch/dh148/dynamics/data/maggie/
rundir=/home/dh148/projects/constantinoplelab/Analysis/david/BiasModel/opto_hetero/
shuff2use=SHUFF
cd $rundir


for epoch in reward coff son
#for epoch in son
do
    #for decodertype in psth2
    for decodertype in psth2 none
    do
	usemlb=true
	dataname=\'$datadir/preprocess_mlb_ofc_${epoch}_${decodertype}.mat\'
	savename=\'$savedir/mlb_allofc_${shuff2use}_${epoch}_${decodertype}.mat\'
	cmd="run_decoding_mlb_allofc($dataname,$savename,$usemlb,$shuff2use)"
	echo $cmd
	matlab -batch "$cmd"

	usemlb=false
	dataname=\'$datadir/preprocess_mlb_ofc_${epoch}_${decodertype}.mat\'
        savename=\'$savedir/mlb2_allofc_${shuff2use}_${epoch}_${decodertype}.mat\'
        cmd="run_decoding_mlb_allofc($dataname,$savename,$usemlb,$shuff2use)"
        echo $cmd
        matlab -batch "$cmd" 
    done
done

	



