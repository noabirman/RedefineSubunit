#!/bin/bash
#SBATCH --job-name=ipsae_score
#SBATCH --exclude=sm-01,sm-16,sm-02,sm-03,sm-04,sm-08

#SBATCH --output=/cs/labs/dina/noabirman/slurms_outs/ipsae/ipsae_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=noa.birman@mail.huji.ac.il
#SBATCH --time=01:00:00        # 1 hour is likely plenty for hundreds of pairs
#SBATCH --mem=4G               # 4GB RAM is safe (script uses very little)
#SBATCH --cpus-per-task=1      # It is single-threaded

# Run the automation wrapper
python3 run_ipsae_on_af_pairs.py