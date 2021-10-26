#!/bin/csh
#SBATCH -A thous_opt
#SBATCH -N 1
#SBATCH -t "30:00"
#SBATCH -p short
#SBATCH -J test_spds
#SBATCH -e test_spds.stderr.%j
#SBATCH -o test_spds.stdout.%j
module purge
source env.csh
unlimit
unlimit coredumpsize
echo
echo "module list"
echo
module list >& tmp_ml.lis
cat tmp_ml.lis
/bin/rm -rf tmp_ml.lis
echo
echo "printenv output:"
echo
printenv
echo
echo "limit output:"
echo
limit
echo
echo "ldd -r output:"
echo
ldd -r ./era_opt
echo
echo era_opt output
echo
mpirun -n 24 ./era_opt f100.in

