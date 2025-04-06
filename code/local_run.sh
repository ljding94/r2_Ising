L=100
T=1.0
sigma=0.1
init="random"
run=0
folder="/Users/ldq/Work/r2_Ising_RandH/data/data_local/data_pool"

#for beta in 0.5 1.0 2.0 4.0 6.0 8.0 10.0 20.0; do
for T in 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7; do
    for sigma in 0.0; do
        nohup ./r2_Ising $L $T $sigma $init $run  $folder&
    done
done