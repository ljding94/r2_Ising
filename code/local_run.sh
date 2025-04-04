L=100
beta=5.0
sigma=0.1
init="random"
run=0
folder="/Users/ldq/Work/r2_Ising_RandH/data/data_local"

for beta in 0.2 0.5 1.0 5.0 10.0 20.0 50.0; do
    for sigma in 0.0 0.3; do
        nohup ./r2_Ising $L $beta $sigma $init $run  $folder&
    done
done