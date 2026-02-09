#!/bin/bash

# Run 2D T Cell experiments for varying cell deformation frequency, omega

DIR="./experiments"

NP="300" # num particles
DT="0.005" # time step
T_FINAL="1500" # run time
WRITEEVERY="10000" # Write system to file every x timesteps
F_REP="20"
F_ADH="0.75"
AMP_RAD="0.15"
L="43"

REALIZATIONS="10"

frequencies="30"

mkdir $DIR

g++ main.cpp # Compile - can remove if you dont wanna compile every time

# Run experiments while varying omega
for freq in $frequencies; do
    for ((seed = 0; seed < $REALIZATIONS; seed++)); do
	echo "Running for omega = $freq"
	
	./a.out --np $NP --L $L --seed $seed --dt $DT --tfinal $T_FINAL --writeEvery $WRITEEVERY --w $freq --f_rep $F_REP --f_adh $F_ADH --amp_rad $AMP_RAD
	
	touch config.yaml
	echo "np: $NP" >> config.yaml
	echo "seed: $SEED" >> config.yaml
	echo "dt: $DT" >> config.yaml
	echo "t_final: $T_FINAL" >> config.yaml
	echo "write_every: $WRITEEVERY" >> config.yaml
	echo "frequency: $freq" >> config.yaml
	echo "f_rep: $F_REP" >> config.yaml
	echo "f_adh: $F_ADH" >> config.yaml
	echo "amp_rad: $AMP_RAD" >> config.yaml
	

	LOCAL_PATH=$DIR/w=$freq/$seed
	
	mkdir -p $LOCAL_PATH
	mv Res* ./$LOCAL_PATH
	mv First* ./$LOCAL_PATH	
	mv config.yaml ./$LOCAL_PATH
    done
done




