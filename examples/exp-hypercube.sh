for run in `seq 1 10000`;
do 
    for h in `seq 0 11`;
    do
        ../build/src/hashingtest $(od -A n -t u -N 4 /dev/urandom) $h | tee -a $HOSTNAME-hypercube.txt
    done;
done;
