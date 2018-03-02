for l in `seq 10 24`;
do
    for run in `seq 1 10000`;
    do 
        for h in `seq 0 11`;
        do
            ../build/src/hashingtest $(od -A n -t u -N 4 /dev/urandom) $h $((2**$l)) | tee -a $HOSTNAME-keys.txt
        done;
    done;
done;
