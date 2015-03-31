while [ 1 ]
do
export CILK_NWORKERS=4; (time -p ./main 1 1000000000 2>&1 | grep "real|user|system") 2>> test_output_$(hostname).out
done
