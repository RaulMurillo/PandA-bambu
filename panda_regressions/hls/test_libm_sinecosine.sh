#!/bin/bash
abs_script=$(readlink -e $0)
dir_script=$(dirname $abs_script)
if test -f output_test_libm_sinecosine/finished; then
   exit 0
fi
if [[ -z "$J" ]]; then
   J="1"
fi
echo "Parallel jobs: $J"
rm -fr output_test_libm_sinecosine
mkdir output_test_libm_sinecosine
cd output_test_libm_sinecosine
gcc -fopenmp -O3 -I$dir_script/../../etc/libbambu/ $dir_script/../../etc/libbambu/libm/hotbm_sine_cosine.c -DCHECK_TRIG_FUNCTIONS -lm -lmpfr -lgmp
export OMP_NUM_THREADS=$J
./a.out
return_value=$?
cd ..
if test $return_value != 0; then
   echo "C based test of softfloat sine/cosine not passed."
   exit $return_value
fi
touch output_test_libm_sinecosine/finished
exit 0
