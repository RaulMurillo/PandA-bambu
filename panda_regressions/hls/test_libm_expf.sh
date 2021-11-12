#!/bin/bash
abs_script=$(readlink -e $0)
dir_script=$(dirname $abs_script)
if test -f output_test_libm_expf/finished; then
   exit 0
fi
if [[ -z "$J" ]]; then
   J="1"
fi
echo "Parallel jobs: $J"
rm -fr output_test_libm_expf
mkdir output_test_libm_expf
cd output_test_libm_expf
gcc -fopenmp -O3 -I$dir_script/../../etc/libbambu/ $dir_script/../../etc/libbambu/libm/hotbm_expf.c -DCHECK_EXP_FUNCTION -lm -lmpfr -lgmp
export OMP_NUM_THREADS=$J
./a.out
return_value=$?
cd ..
if test $return_value != 0; then
   echo "C based test of softfloat expf not passed."
   exit $return_value
fi
touch output_test_libm_expf/finished
exit 0
