#!/bin/bash
script=$(readlink -e $0)
root_dir=$(dirname $script)

mkdir -p icrc
cd icrc
echo "#synthesis of icrc"
timeout 2h bambu $root_dir/spec.c --top-fname=icrc --simulator=VERILATOR --device-name=LFE335EA8FN484C --simulate --generate-tb=$root_dir/test_icrc.xml --channels-type=MEM_ACC_NN --experimental-setup=BAMBU return_value=$?
if test $return_value != 0; then
   exit $return_value
fi
cd ..

mkdir -p main
cd main
echo "#synthesis of main"
timeout 2h bambu $root_dir/spec.c  --simulator=VERILATOR  --device-name=LFE335EA8FN484C --simulate --generate-tb=$root_dir/test.xml --channels-type=MEM_ACC_NN --experimental-setup=BAMBU return_value=$?
if test $return_value != 0; then
   exit $return_value
fi


