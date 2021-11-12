#!/bin/bash
script=$(readlink -e $0)
root_dir=$(dirname $script)

mkdir -p constrained_synth_altera
cd constrained_synth_altera
echo "# Quartus II synthesis and ICARUS simulation"
timeout 2h bambu $root_dir/module.c --generate-tb=$root_dir/test.xml --simulator=ICARUS --device-name=EP2C70F896C6 --evaluation --experimental-setup=BAMBU --generate-interface=WB4 $root_dir/constraints_STD.xml  --cprf=0.9 --skip-pipe-parameter=1 "$@"
return_value=$?
if test $return_value != 0; then
   exit $return_value
fi
cd ..
