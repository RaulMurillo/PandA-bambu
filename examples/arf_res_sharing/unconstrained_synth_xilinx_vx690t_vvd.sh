#!/bin/bash
script=$(readlink -e $0)
root_dir=$(dirname $script)

rm -rf unconstrained_synth_xilinx_vx690t_vvd
mkdir -p unconstrained_synth_xilinx_vx690t_vvd
cd unconstrained_synth_xilinx_vx690t_vvd
echo "# Vivado synthesis and ICARUS simulation"
timeout 2h bambu $root_dir/module.c --generate-tb=$root_dir/test.xml --simulator=ICARUS --evaluation --experimental-setup=BAMBU --generate-interface=WB4 --clock-period=5 --cprf=0.9 --skip-pipe-parameter=1 --device-name=xc7vx690t,-3,ffg1930,VVD "$@"
return_value=$?
if test $return_value != 0; then
   exit $return_value
fi
cd ..


