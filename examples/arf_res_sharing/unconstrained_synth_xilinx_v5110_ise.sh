#!/bin/bash
script=$(readlink -e $0)
root_dir=$(dirname $script)

rm -rf unconstrained_synth_xilinx_v5110_ise
mkdir -p unconstrained_synth_xilinx_v5110_ise
cd unconstrained_synth_xilinx_v5110_ise
echo "# ISE synthesis and ICARUS simulation"
timeout 2h bambu $root_dir/module.c --generate-tb=$root_dir/test.xml --simulator=ICARUS --device-name=xc5vlx110t,-1,ff1136 --evaluation --experimental-setup=BAMBU --generate-interface=WB4 --clock-period=5 --cprf=0.9 --skip-pipe-parameter=1 "$@"
return_value=$?
if test $return_value != 0; then
   exit $return_value
fi
cd ..


