#!/bin/bash
abs_script=$(readlink -e $0)
root_dir=$(dirname $abs_script)/../..
$root_dir/etc/scripts/characterize.py --technology-files=,,,,,,$root_dir/etc/lib/technology/C_FP_IPs_posit.xml,,, $@


