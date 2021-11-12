#!/bin/bash
script_dir="$(dirname $(readlink -e $0))"
BATCH_ARGS=("--simulate" "--experimental-setup=BAMBU" "--expose-globals")
OUT_SUFFIX="bambu_specific_test6"

$script_dir/../../etc/scripts/test_panda.py --tool=bambu  \
   --args="--configuration-name=CLANG11_O0 -O0 --compiler=I386_CLANG11 ${BATCH_ARGS[*]}" \
   --args="--configuration-name=CLANG11_O1 -O1 --compiler=I386_CLANG11 ${BATCH_ARGS[*]}" \
   --args="--configuration-name=CLANG11_O2 -O2 --compiler=I386_CLANG11 ${BATCH_ARGS[*]}" \
   -lbambu_specific_test6_list \
   -o "output_${OUT_SUFFIX}" -b$script_dir \
   --name="${OUT_SUFFIX}" "$@"
exit $?

#  --args="--configuration-name=CLANG10_O0 -O0 --compiler=I386_CLANG10 ${BATCH_ARGS[*]}" \
#  --args="--configuration-name=CLANG10_O1 -O1 --compiler=I386_CLANG10 ${BATCH_ARGS[*]}" \
#  --args="--configuration-name=CLANG10_O2 -O2 --compiler=I386_CLANG10 ${BATCH_ARGS[*]}" \
