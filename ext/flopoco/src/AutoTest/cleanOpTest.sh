rm -R AutotestResults/tmp
nbTests=$(grep -c './flopoco' AutotestResults/$1/report)
nbErrors=$(grep -c ERROR AutotestResults/$1/report)
nbVHDL=$(grep -c "VHDL generated"  AutotestResults/$1/report)
nbSuccess=$((nbVHDL-nbErrors))
rateError=$(((nbErrors*100)/nbTests))
rateSuccess=$(((nbSuccess*100)/nbTests))
rateVHDL=$(((nbVHDL*100)/nbTests))

echo "Operator $1:\n   $nbTests tests; \t VHDL generation OK for $nbVHDL ($rateVHDL%); \t test OK for $nbSuccess ($rateSuccess%) " >> AutotestResults/report
echo "$1: \t $nbTests tests; \t VHDL generation OK for $nbVHDL ($rateVHDL%); \t test OK for $nbSuccess ($rateSuccess%) "

