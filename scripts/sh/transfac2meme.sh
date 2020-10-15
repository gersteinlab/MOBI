######## 
# transfac to meme
# input: 
#     $1: input dir of transfac
#     $2: output dir of meme

module load MEME

mkdirhier $2

for file in $1/*.PWM.txt; do
    transfac2meme -number $file > $2/$(basename $file .PWM.txt).meme
done