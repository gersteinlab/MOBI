######## 
# pwm to meme
# input: 
#     $1: input dir of pwm
#     $2: output dir of meme

module load MEME

mkdirhier $2

for file in $1/*.pwm; do
    sed "1d" $file | cut -f2- | matrix2meme -dna > $2/$(basename $file .pwm).meme
done