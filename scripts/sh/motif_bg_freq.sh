# $1: motif dir
# $2: species, current support "human", "worm", "fly"

if [ $2 == "fly" ]; then 
    for file in $1/*.meme ; do
        sed -i 's/A 0.25 C 0.25 G 0.25 T 0.25/A 0.29131 C 0.20869 G 0.20869 T 0.29131/g' $file
    done
fi

if [ $2 == "worm" ]; then 
    for file in $1/*.meme ; do
        sed -i 's/A 0.25 C 0.25 G 0.25 T 0.25/A 0.32280 C 0.17720 G 0.17720 T 0.32280/g' $file
    done
fi

if [ $2 == "human" ]; then 
    for file in $1/*.meme ; do
        sed -i 's/A 0.25 C 0.25 G 0.25 T 0.25/A 0.29528 C 0.20472 G 0.20472 T 0.29528/g' $file
    done
fi