#!/bin/bash

# remove train and test if they exist
rm -rf ./train ./test
# create directories for training and testing
mkdir ./test
# create copy to remove later of the dataset
cp -r ./recordings ./train
# loop and move 10 files per each speacker in test
digits=$(seq 0 1 9)
speakers=("george" "jackson" "lucas" "nicolas" "theo" "yweweler")
for digit in ${digits[@]}; do
    for speaker in ${speakers[@]}; do
        find ./train/ -iname "*$digit\_$speaker\_*" | sort -R | tail -10 | xargs -Ixx mv xx ./test
    done
done

# reorganize train and test into subfolders per digit spoken
for digit in ${digits[@]}; do
    mkdir -p ./train/$digit ./test/$digit
    find ./train -iname "*$digit\_*" | xargs -Ixx mv xx ./train/$digit
    find ./test -iname "*$digit\_*" | xargs -Ixx mv xx ./test/$digit
done





