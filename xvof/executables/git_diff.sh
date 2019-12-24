#!/bin/bash

ListeRep="$(ls -d */)"
for my_dir in ${ListeRep}; do 
    cd $my_dir
    echo "----------------------------------------------------" >> '//home/marie/Documents/result_git_diff.txt'
    echo "Analyse du repertoire : $my_dir" >> '//home/marie/Documents/result_git_diff.txt'
    echo "----------------------------------------------------" >> '//home/marie/Documents/result_git_diff.txt'
    git diff --src-prefix="echec/" --dst-prefix="reussite/"  aa2a189 815dfb3 -- "./" >> '//home/marie/Documents/result_git_diff.txt'
    cd ..
done
