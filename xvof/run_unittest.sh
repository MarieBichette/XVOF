#!/bin/bash

echo "----------------------------------------------------"
cd cell/test_cell

ListeModules="$(ls .)"
echo "Run UnitTest for cell modules"
for file in ${ListeModules}; do
    module=$(basename $file .py)
    module=$(basename $module .pyc)
    echo $module
    if [ "$module" != "test_variables" ] && [ "$module" != "__init__" ]; then
        echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo $module

        python2.7 -m unittest $module > "//home/marie/Documents/These/UnittestResults/result_$module.log" 2> "//home/marie/Documents/These/UnittestResults/error_$module.log"
        # un seul chevron pour écraser le fichier à chaque fois
        # 2> redirige le flux en erreur = sorties de unittest
    fi
done

cd ../..
cd node/test_node

ListeModules="$(ls .)"
echo "Run UnitTest for node modules"
for file in ${ListeModules}; do
    module=$(basename $file .py)
    module=$(basename $module .pyc)
    echo $module
    if [ "$module" != "test_variables" ] && [ "$module" != "__init__" ]; then
        echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo $module

        python2.7 -m unittest $module > "//home/marie/Documents/These/UnittestResults/result_$module.log" 2> "//home/marie/Documents/These/UnittestResults/error_$module.log"
        # un seul chevron pour écraser le fichier à chaque fois
        # 2> redirige le flux en erreur = sorties de unittest
    fi
done

