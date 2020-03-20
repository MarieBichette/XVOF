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

        python -m unittest $module > "$(pwd)/result_$module.log"
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

        python -m unittest $module > "$(pwd)/result_$module.log" 
    fi
done

