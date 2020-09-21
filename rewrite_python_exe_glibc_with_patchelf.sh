#!/usr/bin/env bash

# TODO edit this line to specify location of new glibc 
#export GLIBC_PATH=/cluster/tufts/hugheslab/miniconda2/envs/ape/custom_libs/
#export GLIBC_LD_PATH=$GLIBC_PATH/lib/x86_64-linux-gnu/ld-2.23.so
export GLIBC_PATH=/home/gomeza/install/glibcpath/glibc-2.27
export GLIBC_LD_PATH=/home/gomeza/install/glibcpath/glibc-2.27/lib/ld-2.27.so

if [[ ! -f $GLIBC_LD_PATH ]]; then
    echo "ERROR: Provided GLIBC_LD_PATH not valid"
    exit
fi

echo "OVERWRITING PYTHON EXECUTABLE:"
python_exe=`which python`
echo $python_exe

IS_CONDA_ENV=`python -c "print('$python_exe'.count('/envs/') > 0)"`
echo "IS_CONDA_ENV: $IS_CONDA_ENV"

if [[ $IS_CONDA_ENV -ne 'True' ]]; then
    echo "ERROR: Current python executable not in conda env. Will not alter to avoid problems."
    exit
fi

CONDA_ENV_LIB=`python -c "print('$python_exe'.replace('/bin/python', '/lib'))"`

echo "CREATING BACKUP PYTHON"
python_tmp_exe=`python -c "print('$python_exe'.replace('python', 'python_backup'))"`
cp $python_exe $python_tmp_exe
echo "$python_tmp_exe"

#rpath=$GLIBC_PATH/lib/x86_64-linux-gnu:$CONDA_ENV_LIB:/usr/lib64:/lib64:/lib
rpath=$GLIBC_PATH/lib:$CONDA_ENV_LIB:/usr/lib64:/lib64:/lib

echo "CALLING PATCHELF on 'python' binary"
patchelf --set-interpreter $GLIBC_LD_PATH --set-rpath $rpath $python_exe
echo "DONE! patchelf complete"

