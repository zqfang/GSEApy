#!/bin/bash

# Installs gseapy and requirements into a fresh Python 2 or 3 environment
# and runs tests.
#
set -e

PY_VERSION=$1

usage="Usage: $0 py_version[2|3]"
: ${PY_VERSION:?$usage}

log () {
    echo
    echo "[`date`] TEST HARNESS: $1"
    echo
}

log "removing existing env pbtpy${PY_VERSION}"
name=pbtpy${PY_VERSION}
conda env list | grep -q $name && conda env remove -y -n $name

log "starting with basic environment"
conda create -y -n $name  python=${PY_VERSION}
source activate $name

log "temporarily install cython"
conda install cython

log "force re-cythonizing"
rm -rf dist build
python setup.py clean
python setup.py build
python setup.py sdist

log "uninstall cython"
conda remove cython

log "test installation of sdist"
set -x
(cd dist && pip install gseapy-*.tar.gz && python -c 'import gseapy')
set +x

python setup.py clean

log "install test requirements"
source deactivate
conda env list | grep -q $name && conda env remove -y -n $name
conda create -y -n $name  python=${PY_VERSION} \
    --file "requirements.txt" \
    --file "test-requirements.txt" 
	
source activate $name

log "run tests"
nosetests
# #(cd docs && make clean && make doctest)

