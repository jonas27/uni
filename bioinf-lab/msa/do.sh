#!/bin/bash
set -e

# Start all tests
function tests {
    printf "_______Start: Run all tests._______\n"
    python test_msa.py 
    printf "_______End: RUn all tests._______\n\n\n"
}

function lint {
    printf "_______Start: lint all files._______\n"
    pylint msa.py | grep -v "doesn't conform to snake_case\|Too many"
    pylint test_msa.py | grep -v "doesn't conform to \|method docstring"
    printf "_______End: linting done._______\n"
}

function help {
    printf "Use for testing and linting.\n"
    printf 'Should be used like this: "./do.sh -t -l"\n'
}

while test $# != 0
do
    case "$1" in
    -t) tests ;;
    --test) tests ;;
    -l) lint ;;
    --lint) lint ;;
    -h) help ;;
    --help) help ;;
    esac
    shift
done
