#!/bin/bash
set -e

# Start all tests
function tests {
    printf "_______Start: Run all tests._______\n"
    python test_gotoh_unit.py 
    python test_gotoh_assignment1_dna.py
    python test_gotoh_assignment2_proteins.py
    printf "_______End: RUn all tests._______\n\n\n"
}

function lint {
    printf "_______Start: lint all files._______\n"
    pylint gotoh.py | grep -v "doesn't conform to snake_case"
    pylint gotoh_helpers.py | grep -v "doesn't conform to snake_case"
    pylint test_gotoh_unit.py | grep -v "doesn't conform to snake_case"
    pylint test_gotoh_assignment1_dna.py | grep -v "doesn't conform to snake_case"
    pylint test_gotoh_assignment2_proteins.py | grep -v "doesn't conform to snake_case"
    printf "_______End: linting done._______\n"
}

while test $# != 0
do
    case "$1" in
    -t) tests ;;
    --test) tests ;;
    -l) lint ;;
    --lint) lint ;;
    # -A) all_into_one=t
    #     unpack_unreachable=--unpack-unreachable ;;
    # -d) remove_redundant=t ;;
    # -q) GIT_QUIET=t ;;
    # -f) no_reuse=--no-reuse-object ;;
    # -l) local=--local ;;
    # --max-pack-size|--window|--window-memory|--depth)
    #     extra="$extra $1=$2"; shift ;;
    # --) shift; break;;
    # *)  usage ;;
    esac
    shift
done
