#!/usr/bin/env bash

testfunction() {
    local filename=${1}

    if [ -f "${filename}" ]; then
        echo "INFO: ${filename} already exists! Moving on..."
        return
    fi

    echo "still here with ${1}"

}
export -f testfunction

parallel -j 1 testfunction ::: $@