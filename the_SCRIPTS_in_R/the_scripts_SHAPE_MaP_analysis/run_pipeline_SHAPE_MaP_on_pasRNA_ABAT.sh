#!/bin/bash
#----------------------------------------------------------------------#
# This file is a part of ShapeMapper, and is licensed under the terms  #
# of the MIT license. Copyright 2018 Steven Busan.                     #
#----------------------------------------------------------------------#

# Run simple ShapeMapper pipeline on a small subset of example data,
# outputting additional information for power users interested in
# running individual shapemapper modules "manually".

set -e # exit on first error (if any)

# Find the parent folder of this script,
# resolving (possibly nested) symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$BASE_DIR/$SOURCE"
done
BASE_DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

export PATH=${BASE_DIR}:${PATH}

cd ${BASE_DIR}

shapemapper \
--name "the_results_pasRNA_ABAT" \
--target example_data_ABAT/ABAT.fa \
--amplicon \
--primers the_primers_SHAPE_MaP.txt \
--overwrite \
--min-depth 1000 \
--modified --folder example_data_ABAT/ABATplus \
--untreated --folder example_data_ABAT/ABATminus \
--verbose \
--output-parsed-mutations \
--output-counted-mutations
