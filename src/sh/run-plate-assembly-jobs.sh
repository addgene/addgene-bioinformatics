#!/bin/bash

# Print usage
usage() {
    cat << EOF

NAME
     run-plate-assembly-job -- Assembles each well of selected plates

SYNOPSIS
     run-plate-assembly-job [-d output-directory] [-st] assembler

DESCRIPTION
     Assembles each well of selected plates using the specified
     assembler, one of 'masurca', 'novoplasty', 'shovill', 'skesa',
     'spades', or 'unicycler'. Output files are placed in the current
     working directory, or in the directory specified as an option.

OPTIONS
-d   Directory in which to place assembler output files. Default is
     selected assembler.
-s   Date and time stamp the ouput directory
-t   Use a test plate containing two wells

EOF
}

# Parse command line options
OUTPUT_DIRECTORY=""
USE_DATE_STAMP=0
USE_TEST_PLATE=0
while getopts ":d:sth" opt; do
    case $opt in
	d)
	    OUTPUT_DIRECTORY=$OPTARG
	    ;;
	s)
	    USE_DATE_STAMP=1
	    ;;
	t)
	    USE_TEST_PLATE=1
	    ;;
	h)
	    usage
	    exit 0
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG" >&2
	    usage
	    exit 1
	    ;;
	\:)
	    echo "Option -$OPTARG requires an argument." >&2
	    usage
	    exit 1
	    ;;
    esac
done

# Parse command line arguments
shift `expr $OPTIND - 1`
if [ "$#" -ne 1 ]; then
    echo "One argument is required: exiting"
    exit 1
fi
ASSEMBLER=$1

# Exit immediately if a simple command exits with a non-zero status
set -e

# Determine the repository base directory
BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed s%/src/sh%% )"

# Activate the virtual environment
source $HOME/.virtualenvs/addgene-bioinformatics/bin/activate

# Work in the specified directory, or if a directory is not specified,
# in a directory with the same name as the assembler
if [ -z "$OUTPUT_DIRECTORY" ]; then
    OUTPUT_DIRECTORY=$ASSEMBLER
fi
if [ $USE_DATE_STAMP == 1 ]; then
    DATESTAMP=`date "+%Y-%m-%d-%H%M%S"`
    OUTPUT_DIRECTORY=$OUTPUT_DIRECTORY-$DATESTAMP
fi
mkdir -p $OUTPUT_DIRECTORY
pushd $OUTPUT_DIRECTORY

# Specify plates ...
PLATES=""
if [ $USE_TEST_PLATE == 1 ]; then

    # ... for testing
    PLATES="$PLATES A11935X_sW0148"

else
    # ... for processing

    PLATES="$PLATES A11935_sW0148"
    PLATES="$PLATES A11938A_sW0144"
    PLATES="$PLATES A11942A_sW0145"
    PLATES="$PLATES A11948_sW0148"
    PLATES="$PLATES A11953A_sW0150"

    PLATES="$PLATES A11956A_sW0152"
    PLATES="$PLATES A11957_sW0148"
    PLATES="$PLATES A11959A_sW0150"
    PLATES="$PLATES A11960A_sW0154"
    PLATES="$PLATES A11966A_sW0152"

    # PLATES="$PLATES A11967A_sW0154"
    # PLATES="$PLATES A11967B_sW0154"
    # PLATES="$PLATES A11970A_sW0154"
    # PLATES="$PLATES A11978B_sW0154"
    # PLATES="$PLATES A11979_sW0154"

    # PLATES="$PLATES A11980A_sW0155"
    # PLATES="$PLATES A11983_sW0157"
    # PLATES="$PLATES A11984_sW0158"
    # PLATES="$PLATES A11985A_sW0155"
    # PLATES="$PLATES A11988A_sW0155"

fi
    
# Process each plate
for PLATE in $PLATES; do

    # Remove the plate assembly job file store, since it exists if an
    # earlier plate assembly job failed
    rm -rf pajfs

    # Run the plate assembly job
    python ${BASE}/src/python/PlateAssemblyJob.py \
	   -s s3 -d addgene-sequencing-data/2018/FASTQ \
	   -p $PLATE \
	   -a $ASSEMBLER \
	   pajfs

    # Archive the plate assembly job output files
    tar -czvf $PLATE.tar.gz $PLATE

done

# Deactivate the virtual environment
deactivate

# Return to the original working directory
if [ -n "$OUTPUT_DIRECTORY" ]; then
    popd
fi

# Archive the plate assembly jobs output directory
tar -czvf $OUTPUT_DIRECTORY.tar.gz \
    $OUTPUT_DIRECTORY --exclude=*.tar.gz
