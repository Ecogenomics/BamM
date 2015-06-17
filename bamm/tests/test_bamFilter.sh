#
###############################################################################
#                                                                             #
#    test_bamFilter.sh                                                        #
#                                                                             #
#    Test output of BamM 'filter' command                                     #
#                                                                             #
#    Copyright (C) Tim Lamberton                                              #
#                                                                             #
#    tim.lamberton@gmail.com                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This library is free software; you can redistribute it and/or            #
#    modify it under the terms of the GNU Lesser General Public               #
#    License as published by the Free Software Foundation; either             #
#    version 3.0 of the License, or (at your option) any later version.       #
#                                                                             #
#    This library is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        #
#    Lesser General Public License for more details.                          #
#                                                                             #
#    You should have received a copy of the GNU Lesser General Public         #
#    License along with this library.                                         #
#                                                                             #
###############################################################################

data_dir=~/git/BamM/bamm/tests/filter_test_data
bamm_exe=~/git/BamM/bin/bamm

### arguments
make=
f=
if [ "$1" == "make" ]
then
  bam_id=$2
  make=1
elif [ "$1" == "fail" ]
then
  bam_id=$2
  f='f'
else
  bam_id=$1
fi
if [ ! "$bam_id" ]
then echo "ERROR: Missing test .bam index" >&2; exit 1;
fi

### Helpers ###
output_dir=$data_dir/$bam_id/output
OUTPUT_FILE=$output_dir/${bam_id}_filtered.bam
FAILED=
v=

run_test() {
  CURRENT_TEST="$bamm_exe filter -b $data_dir/$bam_id$f.bam -o $output_dir $@"
  if [ $v ]
  then echo "RUNNING: $CURRENT_TEST" >&2
  fi
  $CURRENT_TEST
  if [ $make ]
  then
    mv $OUTPUT_FILE $TEST_FILE
  else
    assert_equal_length && pass || fail
    rm $OUTPUT_FILE
  fi
}

assert_equal_length() {
  size1=$(samtools view $OUTPUT_FILE |wc -l)
  size2=$(samtools view $TEST_FILE |wc -l)
  if [ $size1 != $size2 ]
  then return 1
  else return 0
  fi
}

pass() {
  if [ $v ]
  then echo "PASSED: $CURRENT_TEST" >&2
  fi
}

fail() {
  echo "FAILED: $CURRENT_TEST" >&2
  FAILED=1;
}

tests_done() {
  if [ $make ]
  then
    echo "Finished generating data for test set $bam_id." >&2
    exit 0
  elif [ $FAILED ]
  then
    echo "There were failing tests for set $bam_id.">&2
    exit 1
  else
    echo "All tests passed for set $bam_id." >&2
    exit 0
  fi
}

test_file() {
  TEST_FILE=$data_dir/$bam_id/${bam_id}_$1.bam
}

### Test 1 ###
test_file none
run_test --use_secondary --use_supplementary --percentage_aln 0 --percentage_id 0

### Test 2 ###
test_file aln_only_0.9
run_test --use_secondary --use_supplementary --percentage_aln 0.9

### Test 3 ###
test_file aln_only_1.01
run_test --use_secondary --use_supplementary --percentage_aln 1.01

### Test 4 ###
test_file id_only_0.9
run_test --use_secondary --use_supplementary --percentage_id 0.9

### Test 5 ###
test_file id_only_1.01
run_test --use_secondary --use_supplementary --percentage_id 1.01

### Test 6 ###
test_file no_secondary_only
run_test --use_supplementary

### Test 7 ###
test_file no_supp_only
run_test --use_secondary

### Test 8 ###
test_file all_conds
run_test --percentage_aln 0.9 --percentage_id 0.9

### Clean up
rmdir $output_dir
tests_done
