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
pre_bam=$data_dir/pre.bam
output_dir=$data_dir/output
bamm_exe=~/git/BamM/bin/bamm

### Helpers ###
OUTPUT_FILE=$output_dir/$(basename $pre_bam .bam)_filtered.bam
FAILED=

run_test() {
  CURRENT_TEST="$bamm_exe filter -b $pre_bam -o $output_dir $@"
  $CURRENT_TEST
  assert_equal_length
  rm $OUTPUT_FILE
}

assert_equal_length() {
  size1=samtools view $OUTPUT_FILE |wc -l
  size2=samtools view $TEST_FILE |wc -l
  if [ $size1 != $size2 ]; then fail_test; fi
}


fail_test() {
  echo "FAILED: $CURRENT_TEST" >&2
  FAILED=1;
}

tests_done() {
  if [ $FAILED ]
  then
    echo "There were failing test.">&2
    exit 1
  else
    echo "All tests passed." >&2
    exit 0
  fi
}

### Test 1 ###
TEST_FILE=$data_dir/none.bam
run_test --use_secondary --use_supplementary --percentage_aln 0 --percentage_id 0

### Test 2 ###
TEST_FILE=$data_dir/aln_only_0.9.bam
run_test --use_secondary --use_supplementary --percentage_aln 0.9

### Test 3 ###
TEST_FILE=$data_dir/aln_only_1.01.bam
run_test --use_secondary --use_supplementary --percentage_aln 1.01

### Test 4 ###
TEST_FILE=$data_dir/id_only_0.9.bam
run_test --use_secondary --use_supplementary --percentage_id 0.9

### Test 5 ###
TEST_FILE=$data_dir/id_only_1.01.bam
run_test --use_secondary --use_supplementary --percentage_id 1.01

### Test 6 ###
TEST_FILE=$data_dir/no_second_only.bam
run_test --use_supplementary

### Test 7 ###
TEST_FILE=$data_dir/no_supp_only.bam
run_test --use_secondary

### Test 8 ###
TEST_FILE=$data_dir/all_conds.bam
run_test --percentage_aln 0.9 --percentage_id 0.9

### Clean up
rmdir $output_dir
