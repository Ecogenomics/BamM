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

data_dir=~/git/BamM/test/filter_test_data
pre_bam=$data_dir/pre.bam
output_dir=$data_dir/output
bamm_exe=~/git/BamM/bin/bamm

### Helpers ###
output_file=$output_dir/$(basename $pre_bam .bam)_filtered.bam
test_1_file=$data_dir/none.bam
test_2_file=$data_dir/aln_only_0.5.bam
test_3_file=$data_dir/aln_only_0.9.bam
test_4_file=$data_dir/id_only_0.5.bam
test_5_file=$data_dir/id_only_0.9.bam
test_6_file=$data_dir/no_second_only.bam
test_7_file=$data_dir/no_supp_only.bam
test_8_file=$data_dir/all_conds.bam

run_test() {
  $bamm_exe filter -b $pre_bam -o $output_dir $@
}

assert_equal_length() {
  size1=samtools view $output_file |wc -l
  size2=samtools view $"test_$1_file" |wc -l
  if [ $size1 != $size2 ]; then echo "Failed test no. $1" >&2 ; exit 1; fi
}

### Test 1 ###
run_test --use_secondary --use_supplementary
assert_equal_length 1
rm $output_file

### Test 2 ###
run_test --use_secondary --use_supplementary --percentage_aln 0.5
assert_equal_length 2
rm $output_file

### Test 3 ###
run_test --use_secondary --use_supplementary --percentage_aln 0.9
assert_equal_length 3
rm $output_file

### Test 4 ###
run_test --use_secondary --use_supplementary --percentage_id 0.5
assert_equal_length 4
rm $output_file

### Test 5 ###
run_test --use_secondary --use_supplementary --percentage_id 0.9
assert_equal_length 5

### Test 6 ###
run_test --use_supplementary
assert_equal_length 6
rm $output_file

### Test 7 ###
run_test --use_secondary
assert_equal_length 7
rm $output_file

### Test 8 ###
run_test --percentage_aln 0.9 --percentage_id 0.9
assert_equal_length 8
rm $output_file
rmdir $output_dir

