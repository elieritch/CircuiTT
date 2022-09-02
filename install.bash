#!/bin/bash

root_dir=$(dirname "$0");
cd ${root_dir}
wget https://figshare.com/ndownloader/files/34540430 -O allfiles.tar.gz;
tar -xzvf allfiles.tar.gz;

mv -v ${root_dir}/allfiles/supp_files/* ${root_dir}/supp_files/;
mv -v ${root_dir}/allfiles/test_data_input/* ${root_dir}/test_data_input/;
mv -v ${root_dir}/allfiles/test_data_output/* ${root_dir}/test_data_output/;
mv -v ${root_dir}/allfiles/picard.jar ${root_dir}/dependencies/picard.jar;

rm allfiles.tar.gz;
rm -r allfiles;