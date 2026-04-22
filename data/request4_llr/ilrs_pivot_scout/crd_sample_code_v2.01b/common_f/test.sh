#!/bin/sh
# Test wr_test and crd_readf/crd_write modules
cd test_files
../wr_testf 7080_giovea_crd_20080508_09_00.frd 7080_giovea_crd_20080508_09_00.frd_wr
../wr_testf 7080_giovea_crd_20080508_09_00.npt 7080_giovea_crd_20080508_09_00.npt_wr
../wr_testf 7080_giovea_crd_20080508_09_00_v2.frd 7080_giovea_crd_20080508_09_00_v2.frd_wr
../wr_testf 7080_giovea_crd_20080508_09_00_v2na.frd 7080_giovea_crd_20080508_09_00_v2na.frd_wr
../wr_testf 7080_giovea_crd_20080508_09_00_v2.npt 7080_giovea_crd_20080508_09_00_v2.npt_wr
../wr_testf 7080_giovea_crd_20080508_09_00_v2na.npt 7080_giovea_crd_20080508_09_00_v2na.npt_wr
echo Comparing full rate files
echo ------ V1
diff 7080_giovea_crd_20080508_09_00.frd_wr 7080_giovea_crd_20080508_09_00.frd_wr.ref
echo ------ V2
diff 7080_giovea_crd_20080508_09_00_v2.frd_wr 7080_giovea_crd_20080508_09_00_v2.frd_wr.ref
echo ------ V2na
diff 7080_giovea_crd_20080508_09_00_v2na.frd_wr 7080_giovea_crd_20080508_09_00_v2na.frd_wr.ref
echo comparing normal point files
echo ------ V1
diff 7080_giovea_crd_20080508_09_00.npt_wr 7080_giovea_crd_20080508_09_00.npt_wr.ref
echo ------ V2
diff 7080_giovea_crd_20080508_09_00_v2.npt_wr 7080_giovea_crd_20080508_09_00_v2.npt_wr.ref
echo ------ V2na
diff 7080_giovea_crd_20080508_09_00_v2na.npt_wr 7080_giovea_crd_20080508_09_00_v2na.npt_wr.ref

