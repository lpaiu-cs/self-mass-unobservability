cd test_files
../merit_to_crd -i s25y07d336t1600#1155.mrt -o 7080_lageos1_20071202_16_00.frd
../crd_to_merit -i 7080_lageos1_20071202_16_00.frd -o s25y07d336t1600#1155.mrt_c
../cstg_to_crd -i s25y07d336t1600#1155.qld -o 7080_lageos1_20071202_16_00
../cstg_to_crd -i s25y07d336t1600#1155.npt -o 7080_lageos1_20071202_16_00
../crd_to_cstg_np -i 7080_lageos1_20071202_16_00.npt -o s25y07d336t1600#1155.np
../crd_to_cstg_np -i 7080_lageos1_20071202_16_00.np2 -o s25y07d336t1600#1155.np-2
../crd_to_cstg_np -i 7080_lageos1_20071202_16_00.np2-na -o s25y07d336t1600#1155.np-2na
../crd_to_cstg_ql -i 7080_lageos1_20071202_16_00.qlk -o s25y07d336t1600#1155.ql
../crd_to_cstg_ql -i 7080_lageos1_20071202_16_00.ql2-na -o s25y07d336t1600#1155.ql-2na
echo comparing 7080_lageos1_20071202_16_00.frd
diff 7080_lageos1_20071202_16_00.frd 7080_lageos1_20071202_16_00.frd.ref
echo comparing s25y07d336t1600#1155.mrt_c
diff s25y07d336t1600#1155.mrt_c s25y07d336t1600#1155.mrt_c.ref
echo comparing 7080_lageos1_20071202_16_00.npt
diff 7080_lageos1_20071202_16_00.npt 7080_lageos1_20071202_16_00.npt.ref
echo comparing 7080_lageos1_20071202_16_00.qlk
diff 7080_lageos1_20071202_16_00.qlk 7080_lageos1_20071202_16_00.qlk.ref
echo comparing s25y07d336t1600#1155.np
diff s25y07d336t1600#1155.np s25y07d336t1600#1155.np.ref
echo comparing s25y07d336t1600#1155.np CRD v2
diff s25y07d336t1600#1155.np-2 s25y07d336t1600#1155.np-2.ref
echo comparing s25y07d336t1600#1155.np CRD v2-na
diff s25y07d336t1600#1155.np-2na s25y07d336t1600#1155.np-2na.ref
echo comparing s25y07d336t1600#1155.ql
diff s25y07d336t1600#1155.ql s25y07d336t1600#1155.ql.ref
echo comparing s25y07d336t1600#1155.ql CRD v2-na
diff s25y07d336t1600#1155.ql-2na s25y07d336t1600#1155.ql2-na.ref
