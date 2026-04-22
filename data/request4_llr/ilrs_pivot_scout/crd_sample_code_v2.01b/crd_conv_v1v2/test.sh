#!/bin/bash
#
# Test crd_conv_v1v2 programs
#============================
cd test_data
for i in `ls s25*npt s25*frd`
do
echo Running crd_conv on $i
  ../crd_conv $i $i.v2
  echo Compare $i.v2
  diff $i.v2 $i.v2.ref
done
for i in `ls lag*npt lag*frd`
do
echo Running crd_create on $i
  ../crd_create $i $i.v2Y ../create.dat.YARL
  echo Compare $i.v2Y
  diff $i.v2Y $i.v2Y.ref
done
for i in `ls lag*npt lag*frd`
do
echo Running crd_create on $i
  ../crd_create $i $i.v2C ../create.dat.CHAL
  echo Compare $i.v2C
  diff $i.v2C $i.v2C.ref
done
for i in `ls *.v2Y *.v2C`
do
  ../../crd_chk_c/crd_chk $i > $i.chk
  echo Compare $i.chk
  diff $i.chk $i.chk.ref
done
