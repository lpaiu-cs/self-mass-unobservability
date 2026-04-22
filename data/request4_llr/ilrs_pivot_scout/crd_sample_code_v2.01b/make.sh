#!/bin/sh
# Make crd sample code 
# To "make all", type
#	./make.sh
# To "make clean", type 
#	./make.sh clean		# "clean" or anything but blanks
# N.B. make common directories first
cd common_c
make clean
if test -z $1
then
  make all
fi
cd ..
cd crd_conv_slr_c
make clean
if test -z $1
then
  make all
fi
cd ..
# Now go through all directories...
for i in `ls`
  do
    if test -d $i
    then
      cd $i
      if test -e NO_AUTO_COMPILE
      then
        echo "Bypassing $i"
      else
        echo "Doing $i"
	make clean
	if test -z $1
	then
          make all
	fi
      fi
      cd ..
    fi
  done
echo done
