#! /bin/sh
# ファイルを tar.gz に固めて redmagic に送る。
#
make distclean > /dev/null 2>&1 
SUBDIR=`basename $PWD`
TARBALL=$SUBDIR.tar.gz
flist="Makefile *.F90 *.f90 *.h *.sh *.pro *.gp"
mkdir -p $SUBDIR
cp -p $flist $SUBDIR
tar cvfz $TARBALL $SUBDIR
scp $TARBALL redmagic:public_html/konan15/
rm -rf $TARBALL $SUBDIR
