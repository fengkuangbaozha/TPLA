dir=$1
name=$2

##########RNA structure ################
mkdir -p $dir/stru
cd $dir/stru
/psc/program/install/ViennaRNA/bin/RNAfold --MEA -d2 -p < $dir/${name}-CPC_left.lncname.iuox.del.fa > $dir/stru/${name}-CPC_left.lncname.iuox.del.fold
for i in $dir/stru/*_ss.ps; do j=`echo $i |rev |cut -d/ -f1 |rev |sed 's/ps/pdf/g'`; z=`echo $i |rev |cut -d/ -f1 |rev |sed 's/\.ps//g'`; ps2pdf $dir/stru/$i $dir/stru/$j; pdftoppm -png $dir/stru/$j $dir/stru/$z; done

