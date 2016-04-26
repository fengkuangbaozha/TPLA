dir=$1
name=$2


#####RNA structure information##########
mkdir $dir/database/struc_eps
mkdir $dir/database/struc_png
for i in `cut -f1 $dir/database/txt/*.lncrna.oneline.fa |cut -d_ -f1,2`; do cp $dir/stru/${i}_ss.ps $dir/database/struc_eps; cp $dir/stru/${i}_ss-1.png $dir/database/struc_png; done
for i in $dir/database/struc_png/*; do rename '-1' '' $i; done
for i in $dir/database/struc_*/*; do rename 'TCONS' "${name}_TCONS" $i; done
cd $dir/database
#tar czf struc_eps.tar.gz $dir/database/struc_eps/
#tar czf struc_png.tar.gz $dir/database/struc_png/
for i in $dir/database/txt/*; do sed "s/TCONS_/${name}_TCONS_/g" $i |sed "s/XLOC_/${name}_XLOC_/g" > ${i}.new; mv ${i}.new ${i}; done
rm */*.new


