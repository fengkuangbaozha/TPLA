del=$1
file=$2


grep -Fvf $del ${file}.lncname > ${file}.lncname.iuox.del    ####delete TE and hkRNA
grep -Fvf $del ${file}.gtf > ${file}.del.gtf
#grep -Fvf $del ${file}.ej > ${file}.ej.del
cut -f2 ${file}.lncname.iuox.del |sort -u > ${file}.lncname.iuox.del.name
