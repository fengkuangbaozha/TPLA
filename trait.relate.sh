trait=$1
cor=$2     #floral.correlation0.95.txt
gtfname=$3    #cuff85-CPC_left.lncname
lncall=$4
name=$5


grep -Fwf $trait $cor |cut -f 1 |sort -u > ${name}.coexpressed.lncrna   ####get the ripen genes coexpressed lncRNAs
grep -Fwf $trait $gtfname |cut -f 2 |sort -u > ${name}.nearest.lncrna     #####get the ripen genes nearest lncRNAs
comm -1 -2 ${name}.nearest.lncrna ${name}.coexpressed.lncrna |sort > ${name}.coexpressed.nearest.lncrna   ###find ripen gene correlate lnRNAs
grep -Fwf ${name}.coexpressed.nearest.lncrna $lncall > ${name}.coexpressed.nearest.lncrna.info ######extract the full information 

