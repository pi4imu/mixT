f042=$1
ab=$2
t=$3
z=$4
nh=$5
#
cat ttt.dat | LC_NUMERIC="C" awk '{if (NR>=4) print $1,$5,$5*$1}' > tttx.dat
dmstat tttx.dat > ttts.dat
#cat ttts.dat
nch=$(cat ttts.dat | grep good | awk '{if (NR==1) print $2}')
es=$(cat ttts.dat | grep sum | awk '{if (NR==1) print $2}')
e1=$(cat ttts.dat | grep min | awk '{if (NR==1) print $2}')
e2=$(cat ttts.dat | grep max | awk '{if (NR==1) print $2}')
cs=$(cat ttts.dat | grep sum | awk '{if (NR==2) print $2}')
ecs=$(cat ttts.dat | grep sum | awk '{if (NR==3) print $2}')
echo $f042 $ab $t $z $nh $nch $e1 $e2 $es $cs $ecs
rm ttt*.dat
#rm ttts.dat
