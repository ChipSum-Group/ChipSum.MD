#!/bin/bash
if [[ ! -d "combine" ]]; then
mkdir -p combine
fi
# :>combine/sed/XDATCAR
# cp md1/POSCAR combine/sed/POSCAR
for i in {0..10}
do
cp Eckart.py Eckart_pos.py bondVASP_cell.py output_pos.py RDFVASP.py md$i/
cd md$i/
python3 output_pos.py
cd ../
done
cp Eckart.py Eckart_pos.py cal_VACF.py combine/
cp md1/POSCAR combine/
cd combine
for i in Ca Mg Al Si Fe Na O K
do
:>combine/$i.pos
cat md0/$i.pos>>combine/$i.pos
cat md1/$i.pos>>combine/$i.pos
cat md2/$i.pos>>combine/$i.pos
# cat ../md3/$i.pos>>$i.pos
# cat ../md4/$i.pos>>$i.pos
# cat ../md5/$i.pos>>$i.pos
# cat ../md6/$i.pos>>$i.pos
# cat ../md7/$i.pos>>$i.pos
# cat ../md8/$i.pos>>$i.pos
# cat ../md9/$i.pos>>$i.pos
# cat ../md10/$i.pos>>$i.pos
# tail -n +8 md$i/XDATCAR >>combine/sed/XDATCAR
done
# sed -e '/Direct configuration=/d' combine/sed/XDATCAR > combine/sed/XDATCAR0
# mv combine/sed/XDATCAR0 combine/sed/XDATCAR