for i in {1..1000}
do
qsub i_plants_del.sh $i
echo $i
done