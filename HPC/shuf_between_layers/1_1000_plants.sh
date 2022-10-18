for i in {1..1000}
do
qsub i_plants.sh $i
echo $i
done