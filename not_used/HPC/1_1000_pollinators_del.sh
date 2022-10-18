for i in {1..1000}
do
qsub i_pollinator_del.sh $i
echo $i
done