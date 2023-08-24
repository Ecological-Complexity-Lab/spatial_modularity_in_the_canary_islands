for i in {1..1000}
do
qsub i_interactions.sh $i
echo $i
done