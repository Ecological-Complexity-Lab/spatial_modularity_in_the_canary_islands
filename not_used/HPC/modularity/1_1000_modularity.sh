for i in {1..1000}
do
qsub i_modularity.sh $i
echo $i
done