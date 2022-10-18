for i in {1..1000}
do
qsub i_both.sh $i
echo $i
done