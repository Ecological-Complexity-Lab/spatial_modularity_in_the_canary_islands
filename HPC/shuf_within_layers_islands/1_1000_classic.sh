for i in {1..1000}
do
qsub i_classic.sh $i
echo $i
done