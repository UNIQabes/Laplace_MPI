#!/bin/zsh

iterC=10
sum=0
for i in $(seq $iterC)
do
	a=$(eval $1)
	sum=$((a+sum))
done
sum=$((sum/iterC))
echo $sum