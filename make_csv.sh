#!/bin/zsh

iterC=10
echo -n "iterC,"
echo $2
for i in $(seq $iterC)
do
	echo -n "$i,"
	eval $1
done
