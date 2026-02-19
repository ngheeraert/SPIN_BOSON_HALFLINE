npadd=10
delta=0.2
alpha=0.4
tmax=8
merr=0.0000003
dt=0.001
verbose=1
p0=0.0000001


newjobfile="quench"
cp jobfile $newjobfile

echo "./mpol -npini 1 -npadd $npadd -del $delta -al $alpha -wmax 5 -wc 1 -tmax $tmax -dt $dt -merr $merr -p0 $p0 -verbose $verbose" #  >> $newjobfile

echo "mv ../MPOL_\$tpdir/data/* \$PBS_O_WORKDIR/data/" >> $newjobfile
echo "rm -r ../MPOL_\$tpdir" >> $newjobfile

#qsub $newjobfile
rm $newjobfile

