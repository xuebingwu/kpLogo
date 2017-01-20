	for k in 3
	do 
	    for shift in 1
	    do
	    echo $k, $shift
	    rm log output.* data.*; cross_validation -data KO.txt2  -train_cmd "../../bin/PKA data.train.@ -minCount 0.01 -o output.@ -max_k $k -max_shift $shift -weighted 2>> log " -test_cmd "../../bin/PKA data.test.@ -predict output.@ -weighted 2>> log" -nfold 5
		rm tmp
	    cat output.*.score | grep -v inf > output.score
	    Rscript cor.r
	    done 
	done
