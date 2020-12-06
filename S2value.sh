for l in {1..2}  #crit
  do for k in {1..2}  #size
    do for j in {1..2}  #propensity
      do for i in {1..4} #beta
        do sbatch -t 24:00:00 --mem=10000 S0run.sh C21.simulation_run.R $i $j $k $l $1  #$1 = date
      done
    done
  done
done