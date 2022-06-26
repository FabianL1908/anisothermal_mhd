solver_type=fs3by2
discr=bdm
hierarchy=uniform
dt=0.01
Tf=0.1
k=2
output= --output
Ncores=48

ifeq ($(dim),2)
	Ncores=64
	baseN=16
	nref=6
else
ifeq ($(dim),3)
	Ncores=40
	baseN=8
	nref=3
endif
endif

exec=srun --ntasks-per-node $(Ncores) ${VIRTUAL_ENV}/bin/python

ifeq ($(mode),stat)
	Ras2=10 100 300 1000 3000 10000 30000 100000 200000 500000 1000000
	Prs1=1 0.1 0.01 0.003 0.001
	Prs2=1.0 0.1 0.01 0.001
endif

ifeq ($(mode),time)
	Ras=1 100000 1000000
	Prs=1 0.1 0.01 0.001
endif

file_stat=mhd$(dim)d_stationary_upT_BE.py
file_time=mhd$(dim)d_timedep_upT_BE.py

pre:
	mkdir -p logs/
	mkdir -p results/
	mkdir -p dump/

clean_error:
	rm -rf error.txt

allstat: newtonstat picardstat

alltime: newtontime picardtime

newtonstat: pre
	for pr in $(Prs1); do $(exec) $(file_stat) --baseN $(baseN) --k $(k) --nref $(nref) --discr $(discr) --Pr $$pr --Pm 1 --S 1 --Ra 1 --gamma 1000 --gamma2 0 --hierarchy $(hierarchy) --solver-type $(solver_type) --testproblem $(testproblem) --linearisation newton --stab $(output) --checkpoint2 2>&1 | tee -a logs/ldcnewton.log; done
	for pr in $(Prs2);  do $(exec) $(file_stat) --baseN $(baseN) --k $(k) --nref $(nref) --discr $(discr) --Pr $$pr --Pm 1 --S 1 --Ra $(Ras2) --gamma 1000 --gamma2 0 --hierarchy $(hierarchy) --solver-type $(solver_type) --testproblem $(testproblem) --linearisation newton --stab $(output) --checkpoint  2>&1 | tee -a logs/ldcnewton.log; done

picardstat: pre
	$(exec) $(file_stat) --baseN $(baseN) --k $(k) --nref $(nref) --discr $(discr) --Re 1 --Rem $(Rems) --S $(S) --gamma 1000 --gamma2 0 --hierarchy $(hierarchy) --solver-type $(solver_type) --testproblem $(testproblem) --linearisation picard --stab $(output) 2>&1 | tee -a logs/ldcpicard.log
	$(exec) $(file_stat) --baseN $(baseN) --k $(k) --nref $(nref) --discr $(discr) --Re $(Res) --Rem $(Rems) --S $(S) --gamma 1000 --gamma2 0 --hierarchy $(hierarchy) --solver-type $(solver_type) --testproblem $(testproblem) --linearisation picard --stab $(output) --checkpoint 2>&1 | tee -a logs/ldcpicard.log

newtontime: pre
	$(exec) $(file_time) --baseN $(baseN) --k $(k) --nref $(nref) --discr $(discr) --Re $(Res) --Rem $(Rems) --S $(S) --gamma 1000 --gamma2 0 --hierarchy $(hierarchy) --solver-type $(solver_type) --testproblem $(testproblem) --linearisation newton --stab $(output) --dt $(dt) --Tf $(Tf) | tee -a logs/ldcnewton.log

picardtime: pre
	$(exec) $(file_time) --baseN $(baseN) --k $(k) --nref $(nref) --discr $(discr) --Re $(Res) --Rem $(Rems) --S $(S) --gamma 1000 --gamma2 0 --hierarchy $(hierarchy) --solver-type $(solver_type) --testproblem $(testproblem) --linearisation picard --stab $(output) --dt $(dt) --Tf $(Tf) | tee -a logs/ldcpicard.log

test: pre
	mpiexec -n 2 python mhd2d_stationary.py --baseN 10 --k 2 --nref 1 --discr bdm --Re 1 --Rem 1 --S 1 --gamma 1000 --gamma2 0 --hierarchy uniform --solver-type fs2by2 --testproblem ldc --linearisation newton --stab

printresults:
	python3 print_results.py --testproblem $(testproblem)
