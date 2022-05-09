all:
	python3 compute_stability.py
	python3 initial_guess.py
	mpiexec -n 10 python3 rayleigh-benard.py
	defcon gui -p rayleigh-benard.py

clean:
	rm -rf output
	rm -rf initial_guess
	rm -rf tmp

gui:
	defcon gui -p rayleigh-benard.py

branches=0 166 42

stability:
	python3 save_functional.py --branchids $(branches)
	python3 compute_stability_.py --branchids $(branches)
	#cp -r pcitrues

my_vec=1 2 3
test:
	for xx in $(my_vec); do echo $$xx; done

create_latex:
	cp diagram_u/* Tex_RB/data/diagram_u/
	cp diagram_T/* Tex_RB/data/diagram_T/
	cp diagram_B/* Tex_RB/data/diagram_B/
	cp StabilityFigures/* Tex_RB/data/StabilityFigures/
	cd Tex_RB; for branchid in $(branches); do cp diagram_Fabian.tex diagram_$$branchid.tex; sed -i '' "s/branchid{}/branchid{$$branchid}/" diagram_$$branchid.tex; pdflatex -shell-escape diagram_$$branchid; done; cd ..

recompile_latex:
	cd Tex_RB; for branchid in $(branches); do pdflatex -shell-escape diagram_$$branchid; done; cd ..
