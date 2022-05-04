all: clean
	python3 compute_stability.py
	python3 initial_guess.py
	mpiexec -n 12 python3 rayleigh-benard.py
	defcon gui -p rayleigh-benard.py

clean:
	rm -rf output
	rm -rf initial_guess
	rm -rf tmp

gui:
	defcon gui -p rayleigh-benard.py
