
2rvd.data: 2rvd_autopsf.pdb 2rvd_autopsf.psf par_all36_prot.prm chignolin.py in.lammps
	export PYTHONPATH=../../../../package ; python3 chignolin.py

	@lmp -in in.lammps

clean:
	@rm log.lammps 2rvd.data *.jpg
