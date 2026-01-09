### In here we prepare the input files for 1D and 2D US.

### First we need:
- md_eq7.mdp
- md_new_nucl2.gro (from previous step in prep dir)
- plumed_eq.dat (modified based on the correct atoms)
- index.ndx (from previous steps)
- and finally job_eq7.sh for launching the job.

In the end, we will have a md_eq7.gro that will be used in 1D and 2D folders.

In each folder 1D and 2D, there are bash scripts "prepare_file_std_us_1D" and "prepare_file_std_us_2D" to prepare the input files.
Make sure you modify the begining of these script to have correct atoms numbers.
Additionally, due to the limit number of active jobs on Adastra (300 max), there are two scripts: steered_BatchSubmit.sh and us_BatchSubmit.sh that help launching the jobs in batches and keep tracking them on the queue untill all the jobs are finished.
Make sure you run BatchSubmit files on tmux!
