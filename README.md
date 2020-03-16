# GEM-GPU
The GPU version for GEM is undergoing.

to run gem: first 'make clean' and 'make' under the dfftpack directory; then 'make clean' and 'make' in code; finnaly copy 'gem_main' to 'test-summit' and submit the job 'job-summit.h'.

the parameters are set in gem.in, note that the ntube should be the same as that of gem_comm.f90.

