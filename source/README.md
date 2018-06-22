* outfiles/ contains figures as well as output files utilized by other code files. This is also where output files will go when an experiment is run.

* infiles/ contains a directory networkfiles/ that contains the thresholded brain-specific functional interaction networks. It also contains positive and negative sets.

* experiments/ contains bash files for running experiments. A bash file can be run in the command line with the following command:

$ bash filename.bash

* code/ contains our method (predict2.py, fileIO2.py, learners2.py) as well as code to generate figures for the paper (visualize.py, chartmaker.py)