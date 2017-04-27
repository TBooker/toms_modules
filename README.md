These modules contain functions that I use for lots of different scripts

Just so that you know in the future. To add files to the python path you need to make or open a file called .profil in the default directory that terminal or console or whatever Unix commmand line interface you are using is called.

Add the following lines to your .profile to add custom modules to your python path (so they can be called from within python scripts easily):
```sh
PYTHONPATH="path/to/dir/"
export PYTHONPATH
```

This adds a directory to the python path so that you can access custom modules within your scripts.

This is useful as it means that you can store all your functions in one place and you don't have to re-define them for each scripts. Also, since a .py module is compiled before it is used (generating a .pyc file) it speeds things up a bit.


The modules:

  tom.py has a bunch of commands that I use from time to time that make life a little bit easier

  tom_slim.py has files that can be used to handle SLiM output files in the way that I store them (as a single GZIPPED text file that has teh word "None" on a single line separting the individual runs.

  site_frequency_spectrum.py has functions for doing things with the site frequency spectrum (calculating pi, tajima's D, blah blah)

  tom_vcf.py several filters that can be used to fiddle about with VCF lines

  superFreq.py has functions for handling the 2-outgroup .freq format
  
