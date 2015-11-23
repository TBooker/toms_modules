These modules contain functions taht I use for lots of differnt scripts

Just so that you know in the future. To add files to the python path you need to make or open a file called .profil in the default directory that terminal or console or whatever Unix commmand line interface you are using is called.

Write these lines:
PYTHONPATH="path/to/dir/"
export PYTHONPATH

^^ This adds a directory to the python path so that you can access custom modules within your scripts.

This is useful as it means that you can store all your functions in one place and you don't have to re-define them for each scripts. Also, since a .py module is compiled before it is used (generating a .pyc file) it speeds things up a bit.

\\\\\\\\\\\\\\\\

The module: tom.py contains some convenient functions taht I have written for several different things

tom_slim.py has files that can be used to handle slim output files in the way that I store them (as a single GZIPPED text file that has teh word "None" on a single line separting the individual runs.

site_frequency_spectrum.py has functions for doing things with the site frequency spectrum (calculating pi, tajima's D, blah blah)

email_tom.py
a function that allows you to send an email to notify when a certain process is complete

tom_vcf.py several filet that can be used to fiddle about with VCF lines
