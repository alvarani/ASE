#perl addheader.pl header *.sam /bubo/home/h24/alvaj/glob/code/ASE/ASE/GSNAP/bam/header

#!/bin/perl
use strict;
use warnings;
open (FILE,  '<',"header") or die $!;
open (FILE2,  '<', "2_LPS.part_2ab.uniqmappedN.sam") or die $!;
while ($_=<FILE>) {
print $_ ;
}
close (FILE);
while ($_=<FILE2>) {
print $_ ;
}
close (FILE2)



