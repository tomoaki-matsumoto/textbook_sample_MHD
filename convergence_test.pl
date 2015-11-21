#! /usr/bin/perl -w
use strict;

## Hydrodynamical wave
# my $logfile = 'NoneEuler.log';
# my $logfile = 'MUSCL2PC2.log';
# my $logfile = 'MUSCL2RK2.log';
# my $logfile = 'MUSCL3PC2.log';
# my $logfile = 'MUSCL3RK3.log';
# my $logfile = 'MUSCL3K3.log';	# Kutta 3rd order

## scalar advection
# my $logfile = 'AdvNoneEuler.log';
# my $logfile = 'AdvMUSCL2Euler.log';
# my $logfile = 'AdvMUSCL2PC2.log';
# my $logfile = 'AdvSuperBee2PC2.log';
# my $logfile = 'AdvMUSCL2RK2.log';
# my $logfile = 'AdvMUSCL3PC2.log';
# my $logfile = 'AdvMUSCL3RK3.log';
# my $logfile = 'AdvMUSCL3K3.log';


my $logfile = 'test.log';

my $logdir = 'convergence_test';
my $devnull = '>/dev/null 2>&1';
# my $devnull = '';

my @res_pw = ( 5 .. 12 );  # hydro wave 精度が向上すると波が突っ立つので格子点数を 2^12 までとした。
# my @res_pw = ( 5 .. 13 );  # scalar advection
# my @res_pw = ( 8 );
unlink "$logdir/$logfile";
print "$logdir/$logfile\n";
foreach my $p (@res_pw) {
    my $nx = 2 ** $p;
    my $cppflags = "\"CPPFLAGS = -DCONVERGENCE_TEST -DNX=$nx -DNY=1 -DNZ=1 -DNDIM=1\"";
    print $nx, " ";
    system("make clean $devnull");
    system("make $cppflags $devnull");
    system("make $cppflags errornorm $devnull");
    system("./main $devnull");
    system("./errornorm | tail -1 >> $logdir/$logfile");
}
print "\n";

