#! /usr/bin/perl -w
use strict;


# my $logfile = 'MUSCL3RK3.log';
# my $logfile = 'MUSCL3RK3splt.log';  # split
# my $logfile = 'MUSCL3K3.log';	# Kutta 3rd order
# my $logfile = 'MUSCL2RK2.log';
# my $logfile = 'MUSCL2PC2.log';
# my $logfile = 'NoneEuler.log';
# my $logfile = 'MUSCL3PC2.log';

## scalar advection
# my $logfile = 'AdvMUSCL3RK3.log';
# my $logfile = 'AdvMUSCL2RK2.log';
# my $logfile = 'AdvMUSCL3RK2.log';
# my $logfile = 'AdvNoneEuler.log';


my $logfile = 'test.log';

my $logdir = 'convergence_test';
my $devnull = '>/dev/null 2>&1';
# my $devnull = '';

my @res_pw = ( 5 .. 10 );
# my @res_pw = ( 8 );
unlink "$logdir/$logfile";
print "$logdir/$logfile\n";
foreach my $p (@res_pw) {
    my $nx = 2 ** $p;
    my $cppflags = "\"CPPFLAGS = -DCONVERGENCE_TEST -DNX=$nx -DNY=1 -DNZ=1\"";
    print $nx, " ";
    system("make clean $devnull");
    system("make $cppflags $devnull");
    system("make $cppflags errornorm $devnull");
    system("./main $devnull");
    system("./errornorm | tail -1 >> $logdir/$logfile");
}
print "\n";

