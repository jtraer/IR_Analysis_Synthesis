function SummarizeCode(cfl)
%* Finds all commented lines within a file which begin with an asterisk (i.e. '%*') and write them to a .org file.  This file (best viewed with emacs) hopefully provides a summary of the steps followed in the code.

eval(sprintf('! grep "%%\\*" %s.m > tmp.org',cfl))
eval('! sed ''s/^ *//g'' < tmp.org > tmp2.org');  % remove whitespace
eval('! sed ''s/^[ \\t]+//g'' < tmp2.org > tmp.org');  % remove tabs too
eval(sprintf('! sed ''s/^.//'' tmp.org > %s.org',cfl))   % remove '%' so emacs can read the indenting
