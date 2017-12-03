%% Initialize
%profile clear
%profile on

clearvars

r = 2;
numProbs = 1;    % number of problems
filename = 'table_noisy.tex';

opts.verbose      = 1;
opts.minsize      = r+1;
opts.maxsize      = 5*(r+1)+15;

% Small problems, decreasing PP
nn = [2000; 2000;];
mm = [1000; 1000;];
NF = [0.01; 0.02;]; 
PP = [0.35; 0.35;]; % decreasing the density

%% Run tests
run_testsnoisy(r, nn, mm, NF, PP, numProbs, filename, opts);
%profile report
