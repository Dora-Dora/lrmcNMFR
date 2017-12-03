%% Initialize

clearvars
%profile clear
%profile on

r = 3;
numProbs = 1;    % number of problems
filename = 'table_noiseless.tex';

opts.verbose      = false;
opts.minsize      = r+1;
opts.maxsize      = 5*(r+1)+10;

% Small problems, decreasing PP
nn = [4000; 4000;];
mm = [2100; 2100;];
NF = [0.00; 0.00;]; 
PP = [0.40; 0.35;]; % decreasing the density

%% Run tests
run_testsnoiseless(r, nn, mm, [], PP, numProbs, filename, opts);
%profile report
