function [Zorig,Zpart,inds,noiseN] = Zgenerator(params)
% in:
%   params = a struct containing m, n, r, p (see below)
% out:
%   Zpart = a m by n partial matrix with approx percentage p elements known
%   Zorig = the full original random generated matrix of rank r
%   inds = indices according to know entries of Zpart
%   noiseN = true noise level

m = params.m; % number of rows
n = params.n; % number of cols
r = params.r; % rank of Zorig, true rank
p = params.p; % percentage of known entries

if isfield(params,'noiselevel') 
    noiselevel = params.noiselevel;
else
    noiselevel = 0;
end
seed = params.seed;

rng(seed);
saverng = rng; % for redo the same experiments
save('currentrng.mat','saverng')

P = randn(m,r);
Q = randn(r,n);
Zorig = P*Q; % Zorig should be rank r

S = sprand(m,n,p);
inds = find(S);
Zpart = spalloc(m,n,length(inds));

if noiselevel > 0
    magnitude = norm(Zorig(inds),1)/length(inds);
	Noise = noiselevel*magnitude*randn(length(inds),1);
	Zpart(inds) = Zorig(inds) + Noise;
else
	Zpart(inds) = Zorig(inds);
	Noise = 0;
end
noiseN = norm(Noise,1);

end  % end of function Zgenerator

