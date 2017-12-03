function  run_testsnoiseless(r, nn, mm, NF, PP, numProbs, filename, opts)

num_fullcomplete = zeros(length(nn), numProbs);
rank_init        = zeros(length(nn), numProbs);
resid_init       = zeros(length(nn), numProbs);
prerefine_time   = zeros(length(nn), numProbs);
densities        = zeros(length(nn), numProbs);
percentsrecoverd = zeros(length(nn), numProbs);
rv               = zeros(length(nn), numProbs);
rcvrd            = zeros(length(nn), numProbs);
if ~isempty(NF)
    resid_refined    = zeros(length(nn), numProbs);
    rank_refined     = zeros(length(nn), numProbs);
    refined_time     = zeros(length(nn), numProbs);
end

for ll = 1:length(nn)
    
    params.m    = mm(ll);
    m           = params.m;
    params.n    = nn(ll);
    n           = params.n;
    params.r    = r;
    params.tolerrank = max(m,n)*eps;
    tolerrank = params.tolerrank;
    if ~isempty(NF)
        params.noiselevel = NF(ll); % noise level related to the magnitude of the entries
    end
    params.p    = PP(ll); 
    
    %%%%% take average over numProbs %%%%%
    for pp=1:numProbs
        
        params.seed = pp; % parameter for rng(seed)
       
        if opts.verbose
            disp('-------------------------------------------------');
            disp(params);
        end
    
       % generate random problem and solve
       [Zorig,Zpart,indsZ,realNoise]= Zgenerator(params);
       problem.Zorig=Zorig;
       problem.Zpart=Zpart;
       problem.indsZ=indsZ;
       problem.realNoise=realNoise;
       
       % solve the problem
       [flag,resid1,resid2,rank1,rank2,time1,time2,nv,ps] = ...
           CompleteZ(m,n,r,problem,opts,tolerrank,false);
       densities(ll,pp) = length(indsZ)/(numel(Zorig));
       percentsrecoverd(ll,pp)  = ps; 
       fprintf('[ll pp m n rn]  %g %g  %g %g %g \n',ll,pp,m,n,nv)
       fprintf('%g  percent recovered  \n',ps)

        if flag == 0
            num_fullcomplete(ll,pp) = 1;
        else
            num_fullcomplete(ll,pp) = 0;
        end
        rank_init(ll,pp)        = rank1;
        resid_init(ll,pp)       = resid1; 
        prerefine_time(ll,pp)   = time1;
        rv(ll,pp)               = nv;
        rcvrd(ll,pp)            = ps;
        if ~isempty(NF)
            resid_refined(ll,pp)    = resid2;
            rank_refined(ll,pp)     = rank2;
            refined_time(ll,pp)     = time2;
        end
    end   
    
% Latex table formatting
fmt1 = '%5.0f & %5.0f & %5.2f &';

% if allLocalized
%     fmt2 = '';
% else
%     fmt2 = '%9.1f & ';
% end

fmt3 = '%7.2f & %7.2f&  %7.2f & %7.1f & %9.4e';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-3}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|ccc||c|c|c|c|c|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{3}{|c||}{Specifications} & ', ...
    '\multirow{2}{*}{$r_v$} & ', ...
    '\multirow{2}{*}{Rcvrd (\%$Z$)} & ', ...
    '\multirow{2}{*}{Time (s)} & ', ...
    '\multirow{2}{*}{Rank} &', ... 
    '\multirow{2}{*}{Residual (\%$Z$)}']);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$m$', '$n$', 'mean($p$)');
fprintf(fid, hfmt3, '','','','','');
fprintf(fid, fmt4);
for jj = 1:ll   % changed to ll from length nn to save table each case
    %fprintf(fid, fmt1, mm(jj), nn(jj), PP(jj));
    fprintf(fid, fmt1, mm(jj), nn(jj), mean(densities(jj,:)));
    %fprintf(fid, fmt2, mean(num_fullcomplete(jj,:)));
    fprintf(fid, fmt3, ...
        mean(rv(jj,:)), ...
        mean(rcvrd(jj,:)), ...
        mean(prerefine_time(jj,:)), ...
        mean(rank_init(jj,:)), ...
        mean(resid_init(jj,:)));
    fprintf(fid, fmt4);
end

fprintf(fid, '\\end{tabular}\n');
fclose(fid);

system(['cat ', filename]);
system(['cp ', filename, ' ..\..\..\latexfiles.d' ]);

end    % end loop over 1:nn

end    % function end
