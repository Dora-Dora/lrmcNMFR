function  run_testsnoisy(r, nn, mm, NF, PP, numProbs, filename, opts)

num_fullcomplete = zeros(length(nn), numProbs);
rv               = zeros(length(nn), numProbs);
percentrecover   = zeros(length(nn), numProbs);
rank_init        = zeros(length(nn), numProbs);
rank_refined     = zeros(length(nn), numProbs);
resid_init       = zeros(length(nn), numProbs);
resid_refined    = zeros(length(nn), numProbs);
prerefine_time   = zeros(length(nn), numProbs);
refined_time     = zeros(length(nn), numProbs);
densities        = 0*PP;   % for true densities from generation

for ll = 1:length(nn)
    
    params.m    = mm(ll);
    m=params.m;
    params.n    = nn(ll);
    n=params.n;
    params.r    = r; % desired rank for incomplete Z and for complete Y/target rank
    params.tolerrank = max(m,n)*eps+1e-2*NF(ll); % change below after Zpart????
    params.noiselevel   = NF(ll); % noise level related to the magnitude of the entries
    params.p    = PP(ll);   % percent known elements in Z
    
   %%%%% take average over numProbs %%%%%
    for pp=1:numProbs 
        
        params.seed = pp;
       
        if opts.verbose
            disp('-------------------------------------------------');
            disp(params);
        end
        
       % generate random problem and solve
       [Zorig,Zpart,indsZ,~]= Zgenerator(params);
       params.tolerrank = max(m,n)*eps(norm(Zorig))+1e-2*NF(ll); % changed here
       tolerrank=params.tolerrank;
       densities(ll,pp)=length(indsZ)/(numel(Zorig));
       problem.Zorig=Zorig;
       problem.Zpart=Zpart;
       problem.indsZ=indsZ;

       % solve the problem
       [flag,resid1,resid2,rank1,rank2,time1,time2,nv,ps] = ...
            CompleteZ(m,n,r,problem,opts,tolerrank,true);
        if flag == 0    % if it is properly solved
            num_fullcomplete(ll,pp) = 1;
        else
            num_fullcomplete(ll,pp) = 0;
        end
        rv(ll,pp)               = nv;
        percentrecover(ll,pp)   = ps;
        rank_init(ll,pp)        = rank1;
        rank_refined(ll,pp)     = rank2;
        resid_init(ll,pp)       = resid1;
        resid_refined(ll,pp)    = resid2;
        prerefine_time(ll,pp)   = time1;
        refined_time(ll,pp)     = time2;
    end   

fmt1 = '%5.0f & %5.0f & %5.2f & %5.2f &';

% if allLocalized
%     fmt2 = '';
% else
%     fmt2 = '%9.1f & ';
% end

fmt3 = '%7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2e & %7.2e';
fmt4 = ' \\cr\\hline\n';
fmt5 = '\\cr\\cline{1-4}\\cline{6-11}\n';

hfmt1 = regexprep(fmt1, '(\.[0-9])*[def]', 's');
%hfmt2 = regexprep(fmt2, '(\.[0-9])*[def]', 's');
hfmt3 = regexprep(fmt3, '(\.[0-9])*[def]', 's');

fid = fopen(filename, 'w');

fprintf(fid, '%s\n', '\begin{tabular}{|cccc||c|cc|cc|cc|} \hline');
fprintf(fid, '%s', [...
    '\multicolumn{4}{|c||}{Specifications} & ', ...
    '\multirow{2}{*}{Rcvd (\%$Z$)} & ', ...
    '\multicolumn{2}{|c|}{Time (s)} & ', ...
    '\multicolumn{2}{|c|}{Rank} &', ... 
    '\multicolumn{2}{|c|}{Residual (\%$Z$)}']);
fprintf(fid, fmt5);
fprintf(fid, hfmt1, '$m$', '$n$', '\% noise', '$p$');
%fprintf(fid, hfmt2, 'Localized');
fprintf(fid, hfmt3, '', ...
    'initial', 'refine', ...
    'initial', 'refine', ...
    'initial', 'refine');
fprintf(fid, fmt4);

for jj = 1:ll
  indsucc=find(prerefine_time(jj,:)>0);
  fprintf(fid, fmt1, mm(jj), nn(jj), 100*NF(jj), mean(densities(jj,indsucc)));
  fprintf(fid, fmt3, 100*mean(percentrecover(jj,indsucc)), ...
    mean(prerefine_time(jj,indsucc)), mean(refined_time(jj,indsucc)), ...
      mean(rank_init(jj,indsucc)), mean(rank_refined(jj,indsucc)), ...
        mean(resid_init(jj,indsucc)), mean(resid_refined(jj,indsucc)));
  fprintf(fid, fmt4);
end
fprintf(fid, '\\end{tabular}\n');

fclose(fid);

system(['cat ', filename]);

end   % end of for ll

end
