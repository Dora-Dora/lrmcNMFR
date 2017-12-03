function [flag,resid1,resid2,rank1,rank2,time1,time2,nv,ps] ...
                   = CompleteZ(m,n,r,problem,opts, ...
                                tolerrank,noisy,noiselevel)
% function [flag,resid1,resid2,rank1,rank2,time1,time2,nv,ps] = ...
%             completeZ(m,n,r,tolerrank,opts)
% CompleteZ is the main function for completing a partial matrix.
% in: 
%   m = # of rows 
%   n = # of cols
%   r = target rank
%   problem = a struct containing Zorig, Zpart, indsZ
%   tolerrank = tolerance for rank function
%   noisy = true / false for noisy or noiseless
%   opts = a struct containing max, min clique sizes and verbose
%   noiselevel = percentage level of noise

% out:
%   flag = 0 -- solved; 1 -- not solved
%   resid1 = residual before refinement 
%   resid2 = residual after refinement

    Zorig     = problem.Zorig;
    Zpart     = problem.Zpart;
    indsZ     = problem.indsZ;        % sampled - positions in Zorig
    mn        = m+n;
    minsize   = opts.minsize;         % clique size
    maxsize   = opts.maxsize;         % clique size
    verbose   = opts.verbose;         % how much output
    Yfinalm   = zeros(m);
    Yfinaln   = zeros(n);

    if ~exist('tolerrank','var')
        if exist('noiselevel','var')
             tolerrank =  max(m,n)*eps(normest(Zpart))+noiselevel; 
        else
             tolerrank =  max(m,n)*eps(normest(Zpart)); 
        end
    end

    starttic = cputime;
    
    %%%%%% forming original adjacency matrix A %%%%%%  
    Z = spalloc(m,n,length(indsZ));
    Z(indsZ) = 1;
    indskeep = find([Z;sparse(n,n)])+(m+n)*m;  % upper-triangular positions in A
    A=[(true(m))   Z
         Z'     (true(n))];
    A(eye(mn)==1) = 0;  % set diag to 0
    
    %%%%%% finding cliques %%%%%%
    ticcliques = tic;
    y = GrowCliques(A,maxsize,minsize);  % list of cliques (type: cell)
    numCliques = length(y);
    if verbose
         fprintf('time for finding %i cliques is %g \n',...
             numCliques,toc(ticcliques))
    end
  
    %%%%%% memory allocation %%%%%%
    expvctrsr = cell(numCliques,1);   % save exposing vectors rows
    expvctrsc = cell(numCliques,1);   % save exposing vectors cols
    expposr = cell(numCliques,1);     % save indices/positions rows
    expposc = cell(numCliques,1);     % save indices/positions cols
    if (noisy)
            expwtsr = zeros(numCliques,1);    % wts rows
            expwtsc = zeros(numCliques,1);    % wts cols
    end
    successr = 0;
    successc = 0;
    maxc = 0;                       % maxsize of the cliques
    
    %%%%%% main loop for exposing vectors %%%%%%
    for ii=1:numCliques        
        indsY = y{ii};                % indics in the clique
        rowsi = logical(indsY<=m);
        rowsi = indsY(rowsi);         % rows of Zorig
        colsj = logical(indsY>m);
        colsj = indsY(colsj)-m;       % cols of Zorig
        lr = length(rowsi);
        lc = length(colsj);

        if  ~((max(indsY)<=m) || (min(indsY)>=m+1)) ...
                              &&  min(lr,lc)>=r &&  max(lr,lc)>r
            if (noisy)
                pqm = min(lr,lc);
            end
            X = full(Zpart(rowsi,colsj));
            [uXtemp,sXtemp,vXtemp] = svd(X);
            diagsX = diag(sXtemp);
            rankX = find(diagsX >= tolerrank*diagsX(1), 1, 'last');
            
            if rankX < r  % should never happen generically
               fprintf('rankX < r - should not happen generically \n');
               keyboard
            end
            
            % finding the largest clique size
            if lr > r || lc > r
                maxc = max(maxc,lr+lc);
            end
            
            % finding row exposing vector
            if lr>r     
                successr=successr+1;
                if (noisy)
                    expwtsr(ii)=sum(diagsX(r+1:pqm).^2)/(.5*lr*(lr-1));
                end
                expposr{ii}=indsY(1:lr);
                expvctrsr{ii}=uXtemp(:,r+1:lr);  % save only vectors
            end
            
            % finding col exposing vector
            if lc>r
                successc=successc+1;
                if (noisy)
                    expwtsc(ii)=sum(diagsX(r+1:pqm).^2)/(.5*lc*(lc-1));
                end
                expposc{ii}=indsY(lr+1:end);
                expvctrsc{ii}=vXtemp(:,r+1:lc);  % save only vectors
            end
        end
    end
  
    if verbose
        fprintf('maxc =  %i; and %i and %i # successful row and col cliques, resp.  \n',...
                        maxc,successr,successc);
    end
    
    %%%%%% calculating weights in NOISY case %%%%%%
    if (noisy)
        sumwtsr = sum(expwtsr);
        if sumwtsr < 1e-13
            exposwtsr = ones(length(expwtsr),1);
            startri = 1;
            indswtsr = 1:length(expwtsr);
        else
            exposwtsr = 1-(expwtsr/sumwtsr);  % w_X^i in paper - exposing wts
            [sortwtsr,indswtsr] = sort(exposwtsr,'ascend'); 
            startri = find(sortwtsr>0,1);
        end

        sumwtsc = sum(expwtsc);
        if sumwtsc < 1e-13
            exposwtsc = ones(length(expwtsc),1);
            startci = 1;
            indswtsc = 1:length(expwtsc);
        else
            exposwtsc = 1-(expwtsc/sumwtsc);  % w_X^i in paper - exposing wts
            [sortwtsc,indswtsc] = sort(exposwtsc,'ascend'); 
            startci = find(sortwtsc>0,1);
        end
    end
    
    %%%%%% forming final exposing vector Yfinal NOISY case %%%%%%
    if (noisy)
        indsjj =indswtsr(startri:end);  % add up exp. vctrs for cols
        if ~isempty(indsjj)
            for jj=reshape(indsjj,1,length(indsjj))
                clique=expposr{jj};
                if ~isempty(clique)
                    temp=expvctrsr{jj}*expvctrsr{jj}';
                    temp=(temp+temp')/2;         
                    Yfinalm(clique,clique)= ... 
                              Yfinalm(clique,clique) + temp*exposwtsr(jj); 
                end
            end
        end
      
        indsjj =indswtsc(startci:end);  % add up exp. vctrs for cols
        if ~isempty(indsjj)
            for jj=reshape(indsjj,1,length(indsjj))
                clique=expposc{jj}-m;
                if ~isempty(clique)
                    temp=expvctrsc{jj}*expvctrsc{jj}';
                    temp=(temp+temp')/2;         
                    Yfinaln(clique,clique)= ... 
                               Yfinaln(clique,clique) + temp*exposwtsc(jj); 
                end
            end
        end    
      
    %%%%%% forming final exposing vector Yfinal NOISELESS case %%%%%%
    else
        for jj=1:length(expvctrsr)  % add up exp. vctrs for cols
            clique=expposr{jj};
            if ~isempty(clique)
                temp=expvctrsr{jj}*expvctrsr{jj}';
                temp=(temp+temp')/2;         
                Yfinalm(clique,clique)= Yfinalm(clique,clique) + temp; 
            end
        end
        for jj=1:length(expvctrsc)  % add up exp. vctrs for cols
            clique=expposc{jj}-m;
            if ~isempty(clique)
                temp=expvctrsc{jj}*expvctrsc{jj}';
                temp=(temp+temp')/2;         
                Yfinaln(clique,clique)= Yfinaln(clique,clique) + temp; 
            end
        end
    end
  
    %%%%%% checking if Yfinal blocks have 0 on the diagonal %%%%%%
    indszeroYr=find(diag(Yfinalm)==0);
    indszeroYc=find(diag(Yfinaln)==0);
    indszeroY= union(indszeroYr,(indszeroYc+m));
    if ~isempty(indszeroYr)
        if verbose
           fprintf('WARNING:  number of zero rows in Yr  %i >0  , ',...
	                 length(indszeroYr)) 
           fprintf(' shifting out diagonal \n')
        end
        indsdiagYr=sub2ind(size(Yfinalm),indszeroYr,indszeroYr);
        Yfinalm(indsdiagYr)=1e1*(rand(length(indsdiagYr),1)+1); % separate eigs
    end
    if ~isempty(indszeroYc)
        if verbose
           fprintf('WARNING:  number of zero rows in Yc  %i >0  , ',...
	                 length(indszeroYc)) 
           fprintf(' shifting out diagonal \n')
        end
        indsdiagYc = sub2ind(size(Yfinaln),indszeroYc,indszeroYc);
        Yfinaln(indsdiagYc) = 1e1*(rand(length(indsdiagYc),1)+1); % separate eigs
    end

  
    if (noisy)
        ticlindep = tic;
        [rYm,~,~,UYm] = lindep(Yfinalm,tolerrank,r);
        rnp = m-rYm;
        if verbose
            fprintf('time %g  nullity %i  for  Yfinalm \n', toc(ticlindep),rnp)
        end
        ticlindep = tic;
        [rYn,~,~,UYn] = lindep(Yfinaln,tolerrank,r);
        rnq = n-rYn;
        if verbose
        fprintf('time %g  nullity %i  for  Yfinaln \n',...
                   toc(ticlindep),rnq)
        end
    else
        ticnull=tic;
        UYm=null(Yfinalm);
        UYn=null(Yfinaln);
        if verbose
          fprintf('time for null of Yfinal %g \n',toc(ticnull))
        end
    end

  
    %%%% finding the final exposing vector %%%%
    rnp = size(UYm,2);
    rnq = size(UYn,2);
    V = [ [UYm ; sparse(n,rnp)]  [ sparse(m,rnq) ; UYn ]  ]; 
    nv = size(V,2);
    if verbose
      fprintf('size of V is %g %g   \n',size(V))
    end
  
    HZ = true(size(Z)); % for non sampled elements if needed
    ps = 100;           % percent sampled to return in function
    if ~isempty(indszeroY) % zero rows in Yfinal - remove rows/cols of Z data
       indszeroZrows = find(indszeroY <= m);
       if ~isempty(indszeroZrows)
           indszeroZrows = indszeroY(indszeroZrows);
	       Z(indszeroZrows,:) = 0;         % ignore rows not sampled enough
	       HZ(indszeroZrows,:) = false;    % ignore rows not sampled enough
       end
       indszeroZcols = find(indszeroY>m);
       if ~isempty(indszeroZcols)
               indszeroZcols = indszeroY(indszeroZcols)-m;
	       Z(:,indszeroZcols) = 0;         % ignore cols not sampled enough
	       HZ(:,indszeroZcols) = false;    % ignore cols not sampled enough
       end
       ps = 100*sum(sum(HZ))/numel(HZ);
       indskeep = find([Z;sparse(n,n)])+(m+n)*m;
     
       % redo adjacency due to new zeros in Z
       A=[(true(m))   Z
            Z'     (true(n))];
       A(eye(mn)==1)=0;  % set diag to 0
     
    end  % end of if for zero rows in Yfinal
  
    [indsi,indsj] = ind2sub(size(A),indskeep);
    indsi2 = indsi;
    indsj2 = indsj-m;
    matrepszz = zeros(length(indskeep),rnp*rnq);    % top right block only
  
    E2 = zeros(nv);
    E2(1:rnp,rnp+1:end) = 1;            % size of top right block of final FR R
    [mi2,mj2] = find(triu(E2,1));       %  subs for upperblock of smaller R
    indsblku = sub2ind(size(E2),mi2,mj2);

    %%% forming the mapping matrix %%%
    UYmT = UYm';   % short version
    UYnT = UYn';   % short version
    sqrt2 = sqrt(2);
    for zz=1:length(indskeep)
	    Ri2=UYmT(:,indsi2(zz));
	    Rj2=UYnT(:,indsj2(zz));
	    Toffdiag2 = (Ri2(mi2).*Rj2(mj2-rnp))/sqrt2;
	    matrepszz(zz,:)=Toffdiag2; % don't need zero columns - just rt blk
    end
  
    tempY=[sparse(m,m) Zpart;sparse(n,m+n)];
    allz=tempY(indskeep); % vector of known entries of Z
  
    if (~noisy)
        [~,indssub,matrepssubzz,~] = lindep(matrepszz',tolerrank);
        matrepssubzz = matrepssubzz';
        allzsub = allz(indssub);
        if rank(matrepssubzz,1e-9)<min(size(matrepssubzz))
            fprintf('WARNING: error rank cvx too small \n')
        end

        warning off
        cvx_clear
      
        if ~verbose
            cvx_quiet true
        end
        
        cvx_begin sdp
        cvx_precision best
        variable R(nv,nv) symmetric
        minimize trace(R)
        subject to
          matrepssubzz*localsvecupblck(R,indsblku) == allzsub;
          R >= 0
        cvx_end
        warning on

        flag = 0;
        YY = V*R*V';
        YY = (YY+YY')/2; % numerically ensure symmetry
        newZ = YY(1:m,m+1:end);
        ranknewZ = rank(full(newZ));    % using default rank here
        ranknewZorig = ranknewZ;
        time1 = cputime - starttic;
        rank1 = ranknewZorig;
        resid1 = norm(HZ.*(newZ-Zorig))/norm(HZ.*Zorig);
        time2 = time1;
        rank2 = rank1;
        resid2 = resid1;
        
    else
      
        SK = randn(2*size(matrepszz,2),size(matrepszz,1)); % random sketch
        SKa = SK*allz;   % RHS using sketch matrix
        SKM = SK*matrepszz;   % matreps using matrix represention by row
      
        gamma = 1e-5; % regularization parameter
        ttcvx=tic;
        warning('off')
        cvx_clear
        
        if(verbose)
          cvx_quiet false
        else
          cvx_quiet true
        end
        
        cvx_begin sdp
        cvx_precision best
        variable R(nv,nv) symmetric   % size rn is rv in paper
        minimize norm(SKM*localsvecupblck(R,indsblku)-SKa)+gamma*norm(R,'fro')
        subject to
          R >= 0
        cvx_end
 
        if(verbose)
            fprintf('time for cvx for simple constr. least squares %g\n', toc(ttcvx));
            fprintf('%i of %i constraints (nonzeros in Zpart) are kept for cvx \n',...
                       size(SK,1),length(indskeep))
            fprintf('for cvx %i vrbles in R size %i; with %i constraints \n',...
                       nv*(nv+1)/2,nv,size(SK,1))
        end
      
        time1 = cputime - starttic;   % time before refinement
      
        %%%%%% Refinement step %%%%%%
        if strcmp(cvx_status,'Solved')
            flag = 0;
            YY=V*R*V';
            newZ=YY(1:m,m+1:end);
            newZH=newZ.*HZ;
            ps=nnz(HZ)/numel(HZ);
            ZorigH=Zorig.*HZ;
            ranknewZ=lindep(newZ);
            rank1 = ranknewZ;
            resid1 = norm(ZorigH-newZH)/norm(ZorigH);
            resid2 = resid1;   % default for after refinement
            rank2 = rank1;     % default for after refinement
    
            [rankRorig,~,~]=lindep(R);
            traceStart = trace(R);
            reftimes = 20;     % maximum number of times to do refinement+1
            refcount = 1;
            refsteps = linspace(1,.05,reftimes); 
      
            tracetars = refsteps*traceStart;    
            residdeltas = resid1;
            ranksZ = ranknewZ;
            ranksR = rankRorig;
            lamtraces = 0; % dual multipliers: the first (for orig) set to 0
            gamma=0.001; %1e-8;  % for regularizatin
      
            warning('off'); 
            while(refcount < reftimes && lamtraces(end) < 0.1)
                if verbose
                    fprintf('***** Doing refinement: %i *****\n', refcount);
                end
                refcount = refcount + 1;
                tracetar = tracetars(refcount); % tracetar starts from 2

                cvx_clear
                if verbose
                    cvx_quiet false
                else
                    cvx_quiet true
                end

                cvx_begin sdp
                cvx_precision best
                variable R(nv,nv) symmetric
                dual variable lamtrace
                minimize norm(SKM*localsvecupblck(R,indsblku)-SKa)+gamma*norm(R,'fro')
                subject to
                    lamtrace : trace(R) <= tracetar;
                    R >= 0
                cvx_end
        
                YY = V*R*V';
                YY = (YY+YY')/2;
                newZ = YY(1:m,m+1:end);
                newZH = newZ.*HZ;

                [ranknewR,~,~] = lindep(R);
                [ranknewZ,~,~] = lindep(newZH);
                newresZ = norm(newZH-ZorigH)/norm(ZorigH);
                residdeltas = [residdeltas newresZ];
                if verbose
                    fprintf('-----rank newZ and new residual value---%i  %d\n',...
                           ranknewZ,newresZ);
                end
                ranksZ = [ranksZ ranknewZ];
                ranksR = [ranksR ranknewR];
                lamtraces = [lamtraces lamtrace];
        
                % plots mainly for debug
                if(verbose)
                    figure(1)
                    semilogy(tracetars(1:refcount),residdeltas,'x-r')
                    set(gca,'xdir','reverse')
                    title('residuals versus trace')

                    figure(2)
                    plot(tracetars(1:refcount),ranksZ,'xr',tracetars(1:refcount), ...
                                         ranksR,'ok')
                    legend('location','best','rank Z','rank R')
                    set(gca,'xdir','reverse')
                    title('ranks')

                    figure(3)
                    semilogy(tracetars(1:refcount),lamtraces,'x-')
                    set(gca,'xdir','reverse')
                    title('dual variable lamtrace')
                    drawnow
                end
              
                rank2 = min(rank2,ranknewZ);
                resid2 = min(resid2,newresZ);
            end    % end of while
            time2 = cputime - starttic;
        else
            flag = 1;
            resid1 = 0;
            resid2 = 0;
            rank1 = 0;
            rank2 = 0;
            time1 = 0;
            time2 = 0;
        end
    end
  
end % end of function completeZ


function [r,idx,Xsub,Xnull] = lindep(X,tol,trgtr)
% lindep extracts a linearly independent set of columns of a given matrix X
% in:
%   X = input matrix
%   tol = A rank estimation tolerance, default=1e-10
% out:
%   Xsub = extracted columns of X
%   idx =  indices (into X) of the extracted columns 

    if ~nnz(X) % X has no non-zeros and hence no independent columns
        Xsub=[]; 
        idx=[];
        return
    end
    
    if nargin<2, tol=1e-10; end
    
    [Q, R, E] = qr(X,0);
    if ~isvector(R)
        diagr = abs(diag(R));
    else
        diagr = R(1);
    end
    
    r = find(diagr >= tol*diagr(1), 1, 'last'); % rank estimation
    
    if nargin >2
	       r = max(r,size(X,2)-3*trgtr);   % ensure nullspace not too large
    end
    
    if nargout > 1
        idx = sort(E(1:r));
        idx = idx(:);
    end
    
    if nargout > 2
        Xsub = X(:,idx);                      
    end
    
    idx = sort(E(1:r)); 
    Xnull = Q(:,r+1:end);
    
  end % end of function lindep
