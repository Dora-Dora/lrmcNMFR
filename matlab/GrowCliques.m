function listCliques = GrowCliques(A,maxsize,minsize)
% in: 
%   A = an adjacency matrix
% 	maxsize = max size of clique to return
%   minsize = min size of clique to return
% out:
%   A list of cliques (without replications) (type: cell) such that
%   the size is between [minsize,maxsize]

    n = size(A,1);
    II = zeros(n*maxsize,1);
    JJ = zeros(n*maxsize,1);
    k = 0; % count # of edges
    
    for jc = 1:n

        k = k+1; II(k) = jc; JJ(k) = jc;
        NodesToAdd = maxsize-1;
        q = (A(:,jc) > 0);
        
        while NodesToAdd && nnz(q)
            NodesToAdd = NodesToAdd - 1;
            ic = find(q); 
            k = k + 1;  JJ(k) = jc;
            
            % alternating the index that to be added,
            % so that the vertices will gather at
            % all <=m or all > m
            if mod(NodesToAdd,2) == 0
                II(k) = ic(1);
                q = q & A(:,ic(1));
            else
                II(k) = ic(end);
                q = q & A(:,ic(end));
            end
        end
    end
    
    Cliques = sparse(II(1:k),JJ(1:k),1,n,n);   
    Cliques = (unique(Cliques','rows'))';
    
    inds = find(sum(Cliques,1) >= minsize);    
    numCliques = length(inds);    
    listCliques = cell(numCliques,1);
    
    for ii = 1:numCliques
        listCliques{ii} = find(Cliques(:,inds(ii)));
    end

end
