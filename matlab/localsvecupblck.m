function x = localsvecupblck(X,indsblku)
% Use the isometry version  for the str. upper rt blck
% in:
%   X = a symmetric matrix
%   indsblku = indices of upper block
x = sqrt(2) * X(indsblku);
end
