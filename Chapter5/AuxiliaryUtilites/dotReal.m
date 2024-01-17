function res = dotReal( A,B )
% dotReal(A,B) returns the real dot product of the vectors A and B
%
% Input:
% A - vector
% B - vector
%
% Output:
% res - the result of the real dot product
    res=sum(A.*B,1);

end

