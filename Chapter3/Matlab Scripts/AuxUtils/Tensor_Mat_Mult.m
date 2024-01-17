function C = Tensor_Mat_Mult(A,B)
m = size(A,1);
n = size(B,2);
k = size(A,3);
C = zeros(m,n,k);
for j = 1:size(A,2)
  C = C + bsxfun(@times,A(:,j,:),B(j,:,:));
end

end