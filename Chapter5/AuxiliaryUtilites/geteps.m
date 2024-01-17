function epsr = geteps(lambda)
% Input: lambda in micrometer
% Output: epsr the relative electric permittivity

%This is from the work of Johnson and Christy
fileID = fopen('Christy.txt','r');   
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
wavelengths = A(1:3:end);
n = A(2:3:end);
k = A(3:3:end);

Reeps = interp1(wavelengths,n,lambda);
Imeps = interp1(wavelengths,k,lambda);

epsr = (Reeps + 1i*Imeps)^2;
fclose(fileID);
end