function [Unew,V,S]=basefile(region,nsvds,suffix,range)
%function [U,V,S]=basefile(region,nsvds,suffix)
[data]=loadscans([],region);
if nargin<4
    range=300:1831;
end
[U,S,V]=svds(data(range,:),nsvds);
figure; semilogy(diag(S),'*');
figure; for i=1:size(U,2); nsubplot(size(U,2),1,i,1); plot(U(:,i)); end
addzoom
Unew = zeros(length(data),size(U,2));
Unew(range,:)=U;
writebase('sbase.dat',Unew,S,V);
if nargin>=3
    eval(sprintf('!mv sbase.dat sbase%i.%s.dat',region(1),suffix))
    disp(sprintf('Moved sbase file to sbase%i.%s.dat\n',region(1),suffix))
end
