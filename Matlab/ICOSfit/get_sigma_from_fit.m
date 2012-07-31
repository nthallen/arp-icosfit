function [sigma,aveP,cpci]=get_sigma_from_fit(base)

%function sigma=get_sigma_from_fit('ICOSout')

data=load([base '/ICOSsum.dat']);
cpci=data(:,6);
sigma=ones(1,length(cpci))*NaN;
aveP=ones(1,length(cpci))*NaN;

for i =1:length(cpci)
    path = mlf_path(base, cpci(i));
    f=load(path);
    sigma(i)=std((f(:,3)-f(:,4))./f(:,5));
    aveP(i)=mean(f(:,5));
end
