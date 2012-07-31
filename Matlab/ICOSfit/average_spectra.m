function average_spectra(n_avg)
runn=pwd;
runn=runn(max(strfind(runn,'/'))+1:end);
base=['../' runn '.average' num2str(n_avg)];
mkdir('../',[runn '.average' num2str(n_avg)]);
mkdir(base,'cd')
mkdir([base '/cd'],'home')
mkdir([base '/cd/home'],'CR')
mkdir([base '/cd/home/CR'],[runn '.average' num2str(n_avg)])
mkdir([base '/cd/home/CR'],'anal')
mkdir([base '/cd/home/CR/anal'],[runn '.average' num2str(n_avg)])

eval(['!cp cd/home/CR/anal/' getrun(1) '/* ' base '/cd/home/CR/anal/' runn '.average' num2str(n_avg)]);

base=[base '/cd/home/CR/' runn '.average' num2str(n_avg)];

mkdir(base,'Base')
eval(['!cp cd/home/CR/' getrun(1) '/Base/* ' base '/Base']);

mkdir(base,'CPCI');
mkdir([base '/CPCI'],'CPCI14');
baseCPCI=[base '/CPCI/CPCI14/'];
PTE=load('PTE.txt');
l=PTE(end,1);
h=0; m=0;
while (h*3600 + m*60) < l
    mkdir(baseCPCI,sprintf('%02i',h))
    while m < 60 & (h*3600 + m*60) < l
        mkdir(strcat(baseCPCI,sprintf('%02i',h)),sprintf('%02i',m))
        m=m+1;
    end
    h=h+1; m=0;
end
PTEnew=[];
breaks=[1;find(diff(PTE(:,1))>1);size(PTE,1)];
for i = 1:length(breaks)-1
    j=PTE(breaks(i)+1,1);
    while j < (PTE(breaks(i+1)-1,1)-n_avg)
        try [icos,etln]=base2([],j:j+n_avg-1); 
        icos=mean(icos');
        etln=mean(etln');
        writebin(mlf_path(baseCPCI,j),[icos;etln]');
        PTEnew(end+1,:)=PTE(find(PTE(:,1)==j),:);
        catch
        end
        j=j+n_avg;
    end
end
cd(['../' runn '.average' num2str(n_avg)]);
save -ascii 'PTEnew.txt' 'PTEnew'
eval(['!ln -s cd/home/CR/anal/' getrun(1) '/*.mat .'])
eval(['!cp ../sbase* .'])
if exist(['../' runn '/MirrorLoss.mat'])
    eval(['!cp ../' runn '/MirrorLoss.mat .'])
end
if exist(['../' runn '/fitline.dat'])
    eval(['!cp ../' runn '/fitline.dat .'])
end
if exist(['../' runn '/waves.m'])
    eval(['!cp ../' runn '/waves.m .'])
end
eval(['!cp ../' runn '/*_etln.mat .'])


