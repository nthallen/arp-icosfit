function rewrite_ssp
ibase = 'SSP';
obase = 'CPCI14';
index = 1;
while 1
  pi = mlf_path(ibase,index);
  fi = loadbin(pi);
  if isempty(fi) break; end
  po = mlf_path(obase,index);
  fo = fi * [ 1 0 0; 0 1 0; 0 0 -1];
  temp=fo(:,1);
  fo(:,1)=fo(:,3);
  fo(:,3)=temp;
  mlf_mkdir(obase,index);
  writebin( po, fo );
  index = index+1;
end
