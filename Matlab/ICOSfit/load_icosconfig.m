function S = load_icosconfig( base_ );
oldpwd_ = pwd;
cd(base_);
ICOSconfig;
cd(oldpwd_);
vars_ = who;
for i=1:length(vars_)
  var = vars_{i};
  if var(end) ~= '_'
    S.(var) = eval(var);
  end
end
