function hdrs = loadscanhdrs(scans)
% hdrs = loadscanhdrs(scans)
base = find_scans_dir([]);
p = mlf_path(base, scans(1));
[~,hdr] = loadbin(p);
flds = fields(hdr);
z = cell(size(scans));
args = cell(length(flds)*2,1);
args(1:2:end) = flds;
args(2:2:end) = {z};
hdrs(1) = hdr;
for i=2:length(scans)
    p = mlf_path(base,scans(i));
    [~,hdr] = loadbin(p);
    hdrs(i) = hdr;
end
