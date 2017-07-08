function [Tout, rhofout, cfout, qzout, crackedout, tout] = combineout(basename)

%basename = '/Users/tjc/research/crackingfronts/in_out/k215e15/k215e15_stead03*';

files = dir(basename);
i = 1;
for file = files'
  slashloc = findstr('/',basename);
  dirname = basename(1:slashloc(end));
  filename = [dirname, file.name];
  load(filename);
  if i==1
    Tout = T2;
    rhofout = rhof2;
    cfout = cf2;
    qzout = qz2;
    crackedout = cracked;
    tout = t;
  else
    Tout(:,:,i) = T2;
    rhofout(:,:,i) = rhof2;
    cfout(:,:,i) = cf2;
    qzout(:,:,i) = qz2;
    crackedout(:,:,i) = cracked;
    tout(:,:,i) = t;
  end
  i = i+1;
end
