function [Tout] = combineout(basename)

%basename = '/Users/tjc/research/crackingfronts/in_out/k215e15/k215e15_stead03*';

files = dir(basename);
i = 1;
for file = files'
  slashloc = findstr('/',basename);
  dirname = basename(1:slashloc(end));
  filename = [dirname, file.name];
  disp(filename)
  load(filename);
  if i==1
    Tout = T2;
  else
    Tout(:,:,i) = T2;
  end
  i = i+1;
end
