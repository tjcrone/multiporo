function [output] = combineout(basename)
% This function combines the new output files into a merged
% structure for plotting and animations.
%
% Example:
%
% >> output = combineout('~/research/crackingfronts/in_out/k215e16/k215e16_stead03_out_*');

files = dir(basename);
i = 1;
for file = files'
  slashloc = findstr('/',basename);
  dirname = basename(1:slashloc(end));
  filename = [dirname, file.name];
  load(filename);
  if i==1
    output.Tout = T2;
    output.rhofout = rhof2;
    output.cfout = cf2;
    output.qzout = qz2;
    %crackedout = cracked;
    output.tout = t;
  else
    output.Tout(:,:,i) = T2;
    output.rhofout(:,:,i) = rhof2;
    output.cfout(:,:,i) = cf2;
    output.qzout(:,:,i) = qz2;
    %output.crackedout(:,:,i) = cracked;
    output.tout(i) = t;
  end
  i = i+1;
end
