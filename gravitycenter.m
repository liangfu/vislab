function pt = gravitycenter(x)
% GRAVITYCENTER get spatial center of given binary image
  [sy sx]=find(x);
  pt=[mean(sy),mean(sx)];

