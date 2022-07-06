function [xx,yy] = getMesh (N);

  x = linspace(0,1,(2^N)+1);
  y = linspace(0,1,(2^N)+1);
  [xx,yy] = meshgrid(x,y);
endfunction


