## Copyright (C) 2022 Divij
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} CovarianceKernel (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Divij <divij@divij>
## Created: 2022-07-05

function retval = CovarianceKernel (xx,yy,Refn,LenScale,E,Disp,Power)
  N = (2^Refn)+1;
  retval = ones(N^2,N^2);
  xRowMaj = xx'(:)';
  yRowMaj = yy'(:)';
  for i = 1:N^2
    for j = 1:N^2
      r = sqrt(((xRowMaj(i)-xRowMaj(j))^2)+((yRowMaj(i)-yRowMaj(j))^2));
      rByL = r/LenScale;
      retval(i,j) = (1.0+rByL+((rByL^2)/3.0))*exp(-rByL);
      sig_r1 = (exp(-1.0*power(2.0*xRowMaj(i)-1.0-Disp,Power)/E)/(2.0*pi*sqrt(E)))*(exp(-1.0*power(2.0*yRowMaj(i)-1.0-Disp,Power)/E)/(2.0*pi*sqrt(E)));
      sig_r2 = (exp(-1.0*power(2.0*xRowMaj(j)-1.0-Disp,Power)/E)/(2.0*pi*sqrt(E)))*(exp(-1.0*power(2.0*yRowMaj(j)-1.0-Disp,Power)/E)/(2.0*pi*sqrt(E)));
      retval(i,j)*=sig_r1*sig_r2;
    endfor
  endfor
endfunction
