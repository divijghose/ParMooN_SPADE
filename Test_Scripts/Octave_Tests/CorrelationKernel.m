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
  r = sqrt(((xx'(:)-xx'(:)').^2)+((yy'(:)-yy'(:)').^2));
  rByL = r./LenScale;
  Cr = (1+rByL+(rByL.^2)).*(exp(-1.0*rByL));
  sigx = (exp((-1.0/E).*((2.*xx-1-Disp).^Power)))./(2*pi*sqrt(E));
  sigy = (exp((-1.0/E).*((2.*yy-1-Disp).^Power)))./(2*pi*sqrt(E));
  sig = sigx.*sigy;
  sigtot = sig'(:)*sig'(:)';
  retval = Cr.*sigtot;

