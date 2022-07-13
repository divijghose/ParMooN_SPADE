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
## @deftypefn {} {@var{retval} =} DO_Init (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Divij <divij@divij>
## Created: 2022-07-07

function [DO_mean, coeff,DO_mode] = DO_Init (RealznVect, N_R)
  DO_mean = mean(RealznVect')';
  pertVector = RealznVect - DO_mean;
  [Ud,Sd,Vtd] = svd(pertVector);
  dSd = diag(Sd);
  subDim=1;
  while(cumsum(dSd)(subDim)/sum(dSd)<0.999)
    subDim++;
  endwhile
  projVector = pertVector'*Ud;
  coeff = projVector(:,1:subDim);
  DO_mode = Ud(:,1:subDim);
endfunction
