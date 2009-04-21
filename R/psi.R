# Copyright (C) 2009 Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"psi" <-
function(flag, var1, var2, S=35, T=20, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l"){
		buf <- buffer(flag=flag, var1=var1, var2=var2, S=S, T=T, P=P, Pt=Pt, Sit=Sit, pHscale=pHscale, kf=kf, k1k2=k1k2)
		psi <- -buf$PiC/buf$PiD
out <- psi
attr(out,"unit") <- "mol CO2/mol CaCO3"
return(out)
}