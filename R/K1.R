# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"K1" <-
function(S=35,T=25,P=0,k1k2='r',phflag=0){

if (k1k2=='m')
{
k1k2flag=1;
}
if (k1k2=='r')
{
k1k2flag=0;
}
#-------Constantes----------------

#---- issues de equic----
tk = 273.15;           # [K] (for conversion [deg C] <-> [K])
TC = T + tk;           # TC [C]; T[K]
Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
cl3 = Cl^(1/3);   
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl * Cl;   # ionic strength
iom0 = 19.924*S/(1000-1.005*S);
ST = 0.14/96.062/1.80655*S;   # (mol/kg soln) total sulfate


bor = (416.*(S/35.))* 1e-6;   # (mol/kg), DOE94   

	
# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [CO2] = K_1
#
#   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 14)
#   pH-scale: 'total'. mol/kg-soln

tmp1 = 2.83655 - 2307.1266 / TC - 1.5529413 * log(TC);
tmp2 =         - (0.20760841 + 4.0484 / TC) * sqrt(S);
tmp3 =         + 0.08468345 * S - 0.00654208 * S * sqrt(S);   
tmp4 =         + log(1 - 0.001005 * S);

lnK1roy = tmp1 + tmp2 + tmp3 + tmp4;

if (phflag == 0)
{
        K1roy  = exp(lnK1roy);
}
if (phflag == 1)
{
        lnK1roy = lnK1roy-log(total2free);
        K1roy   = exp(lnK1roy);
}
	
		
# --------------------- K1 ---------------------------------------
#   first acidity constant:
#   [H^+] [HCO_3^-] / [H_2CO_3] = K_1
#
#   Mehrbach et al (1973) refit by Lueker et al. (2000).
#
#   pH-scale: 'total'. mol/kg-soln

pK1mehr = 3633.86/TC - 61.2172 + 9.6777*log(TC) - 0.011555*S + 0.0001152*S*S;
	
if (phflag == 0)
{
        K1mehr  = 10^(-pK1mehr);
}
if (phflag == 1)
{
	lnK1mehr = log(10^(-pK1mehr))-log(total2free);
        K1mehr   = exp(lnK1mehr);
}
	
#----------- Roy or Mehrbach. default: Roy 

K1 = K1roy;
	
if (exists('k1k2flag'))
{
	if (k1k2flag == 0)
	{	
	K1 = K1roy;
	}
	if (k1k2flag == 1)
	{
	K1 = K1mehr;
	}
	
}
	
if (P > 0.0)
{
		
	RGAS = 8.314510;        # J mol-1 deg-1 (perfect Gas)  
	R = 83.131;             # mol bar deg-1 
	                        # conversion cm3 -> m3          *1.e-6
        	                #            bar -> Pa = N m-2  *1.e+5
	                        #                => *1.e-1 or *1/10
		
		
	# index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8,
	#        K1P 9, K2P 10, K3P 11
	
	#----- note: there is an error in Table 9 of Millero, 1995.
	#----- The coefficients -b0 and b1
	#----- have to be multiplied by 1.e-3!

	#----- there are some more errors! 
	#----- the signs (+,-) of coefficients in Millero 95 do not
	#----- agree with Millero 79
	
		
		
	a0 = c(-25.5, -15.82, -29.48, -25.60, -18.03, -9.78, -48.76, -46, -14.51, -23.12, -26.57);
	a1 = c(0.1271, -0.0219, 0.1622, 0.2324, 0.0466, -0.0090, 0.5304, 0.5304, 0.1211, 0.1758, 0.2020);
	a2 = c(0.0, 0.0, 2.608*1e-3, -3.6246*1e-3, 0.316*1e-3, -0.942*1e-3, 0.0, 0.0, -0.321*1e-3, -2.647*1e-3, -3.042*1e-3);
	b0 = c(-3.08*1e-3, 1.13*1e-3, -2.84*1e-3, -5.13*1e-3, -4.53*1e-3, -3.91*1e-3, -11.76*1e-3, -11.76*1e-3, -2.67*1e-3, -5.15*1e-3, -4.08*1e-3);
	b1 = c(0.0877*1e-3, -0.1475*1e-3, 0.0, 0.0794*1e-3, 0.09*1e-3, 0.054*1e-3, 0.3692*1e-3, 0.3692*1e-3, 0.0427*1e-3, 0.09*1e-3, 0.0714*1e-3);
	b2 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	
	deltav = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	deltak = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	lnkpok0 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	
	for (ipc in 1:length(a0))
	{
	  deltav[ipc]  =  a0[ipc] + a1[ipc] *T + a2[ipc] *T*T;
	  deltak[ipc]   = (b0[ipc]  + b1[ipc] *T + b2[ipc] *T*T);  
	  lnkpok0[ipc]  = -(deltav[ipc] /(R*TC))*P + (0.5*deltak[ipc] /(R*TC))*P*P;
	}
	
	K1 = K1*exp(lnkpok0[1]);
		
	}
#attr(K1,"unit") = "mol/kg"
#print(K1);
}
