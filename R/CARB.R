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
#
"carb" <-
function(flag=8,var1=8.2,var2=2400,S=35,T=25,P=0,k1k2='r',phflag=0,ini='s'){


if (k1k2=='m')
{
k1k2flag=1;
}
if (k1k2=='r')
{
k1k2flag=0;
}

resultat="carb.out";

# flag = 1      pH-CO2 given
# flag = 2      CO2-HCO3 given
# flag = 3      CO2-CO3 given  
# flag = 4      CO2-ALK given
# flag = 5      CO2-DIC given
# flag = 6      pH and HCO3 given
# flag = 7      pH and CO3 given
# flag = 8      pH and ALK given
# flag = 9      pH and DIC given
# flag = 10     HCO3 and CO3 given
# flag = 11     HCO3 and ALK given
# flag = 12     HCO3 and DIC given
# flag = 13     CO3 and ALK given
# flag = 14     CO3 and DIC given
# flag = 15     ALK and DIC given

if (ini=='s')
{
	if (flag==1)
	{
	PH=var1;
	CO2=var2/1000000;
	}
	
	if (flag==2)
	{
	CO2=var1/1000000;
	HCO3=var2/1000000;
	}
	
	if (flag==3)
	{
	CO2=var1/1000000;
	CO3=var2/1000000;
	}
	
	if (flag==4)
	{
	CO2=var1/1000000;
	ALK=var2/1000000;
	}
	
	if (flag==5)
	{
	CO2=var1/1000000;
	DIC=var2/1000000;
	}
	
	if (flag==6)
	{
	PH=var1;
	HCO3=var2/1000000;
	}
	
	if (flag==7)
	{
	PH=var1;
	CO3=var2/1000000;
	}
	
	if (flag==8)
	{
	PH=var1;
	ALK=var2/1000000;
	}
	
	if (flag==9)
	{
	PH=var1;
	DIC=var2/1000000;
	}
	
	if (flag==10)
	{
	HCO3=var1/1000000;
	CO3=var2/1000000;
	}
	
	if (flag==11)
	{
	HCO3=var1/1000000;
	ALK=var2/1000000;
	}
	
	if (flag==12)
	{
	HCO3=var1/1000000;
	DIC=var2/1000000;
	}
	
	if (flag==13)
	{
	CO3=var1/1000000;
	ALK=var2/1000000;
	}
	
	if (flag==14)
	{
	CO3=var1/1000000;
	DIC=var2/1000000;
	}
	
	if (flag==15)
	{
	ALK=var1/1000000;
	DIC=var2/1000000;
	}
n=1;
}


#phflag =  ;

#phflag=winDialogString("Choise of pH scale          0:Total scale          1:Free scale",phflag );
#k1k2flag <<- "";
#k1k2flag=winDialogString("choose K1 and K2:          0:Roy          1:Mehrbach",k1k2flag);


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


if (ini=='f')
{
fich<-file.choose();
Data <- read.delim(fich, sep="\t", header=TRUE);
l=nchar(fich);
st=l-3;
r=substr(fich,1,st);
resu=paste(r,"out",sep="");
TPH=Data$PH;
TCO2=Data$CO2;
THCO3=Data$HCO3;
TCO3=Data$CO3;
TDIC=Data$DIC;
TALK=Data$ALK;
TS=Data$S;
TT=Data$T;
TP=Data$P;
n=nrow(Data);
}




for (i in 1:n)
{
	if (ini=='f')

	{
	PH=TPH[i];
	CO2=TCO2[i];
	HCO3=THCO3[i];
	CO3=TCO3[i];
	DIC=TDIC[i];
	ALK=TALK[i];
	S=TS[i];
	T=TT[i];
	P=TP[i];
	}
	


	#---------------------------------------------------------------------
	#------------- calcul des K ----------------------------------
	#---------------------------------------------------------------------
		
	#--------------------------------------------------------------
	#------------------ Ks ----------------------------------------
	#       Dickson and Goyet (1994), Chapter 5, p.13
	#       (required for total2free)
	#       Equilibrium constant for HSO4- = H+ + SO4--
	#
	#       K_S  = [H+]free [SO4--] / [HSO4-]
	#       pH-scale: free scale !!!
	#
	#       the term log(1-0.001005*S) converts from
	#       mol/kg H2O to mol/kg soln
	
	
	tmp1 = -4276.1 / TC + 141.328 -23.093*log(TC);
	tmp2 = +(-13856 / TC + 324.57 - 47.986 * log(TC))*sqrt(iom0);
	tmp3 = +(35474 / TC - 771.54 + 114.723 * log(TC))*iom0;
	tmp4 = -2698 / TC *sqrt(iom0)*iom0 + 1776 / TC *iom0 *iom0;
	                                       
	
	lnKs = tmp1 + tmp2 + tmp3 + tmp4 + log(1-0.001005*S);
	
	Ks = exp(lnKs);
	#------- total2free -----------------------------------------------
	#
	#       convert from pH_total ('total`) to pH ('free`):
	#       pH_total = pH_free - log(1+ST/KS(s,tk))
	
	total2free = 1+ST/Ks;
	
	
	#---------------------------------------------------------------------
	# --------------------- Kf  ------------------------------------------
	#  Kf = [H+][F-]/[HF]  
	#
	#   (Dickson and Riley, 1979 in Dickson and Goyet, 
	#   1994, Chapter 5, p. 14)
	#   pH-scale: 'total'   
	
	
	tmp1 = 1590.2/TC - 12.641 + 1.525*sqrt(ION);
	tmp2 = log(1-0.001005*S) + log(1+ST/Ks);
	
	
	lnKf = tmp1 + tmp2;
	if (phflag == 0)
	{
	        Kf  = exp(lnKf);
	}
	if (phflag == 1)
	{
	        lnKf = lnKf-log(total2free);
	        Kf  = exp(lnKf);
	}



	#------- sws2free -----------------------------------------------
	#
	#       convert from pH_sws ('seawater scale`) to pH ('free`):
	#       pH_sws = pH_free - log(1+S_T/K_S(S,T)+F_T/K_F(S,T))
	
	
	FT = 7e-5*(S/35);
	sws2free  = (1+ST/Ks+FT/Kf);
	corr = sws2free/total2free;
	
	
	#-------------------------------------------------------------------
	# --------------------- Kwater -------------------------------------
	#
	#       Millero (1995)(in Dickson and Goyet (1994, Chapter 5, p.18))
	#       $K_w$ in mol/kg-soln.
	#       pH-scale: pH$_{total}$ ('total` scale).
	                                                     
	
	tmp1 = -13847.26/TC + 148.96502 - 23.6521 * log(TC);
	tmp2 = + (118.67/TC - 5.977 + 1.0495*log(TC))*sqrt(S) - 0.01615*S;
	
	lnKw =  tmp1 + tmp2;
	
	if (phflag == 0)
	{
	        Kw  = exp(lnKw);
	}
	if (phflag == 1)
	{
	        lnKw = lnKw-log(total2free);
	        Kw  = exp(lnKw);
	}
	
	
	#---------------------------------------------------------------------
	#---------------------- Kh (K Henry) ---------------------------------
	#
	#               CO2(g) <-> CO2(aq.)
	#               Kh      = [CO2]/ p CO2
	#
	#   Weiss (1974)   [mol/kg/atm]
	#
	#
	tmp = 9345.17 / TC - 60.2409 + 23.3585 * log(TC/100);
	nKhwe74 = tmp + S*(0.023517-0.00023656*TC+0.0047036e-4*TC*TC);
	

	Kh= exp(nKhwe74);
	
	
		
	#---------------------------------------------------------------------------------
	#---------------------- Choix des constantes K1 et K2 ----------------------------
	#---------------------------------------------------------------------------------
	
		
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
	
	
	# --------------------- K2 ----------------------------------------
	#
	#   second acidity constant:
	#   [H^+] [CO_3^--] / [HCO_3^-] = K_2
	#
	#   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 15)
	#   pH-scale: 'total'. mol/kg-soln
	
	tmp1 = -9.226508 - 3351.6106 / TC - 0.2005743 * log(TC);
	tmp2 = (-0.106901773 - 23.9722 / TC) * sqrt(S);
	tmp3 = 0.1130822 * S - 0.00846934 * S^1.5 + log(1 - 0.001005 * S);
	
	lnK2roy = tmp1 + tmp2 + tmp3;
	
	if (phflag == 0)
	{
	        K2roy  = exp(lnK2roy);
	}
	
	if (phflag == 1)
	{
	        lnK2roy = lnK2roy-log(total2free);
	        K2roy   = exp(lnK2roy);
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
	
	
	# --------------------- K2 ----------------------------------------
	#
	#   second acidity constant:
	#   [H^+] [CO_3^--] / [HCO_3^-] = K_2
	#
	#   Mehrbach et al. (1973) refit by Lueker et al. (2000).
	#
	#   pH-scale: 'total'. mol/kg-soln
	
	pK2mehr = 471.78/TC + 25.9290 - 3.16967*log(TC) - 0.01781*S + 0.0001122*S*S;
	
	if (phflag == 0)
	{
	        K2mehr  = 10^(-pK2mehr);
	}
	if (phflag == 1)
	{
		lnK2mehr = log(10^(-pK2mehr))-log(total2free);
	        K2mehr   = exp(lnK2mehr);
	}
	
	
	#----------- Roy or Mehrbach. default: Roy 

	K1 = K1roy;
	K2 = K2roy;
	
	if (exists('k1k2flag'))
	{
		if (k1k2flag == 0)
		{	
		K1 = K1roy;
		K2 = K2roy;
		}
	
		if (k1k2flag == 1)
		{
		K1 = K1mehr;
		K2 = K2mehr;
		}
	
	}



	#---------------------------------------------------------------------
	# --------------------- Kb  --------------------------------------------
	#  Kbor = [H+][B(OH)4-]/[B(OH)3]
	#
	#   (Dickson, 1990 in Dickson and Goyet, 1994, Chapter 5, p. 14)
	#   pH-scale: 'total'. mol/kg-soln
	
	
	tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S^(3/2)-0.0996*S*S);
	tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
	tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S)*log(TC);
	
	lnKb = tmp1 / TC + tmp2 + tmp3 + 0.053105*sqrt(S)*TC;
	
	if (phflag == 0)
	{
	        Kb  = exp(lnKb);
	}
	if (phflag == 1)
	{
	        lnKb = lnKb-log(total2free);
	        Kb  = exp(lnKb);
	}
	



	# --------------------- Phosphoric acid ---------------------
	#
	#   (DOE, 1994)  (Dickson and Goyet): pH_T, mol/(kg-soln)
	#   Ch.5 p. 16
	#
	
	lnK1P = -4576.752 / TC + 115.525 - 18.453*log(TC) + (-106.736 / TC + 0.69171) * sqrt(S) + (-0.65643 / TC - 0.01844) * S;
	lnK2P = -8814.715 / TC + 172.0883 - 27.927 * log(TC) + (-160.34 / TC + 1.3566) * sqrt(S) + (0.37335 / TC - 0.05778) * S;
	lnK3P = -3070.75 / TC - 18.141 + (17.27039 / TC + 2.81197) * sqrt(S) + (-44.99486 / TC - 0.09984) * S;
	
	K1P = exp(lnK1P);
	K2P = exp(lnK2P);
	K3P = exp(lnK3P);
	
	# --------------------- Silicic acid ---------------------------
	#
	#   (DOE, 1994)  (Dickson and Goyet): pH_T, mol/(kg-soln)
	#   Ch.5 p. 17
	#
	
	
	lnKSi = -8904.2 / TC + 117.385 - 19.334*log(TC) + (3.5913-458.79 / TC) * sqrt(iom0) + (188.74 / TC - 1.5998) * iom0 + (0.07871 - 12.1652 / TC) *iom0^2 + log(1-0.001005*S);
	
	KSi = exp(lnKSi);
	
	

	# --------------------- Kspc (calcite) ----------------------------
	#
	# apparent solubility product of calcite
	#
	#  Kspc = [Ca2+]T [CO32-]T
	#
	#  where $[]_T$ refers to the equilibrium total 
	# (free + complexed) ion concentration.
	#
	#  Mucci 1983 mol/kg-soln
	
	tmp1 = -171.9065-0.077993*TC+2839.319/TC+71.595*log10(TC);
	tmp2 = +(-0.77712+0.0028426*TC+178.34/TC)*sqrt(S);
	tmp3 = -0.07711*S+0.0041249*S^1.5;
	log10Kspc = tmp1 + tmp2 + tmp3;
	
	Kspc = 10^(log10Kspc);

	# --------------------- Kspa (aragonite) ----------------------------
	#
	# apparent solubility product of aragonite
	#
	#  Kspa = [Ca2+]T [CO32-]T
	#
	#  where $[]_T$ refers to the equilibrium total 
	# (free + complexed) ion concentration.
	#
	#  Mucci 1983 mol/kg-soln
	
	tmp1 = -171.945-0.077993*TC+2903.293/TC+71.595*log10(TC);
	tmp2 = +(-0.068393+0.0017276*TC+88.135/TC)*sqrt(S);
	tmp3 = -0.10018*S+0.0059415*S^1.5;
	log10Kspa = tmp1 + tmp2 + tmp3;
	
	Kspa = 10^(log10Kspa);


#----------------------------------------------------
# Density of seawater as function of S,T,P.
#
# Millero et al. 1981, Gill, 1982.
#----------------------------------------------------


	#------------ Density of pure water
	
	rhow = 999.842594 + 6.793952e-2*T -9.095290e-3*T^2 + 1.001685e-4*T^3 -1.120083e-6*T^4 + 6.536332e-9*T^5;
	
	#------------ Density of seawater at 1 atm, P=0
	
	A = 8.24493e-1 - 4.0899e-3*T + 7.6438e-5*T^2 - 8.2467e-7*T^3 + 5.3875e-9*T^4; 
	B = -5.72466e-3 + 1.0227e-4*T - 1.6546e-6*T^2; 	
	C = 4.8314e-4;   
	
	rho0 = rhow + A*S + B*S^(3/2) + C*S^2;
	
	
	#-------------- Secant bulk modulus of pure water 
	#
	# The secant bulk modulus is the average change in pressure 
	# divided by the total change in volume per unit of initial volume.
	
	
	Ksbmw = 19652.21 + 148.4206*T - 2.327105*T^2 + 1.360477e-2*T^3 - 5.155288e-5*T^4;
	
	#-------------- Secant bulk modulus of seawater at 1 atm
	
	Ksbm0 = Ksbmw + S*( 54.6746 - 0.603459*T + 1.09987e-2*T^2 - 6.1670e-5*T^3) + S^(3/2)*( 7.944e-2 + 1.6483e-2*T - 5.3009e-4*T^2);
	
	
	#-------------- Secant bulk modulus of seawater at S,T,P
		
	Ksbm = Ksbm0 + P*( 3.239908 + 1.43713e-3*T + 1.16092e-4*T^2 - 5.77905e-7*T^3) + P*S*( 2.2838e-3 - 1.0981e-5*T - 1.6078e-6*T^2) + P*S^(3/2)*1.91075e-4 + P*P*(8.50935e-5 - 6.12293e-6*T + 5.2787e-8*T^2) + P^2*S*(-9.9348e-7 + 2.0816e-8*T + 9.1697e-10*T^2);
	 	
	
	#------------- Density of seawater at S,T,P
	
	rho = rho0/(1-P/Ksbm);
	
	
	#---------------------- Pressure effect on K's (Millero, 95) ----------#
	
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
	K2 = K2*exp(lnkpok0[2]);
	Kb = Kb*exp(lnkpok0[3]);
	Kw = Kw*exp(lnkpok0[4]);
	Ks = Ks*exp(lnkpok0[5]);
	Kf = Kf*exp(lnkpok0[6]);
	Kspc = Kspc*exp(lnkpok0[7]);
	Kspa = Kspa*exp(lnkpok0[8]);
	K1P = K1P*exp(lnkpok0[9]);
	K2P = K2P*exp(lnkpok0[10]);
	K3P = K3P*exp(lnkpok0[11]);
	
	}



#------------------------------------------------------------------#
#------------------------------------------------------------------#
#                            VARIABLES                             #
#------------------------------------------------------------------#
#------------------------------------------------------------------#


	# ----------------- case 1.) PH and CO2 given
	
	PHCO2 <- function(PH, CO2, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	
	{
	
	#disp('flag = 1, pH and CO2 given');
	h=10^(-PH);
	s=CO2;
	DIC = s*(1+K1/h+K1*K2/h/h);
	HCO3 = DIC/(1+h/K1+K2/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;                        
	pCO2 = s*1.e6/Kh;
		
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;

	
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
			
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
				
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
			
			
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
			
		
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
			
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
	
	
	
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}
	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}

	
	




	# ------------ CO2 and HCO3 given ------------------
	
	CO2HCO3 <- function(CO2, HCO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 2, CO2 and HCO3 given');
	s = CO2;
	HCO3 = HCO3;
	p3 = -HCO3/K1;
	p2 = s - HCO3;
	p1 = s*K1 - HCO3*K2;
	p0 = s*K1*K2;
	p = c(p0, p1, p2, p3);
	r = polyroot(p);
	h = max(Re(r));
	pCO2=s*1.e6/Kh;
	h*1e12;
	DIC = s*(1.+K1/h+K1*K2/h/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	ALK = s*(K1/h+2*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	PH = -log10(h);
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}



	# ------------ CO2 and CO3 given ------------------
	
	CO2CO3 <- function(CO2, CO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 3, CO2 and CO3 given');
	s = CO2;
	CO3 = CO3;
	p4 = -CO3/K1/K2;
	p3 = -CO3/K2;
	p2 = s-CO3;
	p1 = s*K1;
	p0 = s*K1*K2;
	p = c(p0, p1, p2, p3, p4);    
	r = polyroot(p);
	h = max(Re(r));
	pCO2=s*1.e6/Kh;
	h*1e12;
	DIC = s*(1.+K1/h+K1*K2/h/h);
	HCO3 = DIC/(1+h/K1+K2/h);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	PH = -log10(h);
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ CO2 and ALK given ------------------

	CO2ALK <- function(CO2, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 4, CO2 and ALK given');
	s = CO2;
	ALK = ALK;
	p4 = 1.;              
	p3 = Kb+ALK;
	p2 = ALK*Kb-s*K1-Kb*bor-Kw;
	p1 = -s*Kb*K1-s*2.*K1*K2-Kw*Kb;
	p0 = -2.*s*Kb*K1*K2;
	p = c(p0, p1, p2, p3, p4);
	r = polyroot(p);
	h = max(Re(r));
	h*1e12;
	pCO2=s*1.e6/Kh;
	DIC = s*(1.+K1/h+K1*K2/h/h);
	HCO3 = DIC/(1+h/K1+K2/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	CO2=s;
	PH=-log10(h);
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}

	
	

	# ------------ CO2 and DIC given ------------------
	
	CO2DIC <- function(CO2, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 5, CO2 and DIC given');
	s = CO2;
	DIC = DIC;
	p2 = DIC - s;
	p1 = -s*K1;
	p0 = -s*K1*K2;
	p = c(p0, p1, p2);
	r = polyroot(p);
	h = max(Re(r));
	PH=-log10(h);
	h*1e12;
	pCO2=s*1.e6/Kh;
	HCO3 = DIC/(1+h/K1+K2/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}


	# ------------ PH and HCO3 given ------------------
	
	PHHCO3 <- function(PH, HCO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 6, pH and HCO3 given');
	h =10^(-PH);
	HCO3 = HCO3;
	DIC = HCO3 * (1+h/K1+K2/h);
	s = DIC / (1.+K1/h+K1*K2/h/h);
	h*1e12;
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	pCO2=s*1.e6/Kh;
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ PH and CO3 given ------------------

	PHCO3 <- function(PH, CO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 7, pH and CO3 given');
	h =10^(-PH);
	CO3 = CO3;
	DIC = CO3 * (1+h/K2+h*h/K1/K2);
	s = DIC / (1.+K1/h+K1*K2/h/h);
	h*1e12;
	HCO3 = DIC/(1+h/K1+K2/h);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	pCO2=s*1.e6/Kh;
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ PH and ALK given ------------------

	PHALK <- function(PH, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 8, pH and ALK given');
	h =10^(-PH);
	ALK = ALK;
	s = (ALK-Kw/h+h-Kb*bor/(Kb+h)) / (K1/h+2.*K1*K2/h/h);
	h*1e12;
	DIC = s*(1.+K1/h+K1*K2/h/h);
	HCO3 = DIC/(1+h/K1+K2/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	CO2=s;
	pCO2=s*1.e6/Kh;
	

	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}






	# ------------ PH and DIC given ------------------
	
	PHDIC <- function(PH, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 9, pH and DIC given');
	h =10^(-PH);
	DIC = DIC;
	s = DIC / (1.+K1/h+K1*K2/h/h);
	h*1e12;
	HCO3 = DIC/(1+h/K1+K2/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	pCO2=s*1.e6/Kh;
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ HCO3 and CO3 given ------------------

	HCO3CO3 <- function(HCO3, CO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
		
	#disp('flag = 10, HCO3 and CO3 given');
	HCO3 = HCO3;  
	CO3 = CO3;
	p3 = -CO3/K1/K2;
	p2 = -CO3/K2 + HCO3/K1;
	p1 = -CO3 + HCO3;
	p0 = HCO3*K2;
	p = c(p0, p1, p2, p3);
	r = polyroot(p);
	h = max(Re(r));
	PH=-log10(h);
	h*1e12;
	DIC = HCO3 * (1+h/K1+K2/h);
	s = DIC / (1.+K1/h+K1*K2/h/h);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	pCO2=s*1.e6/Kh;
	

	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}





	# ------------ HCO3 and ALK given ------------------
	 
	HCO3ALK <- function(HCO3, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 11, HCO3 and ALK given');
	HCO3 = HCO3;
	ALK = ALK;
	p5  = 1.;
	p4  = ALK - HCO3 + K1 + Kb;
	p3  = ALK*(Kb+K1)-HCO3*(K1+Kb+2.*K2)-Kw+K1*Kb+K1*K2-Kb*bor;
	tmp = ALK*(Kb*K1+K1*K2)-HCO3*((Kb+2.*K2)*K1+2.*Kb*K2+K1*K2);
	p2  = tmp +(-K1*Kb*bor-Kw*Kb-K1*Kw+K1*K2*Kb);
	tmp = ALK*Kb*K1*K2-HCO3*(2.*Kb*K1*K2+K2*K1*(Kb+2.*K2));
	p1  = tmp +(-K1*K2*Kb*bor-K1*Kw*Kb-K1*K2*Kw);
	p0  = -HCO3*2.*K2*Kb*K1*K2-K1*K2*Kw*Kb;
	p   = c(p0, p1, p2, p3, p4, p5);
	r = polyroot(p);
	h = max(Re(r));
	PH=-log10(h);
	h*1e12;               
	DIC = HCO3 * (1+h/K1+K2/h);
	s = DIC / (1.+K1/h+K1*K2/h/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	pCO2=s*1.e6/Kh;

	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ HCO3 and DIC given ------------------
	
	HCO3DIC <- function(HCO3, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 12, HCO3 and DIC given');
	HCO3 = HCO3;
	DIC = DIC;
	p2 = HCO3/K1;
	p1 = HCO3-DIC;    
	p0 = HCO3*K2;
	p = c(p0, p1, p2);
	r = polyroot(p);
	h = min(Re(r));       # min instead of max !!!!!
	PH=-log10(h);
	h*1e12;
	s = DIC / (1.+K1/h+K1*K2/h/h);
	DIC = s*(1.+K1/h+K1*K2/h/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	pCO2=s*1.e6/Kh;
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ CO3 and ALK given ------------------
	
	CO3ALK <- function(CO3, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 13, CO3 and ALK given');    
	CO3 = CO3;
	ALK = ALK;
	p5  = -CO3/K2+1.;
	p4  = ALK - CO3*(K1/K2+(Kb+2.*K2)/K2) + Kb + K1;
	tmp = ALK*(Kb+K1)-CO3*(K1+K1*(Kb+2.*K2)/K2+2.*Kb);
	p3  = tmp+(-Kb*bor-Kw+K1*Kb+K1*K2);
	tmp = ALK*(Kb*K1+K1*K2)-CO3*(K1*(Kb+2.*K2)+2.*Kb*K1);
	p2  = tmp+(-Kw*Kb-K1*Kb*bor-K1*Kw+K1*K2*Kb);
	tmp = ALK*Kb*K1*K2-CO3*2.*Kb*K1*K2-K1*Kw*Kb;
	p1  = tmp+(-K1*K2*Kb*bor-K1*K2*Kw);
	p0  = -K1*K2*Kw*Kb;
	p   = c(p0, p1, p2, p3, p4, p5);
	r = polyroot(p);
	h = max(Re(r));
	PH=-log10(h);
	h*1e12;
	DIC = CO3 * (1+h/K2+h^2/K1/K2);
	s = DIC / (1.+K1/h+K1*K2/h/h);
	HCO3 = DIC/(1+h/K1+K2/h);
	CO2=s;
	pCO2=s*1.e6/Kh;
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ CO3 and DIC given ------------------
	
	CO3DIC <- function(CO3, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 14, CO3 and DIC given');
	CO3 = CO3;
	DIC = DIC;
	p2 = CO3/K1/K2;
	p1 = CO3/K2;
	p0 = CO3-DIC;
	p = c(p0, p1, p2);
	r = polyroot(p);
	h = max(Re(r));
	PH=-log10(h);
	h*1e12;
	s = DIC / (1.+K1/h+K1*K2/h/h);
	HCO3 = DIC/(1+h/K1+K2/h);
	ALK = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
	CO2=s;
	pCO2=s*1.e6/Kh;
	
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}




	# ------------ ALK and DIC given ------------------
	
	ALKDIC <- function(ALK, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P)
	{
	#disp('flag = 15, ALK and DIC given');
	ALK = ALK;
	DIC = DIC;
	p5  = -1.;        
	p4  = -ALK-Kb-K1;
	p3  = DIC*K1-ALK*(Kb+K1)+Kb*bor+Kw-Kb*K1-K1*K2;
	tmp = DIC*(Kb*K1+2.*K1*K2)-ALK*(Kb*K1+K1*K2)+Kb*bor*K1;
	p2  = tmp+(Kw*Kb+Kw*K1-Kb*K1*K2);
	tmp = 2.*DIC*Kb*K1*K2-ALK*Kb*K1*K2+Kb*bor*K1*K2;
	p1  = tmp+(+Kw*Kb*K1+Kw*K1*K2);
	p0  = Kw*Kb*K1*K2;
	p   = c(p0, p1, p2, p3, p4, p5);
	r = polyroot(p);
	h = max(Re(r));
	PH=-log10(h);
	h*1e12;
	s = DIC / (1.+K1/h+K1*K2/h/h);
	HCO3 = DIC/(1+h/K1+K2/h);
	CO3 = DIC/(1+h/K2+h*h/K1/K2);
	CO2=s;
	pCO2=s*1.e6/Kh;
	
	Oa = ((0.01028*(S/35))*CO3)/Kspa;
	Oc = ((0.01028*(S/35))*CO3)/Kspc;
	
		
	B=(-1636.75+12.0408*TC-0.0327957*(TC*TC)+0.0000316528*(TC*TC*TC))^-6;
	fCO2= pCO2*exp(((P+1)*100000)*(B+2*((57.7-0.118*TC)^-6))/(8.314*TC))
	
	
	DD=-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DD)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DD)/((h+2*K2)*(h+2*K2));
	PhiD=-1/(h*log(10) * ( B+A+C ) );
	BetaD=-h*log(10)*DIC/CO2*B*PhiD;
			
				
	Q=(h+2*K2);
	V=(Kb*bor)/((h+Kb)*(h+Kb)) + Kw/(h*h)+1;         
			
					
	DB=(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DB)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DB)/((h+2*K2)*(h+2*K2));
	PhiB=-1/(h*log(10) * ( B+A+C ) );
	BetaB=-h*log(10)*DIC/CO2*B*PhiB;
				
				
	DC=2*(( K2*(2*CO3+HCO3)+ Q*V *(h+K2)+(h/K1)*( (2*CO3+HCO3)*Q+2*K2*(2*CO3+HCO3)+h*Q*V))/Q)*(1/(Q-2*(h+K2+h*h/K1)))-((-Kb*bor)/((h+Kb)*(h+Kb)))-(-Kw/((h)*(h)))+1;
	A= (2*K2*(2*CO3+HCO3)+h*(h+2*K2)*DC)/((h+2*K2)*(h+2*K2));
	B=( ( (2*CO3+HCO3) * h)/((h+2*K2)*K1) + (h/K1)* A );
	C= (-K2*(2*CO3+HCO3)+K2*(2*K2+h)*DC)/((h+2*K2)*(h+2*K2));
	PhiC=-1/(h*log(10) * ( B+A+C ) );
	BetaC=-h*log(10)*DIC/CO2*B*PhiC;
				
			
	D1=(K1*(K1*K2-h*h)*DIC)   /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D2=(-K1*K2*(2*h+K1)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	D=D1+2*D2;
	PhiH=1/ (h*log(10)* (D +(-Kb*bor/((h+Kb)*(h+Kb)))  + (-Kw/(h*h))-1))  ;   
			
		
	Pi=(h*K1*(h+2*K2)*DIC)  /  ((h*h+h*K1+K1*K2)*(h*h+h*K1+K1*K2));
	PiH=((-h/Kh)*log(10)*Pi)*PhiH;
	PiB=CO2/(Kh*DIC)*BetaB;
	PiD=CO2/(Kh*DIC)*BetaD;
	PiC=CO2/(Kh*DIC)*BetaC;
		
		
		
	if (ini=='s')
	{
	cat("Salinity:\t\t",S,"\n");
	cat("Temperature:\t\t",T,"oC\n");
	cat("Pressure:\t\t",P,"bar\n");
	cat("pH:\t\t\t",PH,"\n");
	cat("CO2:\t\t\t",CO2,"(mol/kg)\n");
	cat("pCO2:\t\t\t",pCO2,"(uatm)\n");
	cat("fCO2:\t\t\t",fCO2,"(uatm)\n");
	cat("HCO3:\t\t\t",HCO3,"(mol/kg)\n");
	cat("CO3:\t\t\t",CO3,"(mol/kg)\n");
	cat("DIC:\t\t\t",DIC,"(mol/kg)\n");
	cat("ALK:\t\t\t",ALK,"(mol/kg)\n");
	cat("Omega aragonite:\t",Oa,"\n");
	cat("Omega calcite:\t\t",Oc,"\n");
	cat("PhiD:\t\t\t",PhiD,"\n");
	cat("BetaD:\t\t\t",BetaD,"\n");
	cat("PiD:\t\t\t",PiD,"\n");
	cat("PhiB:\t\t\t",PhiB,"\n");
	cat("BetaB:\t\t\t",BetaB,"\n");
	cat("PiB:\t\t\t",PiB,"\n");
	cat("PhiC:\t\t\t",PhiC,"\n");
	cat("BetaC:\t\t\t",BetaC,"\n");
	cat("PiC:\t\t\t",PiC,"\n");
	cat("PhiH:\t\t\t",PhiH,"\n");
	cat("PiH:\t\t\t",PiH,"\n");


	date=date();
		col <- c("date", "flag", "Salinity", "Temperature", "Pressure", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite", "PhiD", "BetaD", "PiD", "PhiB", "BetaB", "PiB", "PhiC", "BetaC", "PiC", "PhiH", "PiH")
	res<-data.frame(date,flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	if (file.exists("carb.out")==TRUE)
		{
		write.table(res,resultat,sep="	",append=TRUE,col.name=FALSE,row.names=FALSE);
		}
		else
		{
		colnames(res) <- col;
		write.table(res,resultat,sep="	",append=TRUE,row.names=FALSE);
		}
	}

	
	if (ini=='f')
	
	{
	Salinity=S;
	Temperature=paste(T,"oC");
	Pressure=paste(P,"bar");
	pH=PH;
	CO2=paste(CO2,"(mol/kg)");
	pCO2=paste(pCO2,"(uatm)");
	fCO2=paste(fCO2,"(uatm)");
	HCO3=paste(HCO3,"(mol/kg)");
	CO3=paste(CO3,"(mol/kg)");
	DIC=paste(DIC,"(mol/kg)");
	ALK=paste(ALK,"(mol/kg)");
	OmegaAragonite=Oa;
	OmegaCalcite=Oc;
	
	res<-data.frame(Salinity,Temperature,Pressure,pH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,OmegaAragonite,OmegaCalcite,PhiD,BetaD,PiD,PhiB,BetaB,PiB,PhiC,BetaC,PiC,PhiH,PiH);
	
	
	if(i==1)
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE);
	}
	else
	{
	write.table(res,resu,sep="	",append=TRUE,row.names=FALSE,col.names=FALSE);
	}
	
	}
	
	}



	if (flag==1)
	{
	PHCO2(PH, CO2, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==2)
	{
	CO2HCO3(CO2, HCO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==3)
	{
	CO2CO3(CO2, CO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==4)
	{
	CO2ALK(CO2, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==5)
	{
	CO2DIC(CO2, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==6)
	{
	PHHCO3(PH, HCO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==7)
	{
	PHCO3(PH, CO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==8)
	{
	PHALK(PH, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==9)
	{
	PHDIC(PH, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==10)
	{
	HCO3CO3(HCO3, CO3, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==11)
	{
	HCO3ALK(HCO3, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==12)
	{
	HCO3DIC(HCO3, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==13)
	{
	CO3ALK(CO3, ALK, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}

	if (flag==14)
	{
	CO3DIC(CO3, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	if (flag==15)
	{
	ALKDIC(ALK, DIC, K1, K2, bor, Kb, Kw, Kh, Kspa, Kspc, S, TC, P);
	}
	
	
	# ======================================================
	
	

	
	
	}

}
