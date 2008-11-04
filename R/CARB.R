# Copyright (C) 2008 Jean-Pierre Gattuso and Héloïse Lavigne and Aurélien Proye
# with a most valuable contribution of Bernard Gentili <gentili@obs-vlfr.fr>
# and valuable suggestions from Jean-Marie Epitalon <epitalon@lsce.saclay.cea.fr>
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
#
carb<-
function(flag, var1, var2, S=35, T=25, P=0, Pt=0, Sit=0){
RES <- data.frame()
n <- max(length(var1), length(var2), length(S), length(T), length(P), length(Pt), length(Sit))
if(length(flag)!=n){ flag <- rep(flag[1],n)}
if(length(var1)!=n){ var1 <- rep(var1[1],n)}
if(length(var2)!=n){ var2 <- rep(var2[1],n)}
if(length(S)!=n){ S <- rep(S[1],n)}
if(length(T)!=n){ T <- rep(T[1],n)}
if(length(P)!=n){ P <- rep(P[1],n)}
if(length(Pt)!=n){ Pt <- rep(Pt[1],n)}
if(length(Sit)!=n){ Sit <- rep(Sit[1],n)}
df <- data.frame(flag, var1, var2, S, T, P, Pt, Sit)
for(i in (1:nrow(df))) {
 flag <- as.numeric(df[i,1])
  var1 <- as.numeric(df[i,2])
  var2 <- as.numeric(df[i,3])
  S <- as.numeric(df[i,4])
  T <- as.numeric(df[i,5])
  P <- as.numeric(df[i,6])
  Pt <- as.numeric(df[i,7])
  Sit <- as.numeric(df[i,8])

	
#-------Constantes----------------

tk = 273.15;           # [K] (for conversion [deg C] <-> [K])

# JME: moved following code block here, after reading imput file

TK = T + tk;           # TK [K]; T[C]

#---- issues de equic----
Cl = S / 1.80655;      # Cl = chlorinity; S = salinity (per mille)
cl3 = Cl^(1/3);
ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl * Cl;   # ionic strength
iom0 = 19.924*S/(1000-1.005*S);
ST = 0.14/96.062/1.80655*S;   # (mol/kg soln) total sulfate
bor = (416.*(S/35.))* 1e-6;   # (mol/kg), DOE94 boron total
fluo = (7*(S/35))*1e-5        # (mol/kg), DOE94 fluoride total

#---------------------------------------------------------------------
#--------------------- calcul des K ----------------------------------
#---------------------------------------------------------------------

K1 <- K1(S=S, T=T, P=P)
K2 <- K2(S=S, T=T, P=P)
Ks <- Ks(S=S, T=T, P=P)
Kf <- Kf(S=S, T=T, P=P)
Kw <- Kw(S=S, T=T, P=P)
Kh <- Kh(S=S, T=T, P=P)
Kb <- Kb(S=S, T=T, P=P)
K1p <- K1p(S=S, T=T, P=P)
K2p <- K2p(S=S, T=T, P=P)
K3p <- K3p(S=S, T=T, P=P)
Ksi <- Ksi(S=S, T=T, P=P)
Kspa <- Kspa(S=S, T=T, P=P)
Kspc <- Kspc(S=S, T=T, P=P)

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
	  lnkpok0[ipc]  = -(deltav[ipc] /(R*TK))*P + (0.5*deltak[ipc] /(R*TK))*P*P;
	}

	K1 = K1*exp(lnkpok0[1]);
	K2 = K2*exp(lnkpok0[2]);
	Kb = Kb*exp(lnkpok0[3]);
	Kw = Kw*exp(lnkpok0[4]);
	Ks = Ks*exp(lnkpok0[5]);
	Kf = Kf*exp(lnkpok0[6]);
	Kspc = Kspc*exp(lnkpok0[7]);
	Kspa = Kspa*exp(lnkpok0[8]);
	K1p = K1p*exp(lnkpok0[9]);
	K2p = K2p*exp(lnkpok0[10]);
	K3p = K3p*exp(lnkpok0[11]);
	}

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

#------------------------------------------------------------------#
#------------------------------------------------------------------#
#                            VARIABLES                             #
#------------------------------------------------------------------#
#------------------------------------------------------------------#

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
# flag = 21     pCO2-pH given
# flag = 22     pCO2-HCO3 given
# flag = 23     pCO2-CO3 given
# flag = 24     pCO2-ALK given
# flag = 25     pCO2-DIC given


	# ------------ case 1.) PH and CO2 given
	if (flag==1)
	{
	PH <- var1
	CO2 <- var2
	h <- 10^(-PH)
	fCO2 <- CO2/Kh
	HCO3 <- (K1*CO2)/h
	CO3 <- (K2*HCO3)/h
	DIC <- CO2 + HCO3 + CO3
	}

	# ------------ case 2.) CO2 and HCO3 given 
	if (flag==2)
	{
	CO2 <- var1
	HCO3 <- var2
	fCO2 <- CO2/Kh
	h <- Kh*K1*fCO2/HCO3
	CO3 <- Kh*K1*K2*fCO2/(h*h)
	DIC <- CO2 + HCO3 + CO3
	PH <- -log10(h)
	}

	# ------------ case 3.) CO2 and CO3 given
	if (flag==3)
	{
	CO2 <- var1
	CO3 <- var2
	fCO2 <- CO2/Kh
	h <- sqrt((Kh*K1*K2*fCO2)/CO3)
	HCO3 <- (Kh*K1*fCO2)/h
	DIC <- CO2 + HCO3 + CO3
	PH <- -log10(h)
	}

	# ------------ case 4.) CO2 and ALK given
	Ks <- Ks(S=S, T=T, P=P)
	if (flag==4)
	{
	CO2 <- var1
	ALK <- var2
	fALK <- function(x)# K1=K1, K2=K2, CO2=CO2, bor=bor, Kb=Kb, Kw=Kw, Pt=Pt, K1p=K1p, K2p=K2p, K3p=K3p, Sit=Sit, Ksi=Ksi, NH3t=NH3t, KNH3=KNH3, H2St=H2St, KH2S=KH2S, ST=ST, Ks=Ks, fluo=fluo, Kf=Kf, ALK=ALK) {
	# composants for ALK
	{DIC <- CO2*(1+K1/x+K1*K2/(x*x))
	hco3 <- DIC*x*K1/(x*x + K1*x + K1*K2)
	co3 <- DIC*K1*K2/(x*x + K1*x + K1*K2)
	boh4 <- bor/(1+x/Kb)
	oh <- Kw/x
	h3po4 <- Pt*x^3/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	hpo4 <- Pt*K1p*K2p*x/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	po4 <- Pt*K1p*K2p*K3p/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	siooh3 <- Sit/(1+x/Ksi)
	h <- x/(1+ST/Ks)
	hso4 <- ST/(1+Ks/h)
	hf <- fluo/(1+Kf/x)
	############
	OUT <- hco3+2*co3+boh4+oh+hpo4+2*po4+siooh3-h-hso4-hf-h3po4-ALK
	OUT}	
	h <- uniroot(fALK,c(1e-12,10^(-3.5)), tol=1e-20)$root
	DIC <- CO2*(1+K1/h+K1*K2/(h*h))
	HCO3 <- (DIC*K1*h)/(h*h+K1*h+K1*K2)
	CO3 <- (DIC*K1*K2)/(h*h+K1*h+K1*K2)
	fCO2 <- CO2/Kh
	PH <- -log10(h)
	}

	# ------------ case 5.) CO2 and DIC given
	if (flag==5)
	{
	CO2 <- var1
	DIC <- var2
	fCO2 <- CO2/Kh
	a <- K1*K2*CO2
	b <- K1*CO2
	c <- CO2 - DIC
	D <- b*b - 4*a*c
	X <- (sqrt(D)-b)/(2*a)  # X = 1/h
	h <- 2*K1*K2*CO2/(sqrt(K1*CO2*K1*CO2 - 4*K1*K2*CO2*(CO2 - DIC))-K1*CO2)
	HCO3 <- Kh*K1*fCO2/h
	CO3 <- DIC - CO2 - HCO3
	PH <- -log10(h)
	}

	# ------------ case 6.) PH and HCO3 given
	if (flag==6)
	{
	PH <- var1
	HCO3 <- var2
	h <- 10^(-PH)
	CO2 <- (HCO3*h)/K1
	CO3 <- K2*HCO3/h
	DIC <- CO2 + HCO3 + CO3
	fCO2 <- CO2/Kh
	}

	# ------------ case 7.) PH and CO3 given	
	if (flag==7)
	{
	PH <- var1
	CO3 <- var2
	h <- 10^(-PH)
	HCO3 <- CO3*h/K2
	CO2 <- HCO3*h/K1
	fCO2 <- CO2/Kh
	DIC <- CO2 + HCO3 + CO3
	}

	# ------------ case 8.) PH and ALK given
	if (flag==8)
	{
	PH <- var1
	ALK <- var2 
	h <- 10^(-PH)
	fALK <- function(x) #h=h, K1=K1, K2=K2, x, bor=bor, Kb=Kb, Kw=Kw, Pt=Pt, K1p=K1p, K2p=K2p, K3p=K3p, Sit=Sit, Ksi=Ksi, NH3t=NH3t, KNH3=KNH3, H2St=H2St, KH2S=KH2S, ST=ST, Ks=Ks, fluo=fluo, Kf=Kf, ALK=ALK) {
	# composants for ALK
	{hco3 <- x*h*K1/(h*h + K1*h + K1*K2)
	co3 <- x*K1*K2/(h*h + K1*h + K1*K2)
	boh4 <- bor/(1+h/Kb)
	oh <- Kw/h
	h3po4 <- Pt*h^3/(h^3+K1p*h^2+K1p*K2p*h+K1p*K2p*K3p)
	hpo4 <- Pt*K1p*K2p*h/(h^3+K1p*h^2+K1p*K2p*h+K1p*K2p*K3p)
	po4 <- Pt*K1p*K2p*K3p/(h^3+K1p*h^2+K1p*K2p*h+K1p*K2p*K3p)
	siooh3 <- Sit/(1+h/Ksi)
	hff <- h/(1+ST/Ks)
	hso4 <- ST/(1+Ks/hff)
	hf <- fluo/(1+Kf/h)
	############
	OUT <- hco3+2*co3+boh4+oh+hpo4+2*po4+siooh3-h-hso4-hf-h3po4-ALK
	OUT}	
	DIC <- uniroot(fALK,c(1e-3,1e-2), tol=1e-20)$root
	CO2 <- DIC/(1+K1/h+K1*K2/(h^2))
	HCO3 <- CO2*K1/h
	CO3 <- HCO3*K2/h
	fCO2 <- CO2/Kh
	}
	
	# ------------ case 9.) PH and DIC given
	if (flag==9)
	{
	PH <- var1
	DIC <- var2
	h <- 10^(-PH)
	HCO3 <- (DIC*K1*h)/(h*h+K1*h+K1*K2)
	CO3 <- (DIC*K1*K2)/(h*h+K1*h+K1*K2)
	CO2 <- h*HCO3/K1
	fCO2 <- CO2/Kh
	}
	
	# ------------ case 10.) HCO3 and CO3 given	
	if (flag==10)
	{
	HCO3 <- var1
	CO3 <- var2
	h <- K2*HCO3/CO3
	CO2 <- h*HCO3/K1
	DIC <- CO2 + HCO3 + CO3
	fCO2 <- CO2/Kh
	PH <- -log10(h)
	}

	# ------------ case 11.) HCO3 and ALK given
	if (flag==11)
	{
	HCO3 <- var1
	ALK <- var2
	
	fALK <- function(x)# K1=K1, K2=K2, HCO3=HCO3, bor=bor, Kb=Kb, Kw=Kw, Pt=Pt, K1p=K1p, K2p=K2p, K3p=K3p, Sit=Sit, Ksi=Ksi, NH3t=NH3t, KNH3=KNH3, H2St=H2St, KH2S=KH2S, ST=ST, Ks=Ks, fluo=fluo, Kf=Kf, ALK=ALK) {
	# composants for ALK
	{DIC <- HCO3*(x^2+K1*x+K1*K2)/(K1*x)
	hco3 <- HCO3
	co3 <- DIC*K1*K2/(x*x + K1*x + K1*K2)
	boh4 <- bor/(1+x/Kb)
	oh <- Kw/x
	h3po4 <- Pt*x^3/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	hpo4 <- Pt*K1p*K2p*x/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	po4 <- Pt*K1p*K2p*K3p/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	siooh3 <- Sit/(1+x/Ksi)
	h <- x/(1+ST/Ks)
	hso4 <- ST/(1+Ks/h)
	hf <- fluo/(1+Kf/x)
	############
	OUT <- hco3+2*co3+boh4+oh+hpo4+2*po4+siooh3-h-hso4-hf-h3po4-ALK
	OUT}

	h <- uniroot(fALK,c(1e-10,10^(-3.5)),tol=1e-20)$root
	CO2 <- h*HCO3/K1
	CO3 <- K2*HCO3/h
	DIC <- CO2 + HCO3 + CO3
	PH <- -log10(h)
	fCO2 <- CO2/Kh
	}

	# ------------ case 12.) HCO3 and DIC given
	if (flag==12)
	{
	HCO3 <- var1
	DIC <- var2
	a <- HCO3
	b <- K1*(HCO3-DIC)
	c <- K1*K2*HCO3
	D <- b*b - 4*a*c
	h <- (-b-sqrt(D))/(2*a)
	CO2 <- h*HCO3/K1
	CO3 <- K2*HCO3/h
	fCO2 <- CO2/Kh
	PH <- -log10(h)
	}

	# ------------ case 13.) CO3 and ALK given
	if (flag==13)
	{
	CO3 <- var1
	ALK <- var2
	
	fALK <- function(x)# K1=K1, K2=K2, HCO3=HCO3, bor=bor, Kb=Kb, Kw=Kw, Pt=Pt, K1p=K1p, K2p=K2p, K3p=K3p, Sit=Sit, Ksi=Ksi, NH3t=NH3t, KNH3=KNH3, H2St=H2St, KH2S=KH2S, ST=ST, Ks=Ks, fluo=fluo, Kf=Kf, ALK=ALK) {
	# composants for ALK
	{DIC <- CO3*(x^2+K1*x+K1*K2)/(K1*K2)
	hco3 <- DIC*K1*x/(x*x + K1*x + K1*K2)
	co3 <- CO3
	boh4 <- bor/(1+x/Kb)
	oh <- Kw/x
	h3po4 <- Pt*x^3/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	hpo4 <- Pt*K1p*K2p*x/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	po4 <- Pt*K1p*K2p*K3p/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	siooh3 <- Sit/(1+x/Ksi)
	h <- x/(1+ST/Ks)
	hso4 <- ST/(1+Ks/h)
	hf <- fluo/(1+Kf/x)
	############
	OUT <- hco3+2*co3+boh4+oh+hpo4+2*po4+siooh3-h-hso4-hf-h3po4-ALK
	OUT}
	h <- uniroot(fALK,c(1e-10,10^(-3.5)),tol=1e-20)$root
	HCO3 <- h*CO3/K2
	CO2 <- h*HCO3/K1
	fCO2 <- CO2/Kh
	DIC <- HCO3+CO2+CO3
	PH <- -log10(h)
	}

	# ------------ case 14.) CO3 and DIC given
	if (flag==14)
	{
	CO3 <- var1
	DIC <- var2
	h <- (-K1*CO3 + sqrt(((K1*CO3)^2)-4*CO3*K1*K2*(CO3-DIC)))/(2*CO3)
	HCO3 <- h*CO3/K2
	CO2 <- h*HCO3/K1
	fCO2 <- CO2/Kh
	PH <- -log10(h)
	}
	
	# ------------ case 15.) ALK and DIC given
	if (flag==15)
	{
	ALK <- var1
	DIC <- var2

	fALK <- function(x) # K1=K1, K2=K2, DIC=DIC, bor=bor, Kb=Kb, Kw=Kw, Pt=Pt, K1p=K1p, K2p=K2p, K3p=K3p, Sit=Sit, Ksi=Ksi, NH3t=NH3t, KNH3=KNH3, H2St=H2St, KH2S=KH2S, ST=ST, Ks=Ks, fluo=fluo, Kf=Kf, ALK=ALK) {
	# composants for ALK
	{hco3 <- DIC*x*K1/(x*x + K1*x + K1*K2)
	co3 <- DIC*K1*K2/(x*x + K1*x + K1*K2)
	boh4 <- bor/(1+x/Kb)
	oh <- Kw/x
	h3po4 <- Pt*x^3/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	hpo4 <- Pt*K1p*K2p*x/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	po4 <- Pt*K1p*K2p*K3p/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	siooh3 <- Sit/(1+x/Ksi)
	h <- x/(1+ST/Ks)
	hso4 <- ST/(1+Ks/h)
	hf <- fluo/(1+Kf/x)
	############
	OUT <- hco3+2*co3+boh4+oh+hpo4+2*po4+siooh3-h-hso4-hf-h3po4-ALK
	OUT}

	h <- uniroot(fALK,c(1e-10,10^(-3.5)),tol=1e-30)$root
	
	HCO3 <- (DIC*K1*h)/(h*h+K1*h+K1*K2)
	CO3 <- (DIC*K1*K2)/(h*h+K1*h+K1*K2)
	CO2 <- h*HCO3/K1
	fCO2 <- CO2/Kh
	PH <- -log10(h)
	}

	# ------------ calculation of pCO2 for cases 1 to 15 
	# JME: corrected fugacity calculation
	if ((flag>=1)&(flag<=15))
	{
	B=(-1636.75+12.0408*TK-0.0327957*(TK*TK)+0.0000316528*(TK*TK*TK))*1e-6;
	pCO2= fCO2*(1/exp(((P+1)*100000)*(B+2*(57.7-0.118*TK)*1e-6)/(8.314*TK)))
	}

	# ------------ calculation of fCO2 for cases 21 to 25
	# JME: corrected fugacity calculation
	if ((flag>=21)&(flag<=25))
	{
	pCO2 <- var1*1e-6
	B=(-1636.75+12.0408*TK-0.0327957*(TK*TK)+0.0000316528*(TK*TK*TK))*1e-6;
	fCO2= pCO2*(exp(((P+1)*100000)*(B+2*(57.7-0.118*TK)*1e-6)/(8.314*TK)))
	}

	# ------------ case 21.) PH and pCO2 given
	if (flag==21)
	{
	PH <- var2
	h <- 10^(-PH)
	CO2 <- Kh*fCO2
	HCO3 <- K1*CO2/h
	CO3 <- K2*HCO3/h
	DIC <- CO2 + HCO3 + CO3
	}

	# ------------ case 22.) HCO3 and pCO2 given
	if (flag==22)
	{
	HCO3 <- var2
	CO2 <- fCO2*Kh
	h <- CO2*K1/HCO3
	CO3 <- HCO3*K2/h
	DIC <- CO2 + HCO3 + CO3
	PH <- -log10(h)
	}

	# ------------ case 23.) CO3 and pCO2 given
	if (flag==23)
	{
	CO3 <- var2
	h <- sqrt(Kh*K1*K2*fCO2/CO3)
	HCO3 <- h*CO3/K2
	CO2 <- h*HCO3/K1
	DIC <- CO2 + HCO3 + CO3
	PH <- -log10(h)
	}

	# ------------ case 24.) ALK and pCO2 given
	if (flag==24)
	{
	ALK <- var2
	CO2 <- fCO2*Kh

	fALK <- function(x)# K1=K1, K2=K2, CO2=CO2, bor=bor, Kb=Kb, Kw=Kw, Pt=Pt, K1p=K1p, K2p=K2p, K3p=K3p, Sit=Sit, Ksi=Ksi, NH3t=NH3t, KNH3=KNH3, H2St=H2St, KH2S=KH2S, ST=ST, Ks=Ks, fluo=fluo, Kf=Kf, ALK=ALK) {
	# composants for ALK
	{DIC <- CO2*(1+K1/x+K1*K2/(x*x))
	hco3 <- DIC*x*K1/(x*x + K1*x + K1*K2)
	co3 <- DIC*K1*K2/(x*x + K1*x + K1*K2)
	boh4 <- bor/(1+x/Kb)
	oh <- Kw/x
	h3po4 <- Pt*x^3/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	hpo4 <- Pt*K1p*K2p*x/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	po4 <- Pt*K1p*K2p*K3p/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	siooh3 <- Sit/(1+x/Ksi)
	h <- x/(1+ST/Ks)
	hso4 <- ST/(1+Ks/h)
	hf <- fluo/(1+Kf/x)
	############
	OUT <- hco3+2*co3+boh4+oh+hpo4+2*po4+siooh3-h-hso4-hf-h3po4-ALK
	OUT}	
	h <- uniroot(fALK,c(1e-10,10^(-3.5)),tol=1e-20)$root
	HCO3 <- K1*CO2/h
	CO3 <- K2*HCO3/h
	PH <- -log10(h)
	DIC <- CO2 + HCO3 + CO3
	}

	# ------------ case 25.) DIC and pCO2 given
	if (flag==25)
	{
	DIC <- var2
	CO2 <- Kh*fCO2
	K <- K1/K2
	HCO3 <- (1/2)*(-K*Kh*fCO2+sqrt((K*Kh*fCO2)^2 - 4*(K*Kh*fCO2)*(Kh*fCO2-DIC)))
	CO3 <- DIC - CO2 - HCO3
	h <- K1*CO2/HCO3
	PH <- -log10(h)
	}

	# ------------ CALCULATION OF ALK in cases 
	Ks <- Ks(S=S, T=T, P=P)
	cases <- c(1, 2, 3, 5, 6, 7, 9, 10, 12, 14, 21, 22, 23, 24, 25)
	if (length(which(cases==flag))==1){
	x <- h
	hco3 <- DIC*x*K1/(x*x + K1*x + K1*K2)
	co3 <- DIC*K1*K2/(x*x + K1*x + K1*K2)
	boh4 <- bor/(1+x/Kb)
	oh <- Kw/x
	h3po4 <- Pt*(x^3)/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	hpo4 <- Pt*K1p*K2p*x/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	po4 <- Pt*K1p*K2p*K3p/(x^3+K1p*x^2+K1p*K2p*x+K1p*K2p*K3p)
	siooh3 <- Sit/(1+x/Ksi)
	h <- x/(1+ST/Ks)
	hso4 <- ST/(1+Ks/h)
	hf <- fluo/(1+Kf/x)
	ALK <- hco3+2*co3+boh4+oh+hpo4+2*po4+siooh3-h-hso4-hf-h3po4
	}

	##########################################################
	# CALCULATION OF ARAGONITE AND CALCITE SATURATION STATE  #
	##########################################################

	Oa  <- ((0.01028*(S/35))*CO3)/Kspa
	Oc  <- ((0.01028*(S/35))*CO3)/Kspc
	
	#PCO2 and fCO2 converted in microatmosphere
	pCO2 <- pCO2*1e6
	fCO2 <- fCO2*1e6	

	res<-data.frame(flag,S,T,P,PH,CO2,pCO2,fCO2,HCO3,CO3,DIC,ALK,Oa,Oc)
	RES<- rbind(RES, res)
	}
	names(RES) <- c("flag", "S", "T", "P", "pH", "CO2", "pCO2", "fCO2", "HCO3", "CO3", "DIC", "ALK", "OmegaAragonite", "OmegaCalcite")
	return(RES)
}


