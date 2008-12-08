# Copyright (C) 2008 Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# sys is 0 if the system is closed and 1 if it is open
# co3 and hco3 are the   amount added in mol/kg


# One adds a certain volume (vol) of HCl (vol  is then negative) or NaOH (vol is then positive) to 1 kg of seawater. 
# vol is in liter
# N is the normality of the HCl or NaOH
# HCl or NaOH are, respectively, a strong acid and a strong base which are therefore fully dissociated. The concentration of H+ or OH- is N mole/l
# HCl 1N: 0.1 ml therefore adds 0.1 e-3 mol of H+, decreasing TA
# NaOH 1N: 0.1 ml therefore adds 0.1 e-3 mol of OH-, increasing TA
# DIC is constant in a closed system
"ppH" <-
function(flag, sys, var1, var2, pCO2a, vol, N, S=35, T=20, P=0, Pt=0, Sit=0){
	if (sys==0) {
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit)
		alkf <- ci$ALK + vol*N  # final alkalinity
		dicf <- ci$DIC	# final dic
		cf <- carb(flag=15, var1=alkf, var2=dicf, S=S, T=T, P=P)
		co <- as.data.frame(c("ppH-closed-initial", rep("ppH-closed-final", nrow(cf)))) 
	}
	if (sys==1) {
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit)
		alkf <- ci$ALK + vol*N # final total alkalinity
		dicf <- ci$DIC	# final dic  before requilibration
		#pHf <- ci$pH + ci$PhiH * (-vol) *N	# final pH using a buffer factor (see Frankignoulle, 1994)
		cc <- carb(flag=15, var1=alkf, var2=dicf, S=S, T=T, P=P) #  before requilibration	
		cf <- carb(flag=24,var1=pCO2a, var2=alkf, S=S, T=T, P=P)		
		co <- as.data.frame(c("ppH-open-initial", rep("ppH-open-final", nrow(cf))))
	}
	out <- rbind(ci, cf)
	out <- cbind(co, out)
	names(out)[1] <- "comment"
	return(out[,1:15])
}

## Example:
## ppH(flag=15, sys=0, var1i=2302e-6, var2i=2050e-6, pCO2a=444.8911, vol=-0.13e-3, N=1, S=35, T=20, P=0, Pt=0, Sit=0)
## Debugging:
#flag <- 15
#sys=0
#var1i=2302e-6
#var2i=2050e-6
#pCO2a=387
#vol=0.13e-3
#N=1
#S=35
#T=20
#P=0
#ci <- carb(data.frame(flag, var1i, var2i, S ,T, P, k1k2='r', phflag=0))
#alkf <- ci$ALK + vol*N  # final alkalinity
#dicf <- ci$DIC	# final dic
#cf <- carb(data.frame(flag=15,alkf, dicf, S, T, P, k1k2='r', phflag=0))		
#co <- as.data.frame(c("ppH-closed-initial", rep("ppH-final-closed", nrow(cf))))
#print(ci)
#print(cf)
#print(co)
#out <- rbind(ci, cf)
#out <- cbind(co, out)
#names(out)[1] <- "comment"
#print(out[,1:15])



#########################################################################
## CLOSED SYSTEM
########################################################################
## One adds 0.1 ml of HCl 1N to 1 kg of seawater. It is a strong acid and is therefore fully dissociated. The concentration of H+ is 1 mole/l
## 0.1 ml therefore adds 0.1 e-3 mol of H+, decreasing TA
## DIC is constant as it is a closed system
##Carbonate chemistry of the initial sea water
#
#rm(out)
#ta <- 2302e-6
#dic <- 2050e-6
#S <- 35 # salinity
#T <- 20 # temperature in ¡C
#c <- carb(data.frame(flag=15,var1=ta,var2=dic,S,T,P=0,k1k2='r',phflag=0))
#out <- c
#rownames(out)[nrow(out)] <- "closed initial"
#acid <-  0.13 # volume of acid added in ml
#dta<- -acid * 1e-3 #mol/kg
#acid_ta <- ta + dta
#acid_dic <- dic
#c_acid <- carb(data.frame(flag=15,var1=acid_ta,var2=acid_dic,S,T,P=0,k1k2='r',phflag=0)) #Now one recalculates the carbonate chemistry with the new ta
#out <- rbind(out, c_acid)
#rownames(out)[nrow(out)] <- "closed final"
#
#
#
########################################################################
## OPEN SYSTEM
########################################################################
##Carbonate chemistry of the initial sea water
#c <- carb(data.frame(flag=15,var1=ta,var2=dic,S,T,P=0,k1k2='r',phflag=0))
#out <- rbind(out, c)
#rownames(out)[nrow(out)] <- "Open Phi initial"
## One adds 0.1 ml of HCl 1N to 1 kg of seawater. It is a strong acid and is therefore fully dissociated. The concentration of H+ is 1 mole/l
## 0.1 ml therefore adds 0.1 e-3 mol of H+, decreasing TA
#acid <-  0.13 # volume of acid added in ml
#dta<- -acid * 1e-3 #mol/kg
#c$pH
#c$PhiH
#acid_pH <- c$pH + c$PhiH * abs(dta)
#acid_ta <- ta + dta
#acid_pH
#acid_ta
#z <- carb(data.frame(flag=8,var1=acid_pH,var2=acid_ta,S,T,P=0,k1k2='r',phflag=0)) #Now one recalculates the carbonate chemistry with the new pH and ALK
#out <- rbind(out, z)
#rownames(out)[nrow(out)] <- "Open Phi bef. equil."
#c_acid_open <- carb(data.frame(flag=24,var1=c$pCO2,var2=z$ALK,S,T,P=0,k1k2='r',phflag=0)) 
#out <- rbind(out, c_acid_open)
#rownames(out)[nrow(out)] <- "Open Phi aft. equilibration"
#
## Compare with DIC and ALK
#c <- carb(data.frame(flag=15,var1=ta,var2=dic,S,T,P=0,k1k2='r',phflag=0))
#out <- rbind(out, c)
#rownames(out)[nrow(out)] <- "Open DIC+TA initial"
#z <- carb(data.frame(flag=15,var1=acid_ta,var2=dic,S,T,P=0,k1k2='r',phflag=0)) # before equilibration
#out <- rbind(out, z)
#rownames(out)[nrow(out)] <- "Open DIC+ALK bef. equil."
#c <- carb(data.frame(flag=24,var1=c$pCO2, var2=acid_ta,S,T,P=0,k1k2='r',phflag=0)) # before equilibration
#out <- rbind(out, c)
#rownames(out)[nrow(out)] <- "Open DIC+ALK aft. equil."
## Note that the simple approach with DIC and ALK just above gives the same results as the
## much more complicated Open-Phi method above it. It is therefore used in the ppH function
#
## Compare with gas bubbling to maintain a pCO2 similar as the one in a closed system
#c_bubbling <- carb(data.frame(flag=24,var1=out["closed final", "pCO2"], var2=ta,S,T,P=0,k1k2='r',phflag=0))
#out <- rbind(out, c_bubbling)
#rownames(out)[nrow(out)] <- "gas bubbling"
#
## Now one wants to replenish the bicarbonate after +acid (closed)
#dHCO3 <- c$HCO3 - c_acid$HCO3 # amount to replenish
#HCO3_pH <- acid_pH + c_acid$PhiB * dHCO3 #new pH
#c_HCO3 <- carb(data.frame(flag=6,var1=HCO3_pH,var2=c$HCO3,S,T,P=0,k1k2='r',phflag=0)) #Now one recalculates the carbonate chemistry with the new pH and the desired bicarbonate concentration
#
#print(out[,1:14])
#