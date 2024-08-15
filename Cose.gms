$Title New
$Offsymxref
$Offsymlist
$inlinecom [ ]
$eolcom //
$include  Sets
$include  Parameters

scalars
IC_SF             /120/
IC_TES            /120/
IC_PB             /910/

IC_C              /440759/
IC_E              /985490/
IC_LPT            /930/
IC_HPT            /121096/
IC_CD             /170000/

OM_SF             /30/
OM_TES            /30/
OM_PB             /110/

OM_E              /19710/
OM_C              /35260/
OM_LPT            /19/
OM_HPT            /2422/
OM_CD             /5100/
CRF               /0.3/
epsilon           /0.5/
;

******* modelling*************
$onText
Parameters D, m_h2, E_e, P_e, m_h2, C_l, C_h, m_p, P_c, LHV, E_c, E_p, CRF1, CRF2, CRF ;
positive variables  A_sf, E_g(t), E_exc(t) ;
variable    cost ;
Equations   obf1 ;

*Daily hydrogen demand:
      D = Cj + Nf + d;

*CRF calculatation
      CRF1 =  (Ir-(DeltaIFF*Iff))/(1+(DeltaIFF* Iff));
      CRF2 =  (1+CRF1)**nn ;
      CRF  =  CRF1*CRF2 /(CRF2 -1)  ;

*electrolyer modelling:
      m_h2 = D/td ;
      P_e  = ((m_h2 )* LHV)/ mu_e;
      E_e  = P_e * td;

*High and low pressure tanks modelling:
      C_l   = m_h2 * tdis;
      C_h   = mp * SF;

*compressor modelling:
      P_c  = ((C_p * T)/ mu_c) * ((P1/P2)^((r-1)/r)-1) * (m_h2);
      E_c  = (P_c * t_d)/ mu_e;

*precooler and dispenser
      E_p  = (0.3/(1.6 * e^(-0.018*T_a))) + ((25*LN*(T_a)-21)/D);


******** CSP POWER *****************

parameter P(t), Q_u(t), CF, Ecsp(t) ;
P(t)    = mu_SF  * DNI(t) * A_sf ;
Q_u(t)  = P(t)   * mu_rec * mu_thf * mu_tes ;
Ecsp(t) = Q_u(t )* mu_PB  * (1-mu_par) * CF ;
$offText

Parameter
DeltaI /1/
DeltaIFF /5/
DeltaD /1/;


Scalar
mu_ACAC /0.95/
mu_sf  /0.9/
mu_rec /0.9/
mu_thf /0.9/
mu_tes /0.9/
mu_pb  /0.9/
mu_par /0.9/
mu_ac  /0.9/
CF     /0.45/
;
*parameters Pt , qu, CF, Ecsp;
parameters Ecsp(t);

*Pt     = mu_sf  * DNI * A_sf ;
*qu     = Pt    * mu_rec * mu_thf * mu_tes ;
*Ecsp   = qu    * mu_pb  * (1-mu_par) * CF ;
*mu_rec  = mu_rec* 1 ; // for sensitivity analysis
Ecsp(t)   = mu_sf * DNI(t) * mu_rec    * mu_pb  * CF/1000;



*parameter sum12;
*sum12=sum(t,Ecsp(t));
*display sum12;
*scalars
*Ee  /9352.6 /
*Ec  /982.33/
*Ep  /0.404 /

*;
*parameter E_L;
*E_L  = (Ee + Ec + Ep)*8750;
*display E_L(t);
*E_L    = E_L * 1; // use this for sensitivity analysis

******* OPTIMIZATION *************
Positive variable A_sf, E_g(t) ;
Positive variable E_exc(t) ;

variable Z ;
Equations obj, con1(t)  ;
obj.. Z =e=  (0.085 * sum(t,E_g(t)) - 0.0231 * sum(t,E_exc(t)) )
           + (( A_Sf*IC_sf * CRF +  IC_PB * sum(t,Ecsp(t))*CRF ) + IC_E + IC_C + IC_HPT + IC_LPT + IC_CD )
           + (  A_sf*OM_sf       +  OM_PB * sum(t,Ecsp(t)) ) + OM_E + OM_C + OM_HPT + OM_LPT + OM_CD ;

con1(t)..  A_sf * Ecsp(t) + E_g(t) - E_exc(t) =e= D(t) ;

Equations con2 ;
con2.. (sum(t,D(t)) - A_sf * sum (t, Ecsp(t)) * mu_ac)/(sum(t,D(t))) =l= epsilon ;

Equation  con3(t) ;
 con3(t)   .. E_g(t)=l= D(t);

*Equation SELL(t);
*SELL(t)..   A_sf * Ecsp(t) =l= D(t) ;


model New /all/;

solve New min z using mip;
parameter Ecspt(t),LCOE,CSPcost,HRScost,LCOH,Tinst,Tmain,hrsinst,hrsmain, Epurchase,Esell,AnnualEnegy,AnnualCSP;
Ecspt(t)= A_sf.l*Ecsp(t)+  0.00000001;
cspcost =  A_Sf.l*IC_sf*CRF    +IC_PB*sum(t,Ecsp(t))*CRF   +A_sf.l*OM_sf    +OM_PB*sum(t,Ecsp(t))   ;
hrscost =  IC_E + IC_C + IC_HPT + IC_LPT + IC_CD +
           OM_E + OM_C + OM_HPT + OM_LPT + OM_CD ;
Tinst = z.l-(OM_E + OM_C + OM_HPT + OM_LPT + OM_CD+A_sf.l*OM_sf    +OM_PB*sum(t,Ecsp(t)));
Tmain = z.l-(IC_E + IC_C + IC_HPT + IC_LPT + IC_CD+A_Sf.l*IC_sf*CRF    +IC_PB*sum(t,Ecsp(t))*CRF   );
hrsinst =IC_E + IC_C + IC_HPT + IC_LPT + IC_CD;
hrsmain =OM_E + OM_C + OM_HPT + OM_LPT + OM_CD;
LCOE = Z.l/sum(t,Ecspt(t));
LCOH =  Z.l/ (4202 *360);
Epurchase=  sum(t, E_g.l(t)) ;
Esell= sum(t, E_exc.l(t)  )  ;

E_g.l(t)$(Not E_g.l(t)) = 0.0000001;
E_exc.l(t)$(Not E_exc.l(t)) = 0.0000001;

AnnualEnegy=  sum(t,Ecspt(t)) + Esell;
AnnualCSP=  sum(t,Ecspt(t));

Display A_SF.l, Z.l;
display AnnualCSP, Esell,Epurchase, AnnualEnegy, lcoe,lcoh,Tinst,Tmain,hrsinst,hrsmain,CSPcost,HRScost, Ecspt, E_g.l, E_exc.l; ;

parameter    EcspJan,  EcspFeb, EcspMar, EcspApr, EcspMay,  EcspJun, EcspJul, EcspAug, EcspSept, EcspOct, EcspNOV, EcspDec;
   EcspJan=  (A_sf.l* sum(t$(ord(t)>=1    and ord(t)<=744),  Ecsp(t)) );
    EcspFeb=  (A_sf.l*sum(t$(ord(t)>=745  and ord(t)<=1416), Ecsp(t)));
    EcspMar=  (A_sf.l*sum(t$(ord(t)>=1417 and ord(t)<=2160), Ecsp(t)));
    EcspApr=  (A_sf.l*sum(t$(ord(t)>=2161 and ord(t)<=2880), Ecsp(t)));
    EcspMay=  (A_sf.l*sum(t$(ord(t)>=2881 and ord(t)<=3624), Ecsp(t)));
    EcspJun=  (A_sf.l* sum(t$(ord(t)>=3625 and ord(t)<=4344), Ecsp(t)));
    EcspJul=  (A_sf.l*sum(t$(ord(t)>=4345 and ord(t)<=5088), Ecsp(t)));
    EcspAug=  (A_sf.l*sum(t$(ord(t)>=5089 and ord(t)<=5832), Ecsp(t)));
    EcspSept=  (A_sf.l*sum(t$(ord(t)>=5833 and ord(t)<=6552), Ecsp(t)));
    EcspOct=  (A_sf.l*sum(t$(ord(t)>=6553 and ord(t)<=7296), Ecsp(t)));
    EcspNOV=  (A_sf.l*sum(t$(ord(t)>=7297 and ord(t)<=8016), Ecsp(t)));
    EcspDec=  (A_sf.l*sum(t$(ord(t)>=8017 and ord(t)<=8760), Ecsp(t)));
Display    EcspJan,  EcspFeb, EcspMar, EcspApr, EcspMay,  EcspJun, EcspJul, EcspAug, EcspSept, EcspOct, EcspNOV, EcspDec;

parameter    E_gJan,  E_gFeb, E_gMar, E_gApr, E_gMay,  E_gJun, E_gJul, E_gAug, E_gSept, E_gOct, E_gNOV, E_gDec;
   E_gJan=  ( sum(t$(ord(t)>=1    and ord(t)<=744),  E_g.l(t)) );
    E_gFeb=  (sum(t$(ord(t)>=745  and ord(t)<=1416), E_g.l(t)));
    E_gMar=  (sum(t$(ord(t)>=1417 and ord(t)<=2160), E_g.l(t)));
    E_gApr=  (sum(t$(ord(t)>=2161 and ord(t)<=2880), E_g.l(t)));
    E_gMay=  (sum(t$(ord(t)>=2881 and ord(t)<=3624), E_g.l(t)));
    E_gJun=  ( sum(t$(ord(t)>=3625 and ord(t)<=4344), E_g.l(t)));
    E_gJul=  (sum(t$(ord(t)>=4345 and ord(t)<=5088), E_g.l(t)));
    E_gAug=  (sum(t$(ord(t)>=5089 and ord(t)<=5832), E_g.l(t)));
    E_gSept=  (sum(t$(ord(t)>=5833 and ord(t)<=6552), E_g.l(t)));
    E_gOct=  (sum(t$(ord(t)>=6553 and ord(t)<=7296), E_g.l(t)));
    E_gNOV=  (sum(t$(ord(t)>=7297 and ord(t)<=8016), E_g.l(t)));
    E_gDec=  (sum(t$(ord(t)>=8017 and ord(t)<=8760), E_g.l(t)));
Display    E_gJan,  E_gFeb, E_gMar, E_gApr, E_gMay,  E_gJun, E_gJul, E_gAug, E_gSept, E_gOct, E_gNOV, E_gDec;

parameter    E_excJan,  E_excFeb, E_excMar, E_excApr, E_excMay,  E_excJun, E_excJul, E_excAug, E_excSept, E_excOct, E_excNOV, E_excDec;
   E_excJan=  ( sum(t$(ord(t)>=1    and ord(t)<=744),  E_exc.l(t)) );
    E_excFeb=  (sum(t$(ord(t)>=745  and ord(t)<=1416), E_exc.l(t)));
    E_excMar=  (sum(t$(ord(t)>=1417 and ord(t)<=2160), E_exc.l(t)));
    E_excApr=  (sum(t$(ord(t)>=2161 and ord(t)<=2880), E_exc.l(t)));
    E_excMay=  (sum(t$(ord(t)>=2881 and ord(t)<=3624), E_exc.l(t)));
    E_excJun=  ( sum(t$(ord(t)>=3625 and ord(t)<=4344), E_exc.l(t)));
    E_excJul=  (sum(t$(ord(t)>=4345 and ord(t)<=5088), E_exc.l(t)));
    E_excAug=  (sum(t$(ord(t)>=5089 and ord(t)<=5832), E_exc.l(t)));
    E_excSept=  (sum(t$(ord(t)>=5833 and ord(t)<=6552), E_exc.l(t)));
    E_excOct=  (sum(t$(ord(t)>=6553 and ord(t)<=7296), E_exc.l(t)));
    E_excNOV=  (sum(t$(ord(t)>=7297 and ord(t)<=8016), E_exc.l(t)));
    E_excDec=  (sum(t$(ord(t)>=8017 and ord(t)<=8760), E_exc.l(t)));
Display    E_excJan,  E_excFeb, E_excMar, E_excApr, E_excMay,  E_excJun, E_excJul, E_excAug, E_excSept, E_excOct, E_excNOV, E_excDec;


parameter    DNIJan,  DNIFeb, DNIMar, DNIApr, DNIMay,  DNIJun, DNIJul, DNIAug, DNISept, DNIOct, DNINOV, DNIDec;
   DNIJan=  ( sum(t$(ord(t)>=1    and ord(t)<=744),  DNI(t)) );
    DNIFeb=  (sum(t$(ord(t)>=745  and ord(t)<=1416), DNI(t)));
    DNIMar=  (sum(t$(ord(t)>=1417 and ord(t)<=2160), DNI(t)));
    DNIApr=  (sum(t$(ord(t)>=2161 and ord(t)<=2880), DNI(t)));
    DNIMay=  (sum(t$(ord(t)>=2881 and ord(t)<=3624), DNI(t)));
    DNIJun=  ( sum(t$(ord(t)>=3625 and ord(t)<=4344), DNI(t)));
    DNIJul=  (sum(t$(ord(t)>=4345 and ord(t)<=5088), DNI(t)));
    DNIAug=  (sum(t$(ord(t)>=5089 and ord(t)<=5832), DNI(t)));
    DNISept=  (sum(t$(ord(t)>=5833 and ord(t)<=6552), DNI(t)));
    DNIOct=  (sum(t$(ord(t)>=6553 and ord(t)<=7296), DNI(t)));
    DNINOV=  (sum(t$(ord(t)>=7297 and ord(t)<=8016), DNI(t)));
    DNIDec=  (sum(t$(ord(t)>=8017 and ord(t)<=8760), DNI(t)));
Display    DNIJan,  DNIFeb, DNIMar, DNIApr, DNIMay,  DNIJun, DNIJul, DNIAug, DNISept, DNIOct, DNINOV, DNIDec;
