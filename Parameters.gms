$Title Parameters
*$include Sets

Parameter  D(t)         Monthly demand
Parameter  DNI(t)         Monthly Radiation
Parameter  Tc(t)        Monthly Radiation
Parameter  Efpv(t)      PV effeciency
Parameter  Vrr(t)       Wind Speed
Parameters OMPpv(n)     Operation and maintenance price of PV system ($)
Parameters OMPwt(n)     Operation and maintenance price of wind turbines $


************************************************************************************************************
************************************************************************************************************

$onecho > Parameters.txt
par=D rng=D!a1 rdim=1
par=DNI rng=DNI!a1 rdim=1
par=Tc rng=Tc!a1 rdim=1
par=Vrr rng=Vrr!a1 rdim=1
par=Efpv rng=Efpv!a1 rdim=1
par=OMPpv rng=OMPpv!a1 rdim=1
par=OMPwt rng=OMPwt!a1 rdim=1

$offecho

$call gdxxrw.exe i=Parameters.xlsx o=Parameters.gdx @Parameters.txt
$gdxin Parameters
$load  D DNI Tc Efpv OMPpv  OMPwt Vrr
$gdxin
