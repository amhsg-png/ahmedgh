$Title Sets

Sets
       t        Months
       n        Years
       mt  months
*      objf     objective functions
       objf     objective functions      / TotalCost, LPSP,  CO2/ ;
*      objf     objective functions      / TotalCost, CO2 / ;


********************************************************************************************************************

$onecho > Sets.txt
set=t rng=t!a1 rdim=1
set=n rng=n!a1 rdim=1
set=mt rng=mt!a1 rdim=1

$offecho

$call gdxxrw.exe i=Sets.xlsx o=Sets.gdx @Sets.txt
$gdxin Sets
$load t n mt
$gdxin

********************************************************************************************************************

* the direction of the objective functions (1 for maximization and -1 for minimization)

$set min -1
$set max +1

Parameter dir(objf) direction of the objective functions / TotalCost    %min%, LPSP %min%,   CO2 %min%/;
*Parameter dir(objf) direction of the objective functions / TotalCost    %min%, CO2 %min%/;
*********************************************************************************************************************
*********************************************************************************************************************
