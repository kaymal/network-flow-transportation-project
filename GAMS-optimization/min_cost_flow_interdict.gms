$TITLE UsOilImport OA4202 2012
$INLINECOM { }
 OPTIONS
   SOLPRINT =     OFF,
   DECIMALS =       1,
   LIMCOL   =      10,
   LIMROW   =      10,
   RESLIM   =      60, {max seconds}
   ITERLIM  =99999999, {max iterations}
   LP       =     CPLEX,
   MIP      =     CPLEX
 ;


$offlisting

$ontext
============================================================================
With the following few lines of code, we are reading in the network specification.
- We read in the set of nodes, from the nodes.csv file
- We read in the set of all arcs, from the arcs_set.csv file
- We define the subset of arcs which we can interdict
- Finally, we read in the arc data, from the arcs_data.csv file
============================================================================
$offtext


SET
    n   nodes  /
$INCLUDE  nodes.csv
 /;

alias(n,i,j)   {(i,j) is for arc distance and flow}

SET
   arcs(i,j)         arclist  /
$ONDELIM
$INCLUDE  arcs.csv
$OFFDELIM
 /;

SET
   attackarcs(i,j)         the arcs on which we can place interdictions  /
$ONDELIM
$INCLUDE  attack_arcs_set.csv
$OFFDELIM
 /;

table arcdata(i,j,*)
$ONDELIM
$INCLUDE  arcs_data2.csv
$OFFDELIM
;

$ontext
============================================================================
We'll need the nC calculation to eliminate arcs completely
============================================================================
$offtext


SCALARS nC, unsatisfied_flow ;
 nC=0.0
;

 LOOP(arcs(i,j),
   IF( nC<ABS(arcdata(i,j,'Cost')),
       nC=ABS(arcdata(i,j,'Cost')) ;
   );
 );

nC = nC * card(n);

display nC;


$ontext
============================================================================
Lets define the primal min cost flow LP.

Our definition includes the parameters xbar(i,j).  We'll use them later, to
re-solve for the min cost flow, given an interdiction plan X.  For now, we'll just
set the xbar(i,j) to zero.
============================================================================
$offtext

 PARAMETERS
   delay
   xbar(i,j)
 ;

 delay = nC;

 LOOP(arcs(i,j),
   xbar(i,j)=0 ;
 );

 VARIABLE
   Zprimal           objective function value
 ;

 POSITIVE VARIABLES
   Y(i,j)            flow on arc i j
 ;

 EQUATIONS
   OBJECTIVE
   NETFLOW(n)
   CAPACITY(i,j)
   LOWERBOUND(i,j)
 ;

 OBJECTIVE..
   Zprimal =e=  SUM(arcs(i,j),(arcdata(i,j,'Cost') + delay*xbar(i,j)$(attackarcs(i,j) ) ) *Y(i,j) )
;

 NETFLOW(j)..
   SUM(arcs(i,j),Y(i,j)) - SUM(arcs(j,i),Y(j,i))      {+ flow in - flow out}
     =e= 0
 ;

 CAPACITY(arcs(i,j))..
   Y(i,j) =l= arcdata(i,j,'Upper')
 ;

 LOWERBOUND(arcs(i,j))..
   Y(i,j) =g= arcdata(i,j,'Lower')
 ;


 MODEL MinCostFlow
 /
   OBJECTIVE
   NETFLOW
   CAPACITY
   LOWERBOUND
 /;

$ontext
============================================================================
Now, lets solve for the min cost flow when there are no interdictions.
============================================================================
$offtext

 SOLVE MinCostFlow USING LP MINIMIZING Zprimal;       { find the uninterdicted min cost flow}

$ontext
============================================================================
Lets print the uninterdicted min cost flow solution to a file called modelOutput.out
============================================================================
$offtext


 file out /modelOutput.out/;
 put out;

 put 'nC is ', nC / ;
 put '****' / ;
 put 'solve MinCostFlow with no interdictions' / ;
 LOOP(arcs(i,j)$(Y.l(i,j)>0),
   put '   flow on arc ',i.tl,' -> ', j.tl, ' is ', Y.l(i,j) / ;
 );
 put 'total cost =',Zprimal.l:10:1 /
;


$ontext
============================================================================
Now, lets define the MIP for computing interdictions.

We'll introduce a scalar, 'nAttacks', that specifies the number of
interdictions to compute.

We want to maximize the min cost flow.  So, the problem we want to solve is a
'max min' problem.  To put this into the solver, we'll take the dual of the
inner problem, so the problem becomes a 'max max'.  That is why in the
following code, you'll see the dual of the min cost flow problem.

The binary variables X(i,j), will specify whether we attack arc i j.
============================================================================
$offtext


SCALARS nAttacks;

 VARIABLES
    Zdual
 ;

 VARIABLES
   Rho(n)              dual variables for the netflow constraints
 ;

 NEGATIVE VARIABLES
   Pi1(i,j)              dual variables for the capacity constraints
;
 POSITIVE VARIABLES
   Pi2(i,j)
;

 BINARY VARIABLES
   X(i,j)            should we attack arc i j
 ;

 EQUATIONS
   DUAL_OBJECTIVE
   DUAL_CONSTRAINTS(i,j)
   ATTACK_LIMIT
 ;

 DUAL_OBJECTIVE..
   Zdual =e= SUM(arcs(i,j), arcdata(i,j,'Upper')*Pi1(i,j)) +  SUM(arcs(i,j), arcdata(i,j,'Lower')*Pi2(i,j))
 ;

 DUAL_CONSTRAINTS(arcs(i,j))..
   Rho(j)-Rho(i)+Pi1(i,j)+ Pi2(i,j)=l= arcdata(i,j,'Cost') + delay*X(i,j)$(attackarcs(i,j) )
 ;

 ATTACK_LIMIT..
   SUM(attackarcs(i,j),X(i,j)) =l= nAttacks
 ;


 MODEL InterdictionMIP
 /
   DUAL_OBJECTIVE
   DUAL_CONSTRAINTS
   ATTACK_LIMIT
 /;

$ontext
============================================================================
Now, lets solve for the best interdiction plan.
============================================================================
$offtext



 InterdictionMIP.optfile = 1;
SET attack /1*4/;

LOOP (attack,
 nAttacks = ord(attack);
 SOLVE InterdictionMIP USING MIP MAXIMIZING Zdual;       { Find the best interdiction plan }

$ontext
============================================================================
Lets print the interdiction plan into the output file we had opened, modelOutput.out
============================================================================
$offtext

display X.l;

 put '****' / ;
 put 'interdiction plan with ', nAttacks, ' nAttacks:' / ;
 LOOP(arcs(i,j)$(X.l(i,j)>0),
   put '   attack arc: ',i.tl,' -> ', j.tl / ;
 );
 put 'cost with interdictions in place =',Zdual.l:10:1 /
;

$ontext
============================================================================
Finally, now that we have computed attack locations, lets re-solve for
the operator's new best flow.
============================================================================
$offtext

 LOOP(arcs(i,j),                        { transfer the interdiction solution to constant parameter }
   xbar(i,j)=X.l(i,j) ;
 );

 SOLVE MinCostFlow USING LP MINIMIZING Zprimal;       { find min cost flow plan with the interdictions in place }


  unsatisfied_flow = sum((i,j), Y.l(i,j)$(xbar(i,j) eq 1))
 ;


 put '****' / ;
 put 'solve MinCostFlow with interdictions in place' / ;
 put 'operator best response:' / ;
 LOOP(arcs(i,j)$(Y.l(i,j)>0),
   put '   flow on arc ',i.tl,' -> ', j.tl, ' is ', Y.l(i,j)/ ;
 );
 put 'total cost =',Zprimal.l:10:1 /
 put 'Unsatisfied Flow :   ;', unsatisfied_flow /
;
);
