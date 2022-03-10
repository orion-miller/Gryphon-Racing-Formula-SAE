function FY = MF52_FY(SA,FZ,IA)
global FZ0 R0 KY
global PCY1 PDY1 PDY2 PDY3 ...
       PEY1 PEY2 PEY3 PEY4  ...
       PKY1 PKY2 PKY3 ...
       PHY1 PHY2 PHY3 ...
       PVY1 PVY2 PVY3 PVY4;

global LFZO LCX LMUX LEX LKX  LHX LVX LCY LMUY LEY LKY LHY LVY ...
       LGAY LTR LRES LGAZ LXAL LYKA LVYKA LS LSGKP LSGAL LGYR
         
SA  =  X(:,1)*pi/180;
FZ     =  abs(X(:,2));
IA  =  X(:,3).*pi/180;

GAMMAY = IA .* LGAY; %31 (%48 lgay=lg
FZ0PR  = FZ0   *  LFZO; %15,  NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ    = (FZ-FZ0PR ) ./ FZ0 ; %14,  (%30)

%c
%c -- lateral force (pure side slip)
%c
SHY     = (PHY1+PHY2 .* DFZ) .* LHY + PHY3 .* GAMMAY; %38,  (%55)
ALPHAY  = SA+SHY;  %30 (%47)
CY      = PCY1 .* LCY;  %32 (%49)
MUY     = (PDY1+PDY2 .* DFZ) .* (1.0-PDY3 .* GAMMAY.^2) .* LMUY; %34 (%51)
DY      = MUY .* FZ; %33 (%50)
KY      = PKY1 .* FZ0 .* sin(2.0 .* atan(FZ ./ (PKY2 .* FZ0 .* LFZO))) .* (1.0-PKY3 .* abs(GAMMAY)) .* LFZO .* LKY; %36 (%53)
BY      = KY ./ (CY .* DY);  %37 (%54)
% NOTE, PER SVEN @TNO: "SIGN(ALPHAY)"IS CORRECT AS IN DOCUMENTATION & BELOW; IT'S NOT SUPPOSED TO BE "SIGN(GAMMAY)"
EY      = (PEY1+PEY2 .* DFZ) .* (1.0-(PEY3+PEY4 .* GAMMAY) .* sign(ALPHAY)) .* LEY; %35 (%52)
% NOTE: LVY MULTIPLIES ONLY PVY1&2 IN DOCUMENTATION; ORIG VERSION MULT ALL TERMS
SVY     = FZ .* ((PVY1+PVY2 .* DFZ) .* LVY+(PVY3+PVY4 .* DFZ) .* GAMMAY) .* LMUY; %39 (%56)
FY0     = DY .* sin(CY .* atan(BY .* ALPHAY-EY .* (BY .* ALPHAY-atan(BY .* ALPHAY))))+SVY; %29 (%46)
FY      = FY0; %28
