classdef MF52_MX
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FZ0 = 660
        R0 = 0
        LMX = 1        
        QSX1 = 0
        QSX2 = 0
        QSX3 = 0
    end
    
    methods
        function MXobj = MF52_MX(TIR)
            %Construct MX object
            MXobj.FZ0 = 660;              
            MXobj.R0 = TIR.UNLOADED_RADIUS;             
            MXobj.LMX = TIR.LMX;            
            MXobj.QSX1 = TIR.QSX1;
            MXobj.QSX2 = TIR.QSX2;            
            MXobj.QSX3 = TIR.QSX3;              
        end
        
        function MX = Calculate_Mx(MXobj,FZ,IA,FY)
            %Calculate MX
            NFY = FY./MXobj.FZ0;
            MX = FZ.*MXobj.R0.*(MXobj.QSX1 - MXobj.QSX2.*IA + MXobj.QSX3.*NFY).*MXobj.LMX;
        end
    end
end

