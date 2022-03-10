classdef CarMotor < handle
    %CarMotor is an object class used in the SAE Lap Sim in conjuction with
    %several other classes to build a Car class object.
    %   The car motor is currently defined by its output curve which is a
    %   function of RPM, its weight, center of gravity, mass moment of
    %   inertia, and the number of motors used in the car.
    
    properties
        OutputCurve  % Output curve of the car, given as a 3 column vector.
                     % The first column is RPM, the second is torque in in
                     % lb, and the third is efficiency.
        Weight       % Weight of the motor system given in lbs
        EffectiveCG  % 3 element vector giving motor cg in inches
        NMotors      % Number of motors used
        Name = '';
    end
    
    methods
        function M = CarMotor(OutputCurve,NMotors,Weight,CG)
            % CarMotor Constructor method
            %
            % This method constructs an object of type CarMotor.  To define
            % an object of this class, input an output curve, the number of
            % motors in the system, the weight of the system, the CG and
            % mass moment of inertia.
            %
            % INPUTS
            % Name          Type          Units   Description            
            %**************************************************************
            % OutputCurve   Nx3 matrix    mixed   Output curve of the
            %                                     motor, given in RPMs, in 
            %                                     lbs, and efficiency out 
            %                                     of 1.
            %
            % NMotors       int           N/A     Number of motors on car
            %
            % Weight        float         lbs     Weight of the motor(s)
            %
            % CG            1x3 array     in      CG of motor system
            %
            %
            % OUTPUTS
            % Name          Type          Units   Description            
            %**************************************************************
            % T             CarTire       N/A     CarTire Object
            %
            % VARIABLES
            % Name          Type          Units   Description            
            %**************************************************************
            % NONE
            %
            % FUNCTIONS
            % Name          Location         Description            
            %**************************************************************
            % NONE    
            
            % Assigns values to tire object properties
            M.OutputCurve = OutputCurve;
            M.NMotors = NMotors;
            M.Weight = Weight;
            M.EffectiveCG = CG;
        end
        
        function [ Torque, Efficiency ] = Output(M,RPM)
            I1 = find(M.OutputCurve(:,1) >= RPM, 1, 'first');
            I2 = find(M.OutputCurve(:,1) <= RPM, 1, 'last');
            Diff1 = M.OutputCurve(I1,1) - RPM;
            Diff2 = RPM - M.OutputCurve(I2,1);
            if Diff1 > Diff2
                Torque = M.OutputCurve(I2,2);
                Efficiency = M.OutputCurve(I2,3);
            else
                Torque = M.OutputCurve(I2,2);
                Efficiency = M.OutputCurve(I2,3);
            end
        end
        
        function Plot(M)
            [AX] = plotyy(M.OutputCurve(:,1),M.OutputCurve(:,2),...
                M.OutputCurve(:,1),M.OutputCurve(:,3));
            title([M.Name, ' Output Curve'])
            xlabel('Engine Speed (RPM)')
            set(get(AX(1),'Ylabel'),'String','Engine Torque (ft*lb)')
            set(get(AX(2),'Ylabel'),'String','Engine Efficiency (out of 1)')
            set(AX(2),'ylim',[0.0 1.0])
            grid on
        end
    end
    
end

