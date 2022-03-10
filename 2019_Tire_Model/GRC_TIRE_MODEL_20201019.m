%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orion Miller            GRYPHON RACING TIRE MODEL          September 2019
%This program creates standard plots for cornering and drive brake test
%files, and also solves MF 5.2 fits of Fx, Fy and Mz.

% KNOWN ISSUES:
% For some files, the program misidentifies it as a longitudinal test
% rather than lateral, this can be "fixed" by defining the value of isCornering manually.
% This was code written to work with 12,10,14 psi run files from the later
% test rounds. Other files may not work as expected or at all.
% The radii show as around zero on the data vs time plots in some files.

% REFERENCES:
% Tire and Vehicle Dynamics 2nd edition, Hans Pacejka.
% Race Car Vehicle Dynamics, Milliken and Milliken.
% The Science of Vehicle Dynamics, Massimo Guiggiani.
% Michael Barnard (Pacejka 2012 Code)
% Bjorn Gustavsson (Scrolling Subplot Code)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STARTUP TASKS
clc 
clearvars
tic

Scheme = ColorSchemeDef(); %Define color scheme used for plotting

%% FILE IMPORT
try
[RunFile, RunPath] = uigetfile('DialogTitle','Select Run File to Open','*.mat*'); %Point to file
cd(RunPath); %Change to that directory
catch
    return
end

load(RunFile); %Load the mat file

%Assign chnls into data struct
RunData.Chnl.ET = ET; %Time
RunData.Chnl.V = V; %Roadway velocity
RunData.Chnl.N = N; %RPM of wheel
RunData.Chnl.SA = SA; %Slip angle
RunData.Chnl.IA = IA; %Inclination angle (camber)
RunData.Chnl.RL = RL; %Loaded radius
RunData.Chnl.RE = RE; %Effective radius
RunData.Chnl.P = P; %Pressure
RunData.Chnl.FX = FX; %Longitudinal force
RunData.Chnl.FY = FY; %Lateral force
RunData.Chnl.FZ = FZ; %Normal force
RunData.Chnl.MX = MX; %Overturning moment
RunData.Chnl.MZ = MZ; %Aligning moment
RunData.Chnl.NFX = NFX; %MuX (FX/FZ)
RunData.Chnl.NFY = NFY; %MuY (FY/FZ)
RunData.Chnl.RST = RST; %Road surface temp
RunData.Chnl.TSTI = TSTI; %Tread temp - inner
RunData.Chnl.TSTC = TSTC; %Tread temp - center
RunData.Chnl.TSTO = TSTO; %Tread temp - outer
RunData.Chnl.AMBTMP = AMBTMP; %Ambient temp
RunData.Chnl.SR = SR; %Slip ratio - effective radius definition
RunData.Chnl.SL = SL; %Slip ratio - loaded radius definition
%Clear the other vars out - keep things clean
clear ET V N SA IA RL RE P FX FY FZ MX MZ NFX NFY RST TSTI TSTC TSTO AMBTMP SR SL RUN
%Set some run info
RunData.Info.FileName = extractBefore(RunFile,'.mat'); %Name of file loaded
RunData.Info.Vars = channel.name; %Store chnl names and units
RunData.Info.Units = channel.units; %Store chnl names and units
RunData.Info.TireID = tireid; %Run test type - i.e. cornering, drive/brake
RunData.Info.TestID = testid; %Run test type - i.e. cornering, drive/brake
clear channel tireid testid

%% GENERATE TEST CONDITIONS VS TIME FIGURE

%Plot Control Variables, Temps, Radii
FIG = figure('Name',sprintf('Test Conditions vs. Time: %s, %s', RunData.Info.TireID, RunData.Info.TestID),'NumberTitle','off','Position',[0 0 1300 670],'Color',[1 1 1]);
set(FIG,'units','normalized','outerposition',[0,0,1,1]);
figure(FIG)
d= waitbar(0,'Processing data');

AX = gobjects(10,1); %Initialize graphics objects

AX(1) = scrollsubplot(3,1,1); %Slips vs time
    if strcmp(RunData.Info.TestID,'Cornering')
        plot(RunData.Chnl.ET,RunData.Chnl.SA,'LineWidth',1)
        title('Slip Angle vs. Time');
        ylabel('SA (Degrees)');
        grid on
    else %Add slip ratio to plot if run is drive brake
        title('Slips vs. Time');        
        xlabel('Time (s)');
        
        yyaxis left
        plot(RunData.Chnl.ET,RunData.Chnl.SA,'LineWidth',1)
        ylabel('Slip Angle (Degrees)','Color','k');
        hold on
        ylim([-15 15])
        
        yyaxis right
        plot(RunData.Chnl.ET,RunData.Chnl.SR,'LineWidth',1)
        ylabel('Slip Ratio','Color','k');
        hold off
        ylim([-0.25 0.25])
        
        legend('SA','SR');
        grid on
    end

AX(2) = scrollsubplot(3,1,2); %Pressure vs time
    plot(RunData.Chnl.ET,RunData.Chnl.P*0.145,'LineWidth',1)
    title('Pressure vs. Time');
    ylabel('P (PSI)');
    ylim([6 16])    
    grid on

AX(3) = scrollsubplot(3,1,3); %IA/Camber vs time
    plot(RunData.Chnl.ET,RunData.Chnl.IA,'LineWidth',1)
    title('Inclination Angle vs. Time');
    ylabel('IA (Degrees)');
    ylim([0 6])    
    grid on
    
AX(4) = scrollsubplot(3,1,4); %Normal load vs time
    plot(RunData.Chnl.ET,RunData.Chnl.FZ,'LineWidth',1)
    title('Normal Force vs. Time');
    ylabel('FZ (N)');
    yticks([-1540 -1100 -880 -660 -440 -220])
    grid on
    
AX(5) = scrollsubplot(3,1,5); %Friction forces vs time
    plot(RunData.Chnl.ET,RunData.Chnl.FX,'LineWidth',1)
    hold on
    plot(RunData.Chnl.ET,RunData.Chnl.FY,'LineWidth',1)
    hold off
    title('Friction Forces vs. Time');
    ylabel('Force (N)');
    legend('FX','FY')
    grid on
    
AX(6) = scrollsubplot(3,1,6); %Mu's vs time
    plot(RunData.Chnl.ET,RunData.Chnl.NFX,'Linewidth',1)
    hold on
    plot(RunData.Chnl.ET,RunData.Chnl.NFY,'Linewidth',1)
    hold off
    title('Normalized Friction Forces vs. Time');
    ylabel('Mu (-)');
    ylim([-4 4])
    legend('NFX','NFY')
    grid on
    
AX(7) = scrollsubplot(3,1,7); %Moments vs time
    plot(RunData.Chnl.ET,RunData.Chnl.MX,'LineWidth',1)
    hold on
    plot(RunData.Chnl.ET,RunData.Chnl.MZ,'LineWidth',1)
    hold off
    title('Moments vs. Time');
    ylabel('Torque (N.m)');
    legend('MX','MZ')
    grid on
    
AX(8) = scrollsubplot(3,1,8); %Speed vs time
    plot(RunData.Chnl.ET,RunData.Chnl.V,'LineWidth',1)
    title('Speed vs. Time');
    ylabel('V (kph)');
    ylim([0 100])
    grid on
    
AX(9) = scrollsubplot(3,1,9); %Radius vs time
    plot(RunData.Chnl.ET,RunData.Chnl.RE,'LineWidth',1)
    hold on
    plot(RunData.Chnl.ET,RunData.Chnl.RL,'LineWidth',1)
    hold off
    title('Radii vs. Time');
    ylabel('Radius (cm)');
    legend('RE', 'RL')
    grid on
    
AX(10) = scrollsubplot(3,1,10); %Temps vs time
    plot(RunData.Chnl.ET,RunData.Chnl.TSTI,'LineWidth',1)
    hold on
    plot(RunData.Chnl.ET,RunData.Chnl.TSTC,'LineWidth',1)
    plot(RunData.Chnl.ET,RunData.Chnl.TSTO,'LineWidth',1)
    plot(RunData.Chnl.ET,RunData.Chnl.AMBTMP,'LineWidth',1)
    plot(RunData.Chnl.ET,RunData.Chnl.RST,'LineWidth',1)
    hold off
    title('Temperatures vs. Time');
    xlabel('Time (s)');
    ylabel('Temperature (Degrees C)');
    legend('TSTI', 'TSTC', 'TSTO','RST','AMBTMP')
    grid on
    
    linkaxes(AX,'x');
    
%% DATA PROCESSING
    
% Splits curves with Run Data
dif = diff(RunData.Chnl.ET);
SWEEP_ENDS = find(dif>1); % creates finish splits
SWEEP_STARTS = [1; SWEEP_ENDS + 1]; % creates beginning splits
SWEEP_ENDS(numel(SWEEP_ENDS)+1) = numel(RunData.Chnl.ET); % adds final finish split
SWEEPS_NUM = numel(SWEEP_STARTS);

% Length = numel(RunData.Chnl.ET);

FZ_LIST = [0 -220 -440 -660 -880 -1100 -1540 -2000 -3000];  %Middle 7 values are the vertical loads tested

%Define custom colormap to be used with temp heatplots
CMTemp1 = parula(64);
CMTemp1 = CMTemp1(1:end-12,:);
CMTemp2 = jet(64);
CMTemp2 = CMTemp2(45:64,:); 
CustomColorMap = [CMTemp1; CMTemp2];
CustomColorMap(52,3) = 0.1;
clear CMTemp1 CMTemp2

for i=1:SWEEPS_NUM-1   %Finds tire spring rate by calculating differences in avg FZ and RL from sweep to sweep
    
stiffnessFZ1 = mean(RunData.Chnl.FZ(SWEEP_STARTS(i):SWEEP_ENDS(i)));
stiffnessFZ2 = mean(RunData.Chnl.FZ(SWEEP_STARTS(i+1):SWEEP_ENDS(i+1)));

stiffnessRL1 = mean(RunData.Chnl.RL(SWEEP_STARTS(i):SWEEP_ENDS(i)));
stiffnessRL2 = mean(RunData.Chnl.RL(SWEEP_STARTS(i+1):SWEEP_ENDS(i+1)));

tireRate(i) = abs(stiffnessFZ1 - stiffnessFZ2)./abs(stiffnessRL1 - stiffnessRL2);
end

tireRate = mean(tireRate);
tireRate = tireRate./10;
stiffnessStatement = sprintf('Tire Spring Rate: %.2f N/mm', tireRate);
disp(stiffnessStatement); %Show stiffness statement in cmd line

RunData.Chnl.CmdP = [];
RunData.Chnl.CmdIA = [];
RunData.Chnl.CmdFZ = [];
RunData.Chnl.CmdSA = [];
RunData.Chnl.CmdV = [];

for i=1:SWEEPS_NUM %Loop thru sweeps to find properties
    SWEEP_P = RunData.Chnl.P(SWEEP_STARTS(i):SWEEP_ENDS(i));
    SWEEP_IA = RunData.Chnl.IA(SWEEP_STARTS(i):SWEEP_ENDS(i)); 
    SWEEP_FZ = RunData.Chnl.FZ(SWEEP_STARTS(i):SWEEP_ENDS(i));    
    SWEEP_SA = RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i));
    SWEEP_V = RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i));        
  
    P_AVG(i) = mean(SWEEP_P);
    IA_AVG(i) = mean(SWEEP_IA);
    FZ_AVG(i) = mean(SWEEP_FZ);    
    SA_AVG(i) = mean(SWEEP_SA);
    V_AVG(i) = mean(SWEEP_SA); 
    
    P_CLASS(i) = round(P_AVG(i)*0.145);
    IA_CLASS(i) = round(IA_AVG(i));
    FZ_CLASS(i) = interp1(FZ_LIST,FZ_LIST,FZ_AVG(i),'nearest');    
    SA_CLASS(i) = round(SA_AVG(i));
    V_CLASS(i) = round(V_AVG(i));    
    
    RunData.Chnl.CmdP(SWEEP_STARTS(i):SWEEP_ENDS(i)) = P_CLASS(i);
    RunData.Chnl.CmdIA(SWEEP_STARTS(i):SWEEP_ENDS(i)) = IA_CLASS(i);
    RunData.Chnl.CmdFZ(SWEEP_STARTS(i):SWEEP_ENDS(i)) = FZ_CLASS(i);    
    RunData.Chnl.CmdSA(SWEEP_STARTS(i):SWEEP_ENDS(i)) = SA_CLASS(i);
    RunData.Chnl.CmdV(SWEEP_STARTS(i):SWEEP_ENDS(i)) = V_CLASS(i);    
    
    %Extra stats - not used currently
%     FZ_STDDEV(i) = std(SWEEP_FZ);
%     P_STDDEV(i) = std(SWEEP_P);
%     IA_STDDEV(i) = std(SWEEP_IA);

%     SWEEP_FX = FX(SWEEP_STARTS(i):SWEEP_ENDS(i));
%     SWEEP_FY = FY(SWEEP_STARTS(i):SWEEP_ENDS(i));
%     SWEEP_MZ = MZ(SWEEP_STARTS(i):SWEEP_ENDS(i));%     FY_PEAK_POS(i) = max(SWEEP_FY); %Finds max grips for each sweep

%     FX_PEAK_POS(i) = max(SWEEP_FX);
%     MZ_PEAK_POS(i) = max(SWEEP_MZ);    
%     FY_PEAK_NEG(i) = min(SWEEP_FY);
%     FX_PEAK_NEG(i) = min(SWEEP_FX);
%     MZ_PEAK_NEG(i) = min(SWEEP_MZ);
        
%     TSTI_AVG(i) = mean(TSTI(SWEEP_STARTS(i):SWEEP_ENDS(i)));
%     TSTC_AVG(i) = mean(TSTC(SWEEP_STARTS(i):SWEEP_ENDS(i)));
%     TSTO_AVG(i) = mean(TSTO(SWEEP_STARTS(i):SWEEP_ENDS(i)));
end

%% FITTING CALLS
%Check if user loading tire model
TIRFlag = questdlg('Load TIR File?', ...
	'TIR Load', ...
	'No','Yes','Yes');

if strcmp(TIRFlag,'Yes')
    global FitData
    FitData = TIR_Reader();
    FitData.Output = []; %Add field to store calculated model output    
    %Check if user running fitting or just visualizing model as is
    FittingFlag = questdlg('Run model fitting process?', ...
        'Model Fitting', ...
        'No','Yes','No');
    
    if strcmp(FittingFlag,'Yes')
        FittingFlag = 1;
    else
        FittingFlag = 0;
    end
else
    FittingFlag = 0;
end

if FittingFlag
    if strcmp(RunData.Info.TestID,'Cornering') %If cornering file
        waitbar(0.1,d,'FY Fitting');       
%         [FitData] = ModelOptimizer(RunData,FitData,'FY PURE',d);
        [FitData] = ModelOptimizer(RunData,FitData,'MX PURE',d);        
%         [FitData] = ModelOptimizer(RunData,FitData,'MZ PURE',d);         
    else %If drive brake file
        waitbar(0.1,d,'FX Fitting');    
        [FitData] = ModelOptimizer(RunData,FitData,'FX PURE',d);
        [FitData] = ModelOptimizer(RunData,FitData,'FX COMB',d);        
        [FitData] = ModelOptimizer(RunData,FitData,'FY COMB',d); 
%         [FitData] = ModelOptimizer(RunData,FitData,'MZ COMB',d);         
    end
end

close(d) %Close progress bar
Time = toc; %Processing time
fprintf(sprintf('Processing Time: %.1f seconds\n', Time))
clear Time

%% PLOTTING LOOP
while(1) %"Dummy loop" This is done to make a loop that automatically continues until broken by user
    disp('Enter a value:');
    
    if strcmp(RunData.Info.TestID,'Cornering') %If cornering file
        disp('1 - Normal Force Variation, 2 - Pressure Variation, 3 - Camber Variation');
        disp('5 - Close Figures, 6 - Stop');
    else %If drive brake file
        disp('1 - Normal Force Variation, 2 - Pressure Variation, 3 - Camber Variation, 4 - Slip Angle Variation');
        disp('5 - Close Figures, 6 - Stop');
    end    
    VarChoice = input(''); %Stores the user response
    
    %% PARSING SETUP
    %There is functionality to vary 4 different
    %parameters, holding the others constant. "Fixed 1" and "Fixed 2" are
    %used in each situation to define the things being held constant. "Class 1"
    %and "Class 2" bring in the classifications of those variables for each
    %sweep so that the plotting code can then iterate through each loop and
    %check whether these variables match with the specified conditions. "Class
    %3" is used to seperate out the different values for the parameter being
    %varied so that they can be plotted by colour. For drive brake there is an
    %extra parameter because of the different SA conditions, so there are 3
    %fixed parameters and one varied (CLASS 4).
    Legend = {};
    %-----------------------------------------------------------------------------------------------------%
    % CORNERING PARSING SETUP
    if strcmp(RunData.Info.TestID,'Cornering') %If cornering file
        %-----------------------------------------------------------------------------------------------------%
        switch VarChoice %Run parsing setup according to parameter being varied
            %-----------------------------------------------------------------------------------------------------%
            case 1 %FZ variation
                P_Choice = input('Choose Pressure: ');
                IA_Choice = input('Choose Inclination Angle: ');
                Fixed_1 = P_Choice;
                Fixed_2 = IA_Choice;
                Class_1 = P_CLASS;
                Class_2 = IA_CLASS;
                
                Class_3 = FZ_CLASS;
                Class_3(Class_3 == -220) = 1;
                Class_3(Class_3 == -440) = 2;
                Class_3(Class_3 == -660) = 3;
                Class_3(Class_3 == -880) = 4;
                Class_3(Class_3 == -1100) = 5;
                Class_3(Class_3 == -1540) = 6;
                Class_3(Class_3 == -2000) = 7;
                
                Legend{1} = 'FZ 220'; %Set legend entries
                Legend{2} = 'FZ 440';
                Legend{3} = 'FZ 660';
                Legend{4} = 'FZ 880';
                Legend{5} = 'FZ 1100';
                Legend{6} = 'FZ 1540';
                Legend{7} = 'FZ 2000';
                
                VariedParameter = ('Normal Force Variation'); %This info in this var. is used for the label in the top of each window
            %-----------------------------------------------------------------------------------------------------%
            case 2 %Pressure variation
                FZ_CHOICE = input('Choose Vertical Force: ');
                IA_Choice = input('Choose Inclination Angle: ');
                Fixed_1 = -FZ_CHOICE;
                Fixed_2 = IA_Choice;
                Class_1 = FZ_CLASS;
                Class_2 = IA_CLASS;
                
                Class_3 = P_CLASS;
                Class_3(Class_3 == 8) = 1;
                Class_3(Class_3 == 10) = 2;
                Class_3(Class_3 == 12) = 3;
                Class_3(Class_3 == 14) = 4;
                
                Legend{1} = 'P 8'; %Set legend entries
                Legend{2} = 'P 10';
                Legend{3} = 'P 12';
                Legend{4} = 'P 14';
                
                VariedParameter = ('Pressure Variation');
            %-----------------------------------------------------------------------------------------------------%
            case 3 %Camber variation
                FZ_CHOICE = input('Choose Vertical Force: ');
                P_Choice = input('Choose Pressure: ');
                Fixed_1 = -FZ_CHOICE;
                Fixed_2 = P_Choice;
                Class_1 = FZ_CLASS;
                Class_2 = P_CLASS;
                
                Class_3 = IA_CLASS + 1;
                
                Legend{1} = 'IA 0'; %Set legend entries
                Legend{2} = 'IA 1';
                Legend{3} = 'IA 2';
                Legend{4} = 'IA 3';
                Legend{5} = 'IA 4';
                Legend{6} = 'IA 5';
                Legend{7} = 'IA 6';
                
                VariedParameter = ('Camber Variation');
            %-----------------------------------------------------------------------------------------------------%
            case 4
                warndlg('Slip angle selection only available for drive/brake files','error')
            %-----------------------------------------------------------------------------------------------------%
            case 5
                close all
                continue
            %-----------------------------------------------------------------------------------------------------%
            case 6
                break
            %-----------------------------------------------------------------------------------------------------%
        end
        %-----------------------------------------------------------------------------------------------------%
        % D/B PARSING SETUP
    else %Else if drive brake file
        %-----------------------------------------------------------------------------------------------------%
        switch VarChoice %Run parsing setup according to parameter being varied
            %-----------------------------------------------------------------------------------------------------%
            case 1 %FZ variation
                P_Choice = input('Choose Pressure: ');
                IA_Choice = input('Choose Inclination Angle: ');
                SA_Choice = input('Choose Slip Angle: ');
                Fixed_1 = P_Choice;
                Fixed_2 = IA_Choice;
                Fixed_3 = -SA_Choice;
                Class_1 = P_CLASS;
                Class_2 = IA_CLASS;
                Class_3 = SA_CLASS;
                
                Class_4 = FZ_CLASS;
                Class_4(Class_4 == -220) = 1;
                Class_4(Class_4 == -440) = 2;
                Class_4(Class_4 == -660) = 3;
                Class_4(Class_4 == -880) = 4;
                Class_4(Class_4 == -1100) = 5;
                Class_4(Class_4 == -1540) = 6;
                Class_4(Class_4 == -2000) = 7;
                
                Legend{1} = 'FZ 220'; %Set legend entries
                Legend{2} = 'FZ 440';
                Legend{3} = 'FZ 660';
                Legend{4} = 'FZ 880';
                Legend{5} = 'FZ 1100';
                Legend{6} = 'FZ 1540';
                Legend{7} = 'FZ 2000';
                
                VariedParameter = ('Normal Force Variation'); %This info in this var. is used for the label in the top of each window
            %-----------------------------------------------------------------------------------------------------%
            case 2 %Pressure variation
                FZ_CHOICE = input('Choose Vertical Force: ');
                IA_Choice = input('Choose Inclination Angle: ');
                SA_Choice = input('Choose Slip Angle: ');
                Fixed_1 = -FZ_CHOICE;
                Fixed_2 = IA_Choice;
                Fixed_3 = -SA_Choice;
                Class_1 = FZ_CLASS;
                Class_2 = IA_CLASS;
                Class_3 = SA_CLASS;
                
                Class_4 = P_CLASS;
                Class_4(Class_4 == 8) = 1;
                Class_4(Class_4 == 10) = 2;
                Class_4(Class_4 == 12) = 3;
                Class_4(Class_4 == 14) = 4;
                
                Legend{1} = 'P 8'; %Set legend entries
                Legend{2} = 'P 10';
                Legend{3} = 'P 12';
                Legend{4} = 'P 14';
                
                VariedParameter = ('Pressure Variation');
            %-----------------------------------------------------------------------------------------------------%
            case 3 %Camber variation
                FZ_CHOICE = input('Choose Vertical Force: ');
                P_Choice = input('Choose Pressure: ');
                SA_Choice = input('Choose Slip Angle: ');
                Fixed_1 = -FZ_CHOICE;
                Fixed_2 = P_Choice;
                Fixed_3 = -SA_Choice;
                Class_1 = FZ_CLASS;
                Class_2 = P_CLASS;
                Class_3 = SA_CLASS;
                
                Class_4 = IA_CLASS + 1;
                
                Legend{1} = 'IA 0'; %Set legend entries
                Legend{2} = 'IA 1';
                Legend{3} = 'IA 2';
                Legend{4} = 'IA 3';
                Legend{5} = 'IA 4';
                Legend{6} = 'IA 5';
                Legend{7} = 'IA 6';
                
                VariedParameter = ('Camber Variation');
            %-----------------------------------------------------------------------------------------------------%
            case 4 %Slip angle variation
                FZ_CHOICE = input('Choose Vertical Force: ');
                P_Choice = input('Choose Pressure: ');
                IA_Choice = input('Choose Inclination Angle: ');
                Fixed_1 = -FZ_CHOICE;
                Fixed_2 = P_Choice;
                Fixed_3 = IA_Choice;
                Class_1 = FZ_CLASS;
                Class_2 = P_CLASS;
                Class_3 = IA_CLASS;
                
                Class_4 = -SA_CLASS + 1;
                
                Legend{1} = 'SA 0'; %Set legend entries
                Legend{2} = 'SA 1';
                Legend{3} = 'SA 2';
                Legend{4} = 'SA 3';
                Legend{5} = 'SA 4';
                Legend{6} = 'SA 5';
                Legend{7} = 'SA 6';
                
                VariedParameter = ('Slip Angle Variation');
            %-----------------------------------------------------------------------------------------------------%
            case 5
                close all
                continue
            %-----------------------------------------------------------------------------------------------------%
            case 6
                break
            %-----------------------------------------------------------------------------------------------------%
        end
    end
    %-----------------------------------------------------------------------------------------------------%
    %% PLOTTING
    P = gobjects(7,1); %Initialize line objects to pull legend names from
    LegLabels = {}; %Stores legend labels
    %-----------------------------------------------------------------------------------------------------%
    % CORNERING PLOTTING
    if strcmp(RunData.Info.TestID,'Cornering')
        %-----------------------------------------------------------------------------------------------------%
        %FY VS SA PLOT
        FIG = figure('Name',sprintf('FY vs. Slip Angle | %s | %s', VariedParameter, RunData.Info.TireID), ...
            'Units','normalized','OuterPosition',[0,0,1,1],'NumberTitle','off','Color',[1 1 1]);
        AX = axes(FIG);
        hold(AX,'on')
        for i = 1:SWEEPS_NUM
            if Class_1(i)==Fixed_1&&Class_2(i)==Fixed_2
                P(Class_3(i)) = scatter(AX, RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i)), RunData.Chnl.FY(SWEEP_STARTS(i):SWEEP_ENDS(i)), 8, Scheme.LineColor{Class_3(i)}, 'DisplayName',[RunData.Info.FileName,' ',Legend{Class_3(i)}]);
                if FittingFlag
                    plot(AX, RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i)), FitData.Output.FY(SWEEP_STARTS(i):SWEEP_ENDS(i)), 'Color',Scheme.LineColor{Class_3(i)}, 'LineStyle','-', 'LineWidth',1.5);
                end
            end
        end
        %Axes props
        AX.Title.String = {'FY vs. SA', sprintf('%s', RunData.Info.TireID)};
        AX.YLabel.String = 'FY (N)';        
        AX.XLabel.String = 'SA (degrees)';
        AX.FontSize = 14;
        AX.XTick = -12.5:2.5:12.5;
        AX.YTick = -4000:1000:4000;
        grid on
        set(AX,'xminorgrid','on','yminorgrid','on')
        %Legend setup
        DeleteIdx = [];
        for iLines = 1:length(P) %Remove entries for unused conditions
            if ~strcmp(class(P(iLines)),'matlab.graphics.chart.primitive.Scatter')
                DeleteIdx = [DeleteIdx,iLines];
            end            
        end
        P(DeleteIdx) = [];
        legend(AX,P)
        %-----------------------------------------------------------------------------------------------------%
        %FY VS SA, TEMP PLOT
        FIG = figure('Name',sprintf('FY vs. Slip Angle, Tread Temperature | %s | %s', VariedParameter, RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
        set(FIG,'units','normalized','outerposition',[0,0,1,1]);
        plotinner = subplot(1,3,1);
        for i = 1:SWEEPS_NUM
            if Class_1(i)==Fixed_1&&Class_2(i)==Fixed_2
                scatter(plotinner,RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.FY(SWEEP_STARTS(i):SWEEP_ENDS(i)), ...
                    20,RunData.Chnl.TSTI(SWEEP_STARTS(i):SWEEP_ENDS(i)),'filled')
                colormap(CustomColorMap)
                colorbar
                hold on
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                title('FY vs. SA, TSTI')
                
                plotcenter = subplot(1,3,2);
                scatter(plotcenter,RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.FY(SWEEP_STARTS(i):SWEEP_ENDS(i)), ...
                    20,RunData.Chnl.TSTC(SWEEP_STARTS(i):SWEEP_ENDS(i)),'filled')
                colormap(CustomColorMap)
                colorbar
                hold on
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                title('FY vs. SA, TSTC')
                
                plotouter = subplot(1,3,3);
                scatter(plotouter,RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.FY(SWEEP_STARTS(i):SWEEP_ENDS(i)), ...
                    20,RunData.Chnl.TSTO(SWEEP_STARTS(i):SWEEP_ENDS(i)),'filled')
                colormap(CustomColorMap)
                colorbar
                hold on
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                title('FY vs. SA, TSTO')
            end
        end
        %-----------------------------------------------------------------------------------------------------%
        %MZ VS SA PLOT
        FIG = figure('Name',sprintf('MZ vs. Slip Angle | %s | %s', VariedParameter, RunData.Info.TireID), ...
            'Units','normalized','OuterPosition',[0,0,1,1],'NumberTitle','off','Color',[1 1 1]);
        AX = axes(FIG);
        hold(AX,'on')
        for i = 1:SWEEPS_NUM
            if Class_1(i)==Fixed_1&&Class_2(i)==Fixed_2
                P(Class_3(i)) = scatter(AX, RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i)), RunData.Chnl.MZ(SWEEP_STARTS(i):SWEEP_ENDS(i)), 8, Scheme.LineColor{Class_3(i)}, 'DisplayName',[RunData.Info.FileName,' ',Legend{Class_3(i)}]);
                if FittingFlag
                    plot(AX, RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i)), FitData.Output.MZ(SWEEP_STARTS(i):SWEEP_ENDS(i)), 'Color',Scheme.LineColor{Class_3(i)}, 'LineStyle','-', 'LineWidth',1.5);
                end
            end
        end
        %Axes props
        AX.Title.String = {'MZ vs. SA', sprintf('%s', RunData.Info.TireID)};
        AX.YLabel.String = 'MZ (N)';
        AX.XLabel.String = 'SA (degrees)';        
        AX.FontSize = 14;
        AX.XTick = -12.5:2.5:12.5;
        AX.YTick = -160:20:160;
        grid on
        set(AX,'xminorgrid','on','yminorgrid','on')
        %Legend setup
        DeleteIdx = [];
        for iLines = 1:length(P) %Remove entries for unused conditions
            if ~strcmp(class(P(iLines)),'matlab.graphics.chart.primitive.Scatter')
                DeleteIdx = [DeleteIdx,iLines];
            end            
        end
        P(DeleteIdx) = [];
        legend(AX,P)
        %-----------------------------------------------------------------------------------------------------%
        %MU PLOTS
        switch VariedParameter %Make mu plot depending on the parameter being varied as independent axis
            %-----------------------------------------------------------------------------------------------------%
            case 'Normal Force Variation' %Mu vs FZ plot
                FIG = figure('Name',sprintf('NFY vs. Normal Force | %s | %s', RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2
                        plot(abs(RunData.Chnl.FZ(SWEEP_STARTS(i):SWEEP_ENDS(i))),RunData.Chnl.NFY(SWEEP_STARTS(i):SWEEP_ENDS(i)),'.b');
                        hold on
                    end
                end
                title('NFY vs. FZ ','FontSize',20);
                xlabel('FZ (N)','FontSize',16);
                ylabel('NFY (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',220:220:1540);
                set(gca,'YTick',-5:1:5);
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
            %-----------------------------------------------------------------------------------------------------%                
            case 'Pressure Variation' %Mu vs P plot
                FIG = figure('Name',sprintf('NFY vs. Pressure | %s | %s', RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2
                        plot(0.145*RunData.Chnl.P(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.NFY(SWEEP_STARTS(i):SWEEP_ENDS(i)),'.r');
                        hold on
                    end
                end
                title('NFY vs. P ','FontSize',20);
                xlabel('P (psi)','FontSize',16);
                ylabel('NFY (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',8:2:14);
                set(gca,'YTick',-5:1:5);
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
            %-----------------------------------------------------------------------------------------------------%                
            case 'Camber Variation' %Mu vs IA plot
                FIG = figure('Name',sprintf('NFY vs. Camber | %s | %s', RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2
                        plot(RunData.Chnl.IA(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.NFY(SWEEP_STARTS(i):SWEEP_ENDS(i)),'.g');
                        hold on
                    end
                end
                title('NFY vs. IA ','FontSize',20);
                xlabel('IA (Degrees)','FontSize',16);
                ylabel('NFY (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',0:2:6);
                set(gca,'YTick',-5:1:5);
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                %-----------------------------------------------------------------------------------------------------%
        end
        %-----------------------------------------------------------------------------------------------------%
        % DRIVE/BRAKE PLOTTING
    else
        %-----------------------------------------------------------------------------------------------------%
        % FX vs SR PLOT
        FIG = figure('Name',sprintf('FX vs. Slip Ratio | %s | %s', VariedParameter, RunData.Info.TireID), ...
            'Units','normalized','OuterPosition',[0,0,1,1],'NumberTitle','off','Color',[1 1 1]);
        AX = axes(FIG);
        hold(AX,'on')
        for i = 1:SWEEPS_NUM
            if Class_1(i)==Fixed_1&&Class_2(i)==Fixed_2
                P(Class_3(i)) = scatter(AX, RunData.Chnl.SR(SWEEP_STARTS(i):SWEEP_ENDS(i)), RunData.Chnl.FX(SWEEP_STARTS(i):SWEEP_ENDS(i)), 8, Scheme.LineColor{Class_3(i)}, 'DisplayName',[RunData.Info.FileName,' ',Legend{Class_3(i)}]);
                if FittingFlag
                    plot(AX, RunData.Chnl.SR(SWEEP_STARTS(i):SWEEP_ENDS(i)), FitData.Output.FX(SWEEP_STARTS(i):SWEEP_ENDS(i)), 'Color',Scheme.LineColor{Class_3(i)}, 'LineStyle','-', 'LineWidth',1.5);
                end
            end
        end
        %Axes props
        AX.Title.String = {'FX vs. SR', sprintf('%s', RunData.Info.TireID)};
        AX.YLabel.String = 'FX (N)';
        AX.XLabel.String = 'SR (-)';        
        AX.FontSize = 14;
        AX.XTick = -1.0:0.04:1.0;
        AX.YTick = -5000:1000:5000;
        grid on
        set(AX,'xminorgrid','on','yminorgrid','on')
        %Legend setup
        DeleteIdx = [];
        for iLines = 1:length(P) %Remove entries for unused conditions
            if ~strcmp(class(P(iLines)),'matlab.graphics.chart.primitive.Scatter')
                DeleteIdx = [DeleteIdx,iLines];
            end            
        end
        P(DeleteIdx) = [];
        legend(AX,P)
        %-----------------------------------------------------------------------------------------------------%        
        %% FX vs SR, TEMP PLOT
        FIG = figure('Name',sprintf('FX vs. Slip Ratio, Tire Temperature | %s | %s', VariedParameter, RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
        set(FIG,'units','normalized','outerposition',[0,0,1,1]);
        plotinner = subplot(1,3,1);
        for i = 1:SWEEPS_NUM
            if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2 && Class_3(i)==Fixed_3
                scatter(plotinner,RunData.Chnl.SR(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.FX(SWEEP_STARTS(i):SWEEP_ENDS(i)), ...
                    20,RunData.Chnl.TSTI(SWEEP_STARTS(i):SWEEP_ENDS(i)),'filled')
                colormap(CustomColorMap)
                colorbar
                hold on
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                title('FX vs. SR, TSTI')
                
                plotcenter = subplot(1,3,2);
                scatter(plotcenter,RunData.Chnl.SR(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.FX(SWEEP_STARTS(i):SWEEP_ENDS(i)), ...
                    20,RunData.Chnl.TSTC(SWEEP_STARTS(i):SWEEP_ENDS(i)),'filled')
                colormap(CustomColorMap)
                colorbar
                hold on
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                title('FX vs. SR, TSTC')
                
                plotouter = subplot(1,3,3);
                scatter(plotouter,RunData.Chnl.SR(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.FX(SWEEP_STARTS(i):SWEEP_ENDS(i)), ...
                    20,RunData.Chnl.TSTO(SWEEP_STARTS(i):SWEEP_ENDS(i)),'filled')
                colormap(CustomColorMap)
                colorbar
                hold on
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                title('FX vs. SR, TSTO')
            end
        end
        %-----------------------------------------------------------------------------------------------------%        
        %% MU PLOTS
        switch VariedParameter %Make mu plot depending on the parameter being varied as independent axis
            %-----------------------------------------------------------------------------------------------------%            
            case 'Normal Force Variation' %Mu vs FZ plot
                FIG = figure('Name',sprintf('NFX vs. Normal Force | %s | %s', RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2 && Class_3(i)==Fixed_3
                        plot(abs(RunData.Chnl.FZ(SWEEP_STARTS(i):SWEEP_ENDS(i))),RunData.Chnl.NFX(SWEEP_STARTS(i):SWEEP_ENDS(i)),'.b');
                        hold on
                    end
                end
                title('NFX vs. FZ ','FontSize',20);
                xlabel('FZ (N)','FontSize',16);
                ylabel('NFX (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',220:220:1540);
                set(gca,'YTick',-5:1:5);
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
            %-----------------------------------------------------------------------------------------------------%                
            case 'Pressure Variation' %Mu vs P plot
                FIG = figure('Name',sprintf('NFX vs. Pressure | %s | %s', RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2 && Class_3(i)==Fixed_3
                        plot(0.145*RunData.Chnl.P(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.NFX(SWEEP_STARTS(i):SWEEP_ENDS(i)),'.r');
                        hold on
                    end
                end
                title('NFX vs. P ','FontSize',20);
                xlabel('P (psi)','FontSize',16);
                ylabel('NFX (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',8:2:14);
                set(gca,'YTick',-5:1:5);
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
            %-----------------------------------------------------------------------------------------------------%                
            case 'Camber Variation' %Mu vs IA plot
                FIG = figure('Name',sprintf('NFX vs. Camber | %s | %s', RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2 && Class_3(i)==Fixed_3
                        plot(RunData.Chnl.IA(SWEEP_STARTS(i):SWEEP_ENDS(i)),RunData.Chnl.NFX(SWEEP_STARTS(i):SWEEP_ENDS(i)),'.g');
                        hold on
                    end
                end
                title('NFX vs. IA ','FontSize',20);
                xlabel('IA (Degrees)','FontSize',16);
                ylabel('NFX (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',0:2:6);
                set(gca,'YTick',-5:1:5);
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
            %-----------------------------------------------------------------------------------------------------%    
            case 'Slip Angle Variation' %Mu vs SA plot
                FIG = figure('Name',sprintf('NFX vs. Slip Angle | %s | %s', RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1 && Class_2(i)==Fixed_2 && Class_3(i)==Fixed_3
                        plot(abs(RunData.Chnl.SA(SWEEP_STARTS(i):SWEEP_ENDS(i))),RunData.Chnl.NFX(SWEEP_STARTS(i):SWEEP_ENDS(i)),'.m');
                        hold on
                    end
                end
                title('NFX vs. SA ','FontSize',20);
                xlabel('\alpha (Degrees)','FontSize',16);
                ylabel('NFX (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',0:1:6);
                set(gca,'YTick',-5:1:5);
                grid on
                set(gca,'xminorgrid','on','yminorgrid','on')
                
                %Traction circle plot
                FIG = figure('Name',sprintf('Traction Circle | %s | %s', VariedParameter, RunData.Info.TireID),'NumberTitle','off','Color',[1 1 1]);
                set(FIG,'units','normalized','outerposition',[0,0,1,1]);
                AX = axes(FIG);
                hold(AX,'on')
                for i = 1:SWEEPS_NUM
                    if Class_1(i)==Fixed_1&&Class_2(i)==Fixed_2
                        P(Class_3(i)) = scatter(AX, RunData.Chnl.NFY(SWEEP_STARTS(i):SWEEP_ENDS(i)), RunData.Chnl.NFX(SWEEP_STARTS(i):SWEEP_ENDS(i)), 8, Scheme.LineColor{Class_3(i)}, 'DisplayName',[RunData.Info.FileName,' ',Legend{Class_3(i)}]);
                        if FittingFlag
                            plot(AX, FitData.Output.FY(SWEEP_STARTS(i):SWEEP_ENDS(i)), FitData.Output.FX(SWEEP_STARTS(i):SWEEP_ENDS(i)), 'Color',Scheme.LineColor{Class_3(i)}, 'LineStyle','-', 'LineWidth',1.5);
                        end
                    end
                end
                title('NFX vs. NFY','FontSize',20);
                xlabel('NFY (-)','FontSize',16);
                ylabel('NFX (-)','FontSize',16);
                set(gca,'FontSize',14);
                set(gca,'XTick',-5:1:5);
                set(gca,'YTick',-5:1:5);
                xlim([-4,4]);
                ylim([-4,4]);
                grid on
                set(AX,'xminorgrid','on','yminorgrid','on')
                pbaspect([1 1 1])               
                %Legend setup
                DeleteIdx = [];
                for iLines = 1:length(P) %Remove entries for unused conditions
                    if ~strcmp(class(P(iLines)),'matlab.graphics.chart.primitive.Scatter')
                        DeleteIdx = [DeleteIdx,iLines];
                    end            
                end
                P(DeleteIdx) = [];
                legend(P)          
            %-----------------------------------------------------------------------------------------------------%                
        end
        %-----------------------------------------------------------------------------------------------------%        
    end
end

function theAxis = scrollsubplot(nrows, ncols, thisPlot)
%SCROLLSUBPLOT Create axes in tiled positions.
%   SCROLLSUBPLOT(m,n,p), breaks the Figure window into
%   an m-by-n matrix of small axes, selects the p-th axes for 
%   for the current plot, and returns the axis handle.  The axes 
%   are counted along the top row of the Figure window, then the
%   second row, etc. 
%   Example:
% 
%       scrollsubplot(3,1,-1), plot(income)
%       scrollsubplot(3,1,1), plot(declared_income)
%       scrollsubplot(3,1,2), plot(tax)
%       scrollsubplot(3,1,3), plot(net_income)
%       scrollsubplot(3,1,4), plot(the_little_extra)
% 
%   plots declared_income on the top third of the window, tax in
%   the middle, and the net_income in the bottom third. Above the
%   top of the figure income is ploted and below the lower edge
%   the_little_extra is to be found. To navigate there is a slider
%   along the right figure edge.
% 
%   SCROLLSUBPLOT now also work with less regular subplot-layouts:
%
%   axs3(1) = scrollsubplot(3,3,1);
%   axs3(2) = scrollsubplot(3,3,3);
%   axs3(3) = scrollsubplot(3,3,[4,7]);
%   axs3(4) = scrollsubplot(3,3,5);
%   axs3(5) = scrollsubplot(3,3,[3,6]);
%   axs3(6) = scrollsubplot(3,3,10);
%   axs3(7) = scrollsubplot(3,3,[8,9,11,12]);
%   axs3(8) = scrollsubplot(3,3,[13,14,15]);
%   for i1 = 1:8,
%     if ishandle(axs3(i1))
%       axes(axs3(i1))
%       imagesc(randn(2+3*i1))
%     end
%   end
%
%   The function works well for regular grids where m,n is constant
%   for all p. When m,n varies there is no guarantee that the steps
%   of the slider is nicely adjusted to the sizes of the
%   subplots, further the slider also responds to mouse-wheel scrolling.
%
%   Differences with SUBPLOT: SCROLLSUBPLOT requires 3 input
%   arguments, no compatibility with subplot(323), no handle as
%   input. Further  PERC_OFFSET_L is decreased from 2*0.09 to 0.07
%   and PERC_OFFSET_R is decreased from 2*0.045 to 0.085. This
%   leaves less space for titles and labels, but give a plaid grid
%   of subplots even outside the visible figure area.
%   
%   Bug/feature when the slider is shifted from its initial
%   position and then extra subplots is added, they get
%   mis-positioned.
%   
%   See also SUBPLOT,

%   Copyright © Bjorn Gustavsson 20050526-2014, Modification/extension of
%   Mathworks subplot.
%   Version 2, modified from version 1 in that scroll now is a subfunction,
%   together with 

persistent maxrownr minrownr

%%% This is a matlab bug(?) that blanks axes that have been rescaled to
%%% accomodate a colorbar. I've not tried to fix it. BG
% we will kill all overlapping siblings if we encounter the mnp
% specifier, else we won't bother to check:
narg = nargin;
% kill_siblings = 0;
create_axis = 1;
delay_destroy = 0;
tol = sqrt(eps);
if narg ~= 3 % not compatible with 3.5, i.e. subplot ==
             % subplot(111) errors out
  error('Wrong number of arguments')
end

%check for encoded format
handle = '';
position = '';

kill_siblings = 1;

% if we recovered an identifier earlier, use it:
if(~isempty(handle))
  
  set(get(0,'CurrentFigure'),'CurrentAxes',handle);
  
elseif(isempty(position))
  % if we haven't recovered position yet, generate it from mnp info:
  if (min(thisPlot) < 1)&&0 % negative Thisplot corresponds to
                            % panels above top row, that is OK
    error('Illegal plot number.')
  else
    % This is the percent offset from the subplot grid of the plotbox.
    PERC_OFFSET_L = 0.07;
    PERC_OFFSET_R = 0.085;
    PERC_OFFSET_B = PERC_OFFSET_L;
    PERC_OFFSET_T = PERC_OFFSET_R;
    if nrows > 2
      PERC_OFFSET_T = 1.1*PERC_OFFSET_T;
      PERC_OFFSET_B = 1.4*PERC_OFFSET_B;
    end
    if ncols > 2
      PERC_OFFSET_L = 0.9*PERC_OFFSET_L;
      PERC_OFFSET_R = 0.9*PERC_OFFSET_R;
    end

    % Subplots version:
    % row = (nrows-1) -fix((thisPlot-1)/ncols)
    % col = rem (thisPlot-1, ncols);
    % Slightly modified to allow for having negative thisPlot
    row = (nrows-1) -floor((thisPlot-1)/ncols);
    col = mod (thisPlot-1, ncols);
    
    % From here on to line 190 essentially identical to SUBPLOT (==untouched)
    
    % For this to work the default axes position must be in normalized coordinates
    if ~strcmp(get(gcf,'defaultaxesunits'),'normalized')
      warning('DefaultAxesUnits not normalized.')
      tmp = axes;
      set(axes,'units','normalized')
      def_pos = get(tmp,'position');
      delete(tmp)
    else
      def_pos = get(gcf,'DefaultAxesPosition')+[-.05 -.05 +.1 +.05];
    end
    col_offset = def_pos(3)*(PERC_OFFSET_L+PERC_OFFSET_R)/ ...
	(ncols-PERC_OFFSET_L-PERC_OFFSET_R);
    row_offset = def_pos(4)*(PERC_OFFSET_B+PERC_OFFSET_T)/ ...
	(nrows-PERC_OFFSET_B-PERC_OFFSET_T);
    totalwidth = def_pos(3) + col_offset;
    totalheight = def_pos(4) + row_offset;
    width = totalwidth/ncols*(max(col)-min(col)+1)-col_offset;
    height = totalheight/nrows*(max(row)-min(row)+1)-row_offset;
    position = [def_pos(1)+min(col)*totalwidth/ncols ...
		def_pos(2)+min(row)*totalheight/nrows ...
		width height];
    if width <= 0.5*totalwidth/ncols
      position(1) = def_pos(1)+min(col)*(def_pos(3)/ncols);
      position(3) = 0.7*(def_pos(3)/ncols)*(max(col)-min(col)+1);
    end
    if height <= 0.5*totalheight/nrows
      position(2) = def_pos(2)+min(row)*(def_pos(4)/nrows);
      position(4) = 0.7*(def_pos(4)/nrows)*(max(row)-min(row)+1);
    end
  end
end

% kill overlapping siblings if mnp specifier was used:
nextstate = get(gcf,'nextplot');
if strncmp(nextstate,'replace',7), nextstate = 'add'; end
if(kill_siblings)
  if delay_destroy
    if nargout %#ok<UNRCH>
      error('Function called with too many output arguments')
    else
      set(gcf,'NextPlot','replace'); return,
    end
  end
  sibs = get(gcf, 'Children');
  got_one = 0;
  for i1 = 1:length(sibs)
    if(strcmp(get(sibs(i1),'Type'),'axes'))
      units = get(sibs(i1),'Units');
      set(sibs(i1),'Units','normalized')
      sibpos = get(sibs(i1),'Position');
      set(sibs(i1),'Units',units);
      intersect = 1;
      if(     (position(1) >= sibpos(1) + sibpos(3)-tol) || ...
	      (sibpos(1) >= position(1) + position(3)-tol) || ...
	      (position(2) >= sibpos(2) + sibpos(4)-tol) || ...
	      (sibpos(2) >= position(2) + position(4)-tol))
	intersect = 0;
      end
      if intersect
	if got_one || any(abs(sibpos - position) > tol)
	  delete(sibs(i1));
	else
	  got_one = 1;
	  set(gcf,'CurrentAxes',sibs(i1));
	  if strcmp(nextstate,'new')
	    create_axis = 1;
	  else
	    create_axis = 0;
	  end
	end
      end
    end
  end
  set(gcf,'NextPlot',nextstate);
end

% create the axis:
if create_axis
  if strcmp(nextstate,'new'), figure, end
  ax = axes('units','normal','Position', position);
  set(ax,'units',get(gcf,'defaultaxesunits'))
else 
  ax = gca; 
end


% return identifier, if requested:
if(nargout > 0)
  theAxis = ax;
end

%% SCROLLSUBPLOT modification part
%
% From here on out set up scrollbar if needed
scroll_hndl = findall(gcf,'Type','uicontrol','Tag','scroll');

ax_indx = findall(gcf,'Type','axes');

if length(ax_indx)==1
  maxrownr = -inf;
  minrownr = inf;
end

%% This is setting up the scroll-slider (invisible)
if isempty(scroll_hndl)
  uicontrol('Units','normalized',...
            'Style','Slider',...
            'Position',[.98,0,.02,1],...
            'Min',0,...
            'Max',1,...
            'Value',1,...
            'visible','off',...
            'Tag','scroll',...
            'Callback',@(scr,event) scroll(1));
  set(gcf,'WindowScrollWheelFcn',@wheelScroll)      
end
% making it visible when needed
if ( nrows*ncols < thisPlot || thisPlot < 1 )
  set(scroll_hndl,'visible','on')
end
scroll(1)


maxrownr = max(maxrownr(:),max(-row(:)));
minrownr = min(minrownr(:),min(-row(:)));

% Adjust the slider step-sizes to account for the number of rows of
% subplots: 
set(scroll_hndl,...
    'sliderstep',[1/nrows 1]./0.8)
%     'sliderstep',[1/nrows
%     1]/(1/((nrows)/max(1,1+maxrownr(:)-minrownr(:)-nrows))))   ORIGINAL
set(scroll_hndl,...
    'value',1)

   function scroll(old_val)
   % SCROLL - Scroll subplots vertically
   %   Used by scrollsubplot.
   % Calling:
   %   scroll(old_val)
   % Set as callback function handle:
   %   
   %   See also SCROLLSUBPLOT
   
   % Copyright Bjorn Gustavsson 20050526
   % Version 2, modified to become a subfunction of scrollsubplot.
   
   %% Scroll the subplot axeses
   %  That is change their position some number of steps up or down
   
   % Get the handle of the scroll-slider handle and all the handle
   % to all the subplot axes
   clbk_ui_hndl = findall(gcf,'Type','uicontrol','Tag','scroll');
   ax_hndl = findall(gcf,'Type','axes');
   
   for i2 = length(ax_hndl):-1:1
     a_pos(i2,:) = get(ax_hndl(i2),'position');
   end
   
   pos_y_range = [min(.07,min(a_pos(:,2))) max(a_pos(:,2) + a_pos(:,4) )+.07-.9];
   
   val = get(clbk_ui_hndl,'value');
   step = ( old_val - val) * diff(pos_y_range);
   
   
   for i2 = 1:length(ax_hndl)
     set(ax_hndl(i2),'position',get(ax_hndl(i2),'position') + [0 step 0 0]);
   end
   
   set(clbk_ui_hndl,'callback',@(scr,event) scroll(val));
   
   end % end scroll

   function wheelScroll(src,evnt) %#ok<INUSL>
   % wheelScroll - mouse-wheel wrapper to scroll-function
   % Calling:
   %   wheelScroll(src,evnt)
   % Set as calback
   %   set(gcf,'WindowScrollWheelFcn',@wheelScroll)
   try
     C = findobj(gcf,'type','uicontrol');
     MaxVal = get(C,'Max');%  = [1]
	 minval = get(C,'Min');% = [0]
	 sliderstep = get(C,'SliderStep'); % = [0.0555556 0.222222]
     sliderstep = sliderstep.*3.367; %Scroll increment adjustment *****OMM     
     oldval = get(C,'value');
     if evnt.VerticalScrollCount > 0 
       set(C,'value',max(minval,oldval-sliderstep(1)/10))
     elseif evnt.VerticalScrollCount < 0 
       set(C,'value',min(MaxVal,oldval+sliderstep(1)/10))
     end
     scroll(oldval)
   catch
   end
   end %wheelScroll

end

function Scheme = ColorSchemeDef()

Scheme = [];
Scheme.CurrentScheme = 'Bright';
Scheme.LineScheme = 'Standard';

switch Scheme.CurrentScheme %Sets theme depending on the scheme chosen - this is already updated before the switch
%-----------------------------------------------------------------------------------------------------%                        
    case 'Bright'
        Scheme.BackColor = [1.00,1.00,1.00]; %Background color - used for figure/axes background, legend color
        Scheme.ForeColor = [0.00,0.00,0.00]; %Foreground color - used for titles, labels, axes ticks
%-----------------------------------------------------------------------------------------------------%                                          
    case 'Dark'
        Scheme.BackColor = [0.15,0.15,0.15];
        Scheme.ForeColor = [1.00,1.00,1.00];           
%-----------------------------------------------------------------------------------------------------%         
end

switch Scheme.LineScheme
%-----------------------------------------------------------------------------------------------------%                    
    case 'Standard'
        %Set plot line colors - Translucent colors only called by line
        %plots. Scatter does not accept a 4th val for alpha, uses separate
        %properties instead for alpha.
        Scheme.LineColor{1,1}  = [0.07,0.62,1.00]; %Blue Opaque
        Scheme.LineColor{1,2}  = [0.07,0.62,1.00,0.40]; %Blue Translucent
        Scheme.LineColor{2,1}  = [1.00,0.00,0.00]; %Red Opaque
        Scheme.LineColor{2,2}  = [1.00,0.00,0.00,0.40]; %Red Translucent
        Scheme.LineColor{3,1}  = [1.00,0.60,0.16]; %Orange Opaque
        Scheme.LineColor{3,2}  = [1.00,0.60,0.16,0.40]; %Orange Translucent
        Scheme.LineColor{4,1}  = [0.72,0.27,1.00]; %Purple Opaque
        Scheme.LineColor{4,2}  = [0.72,0.27,1.00,0.40]; %Purple Translucent
        Scheme.LineColor{5,1}  = [0.00,1.00,0.53]; %Green Opaque
        Scheme.LineColor{5,2}  = [0.00,1.00,0.53,0.40]; %Green Translucent
        Scheme.LineColor{6,1}  = [0.74,0.08,0.20,]; %Burgundy Opaque
        Scheme.LineColor{6,2}  = [0.74,0.08,0.20,0.40]; %Burgundy Translucent 
        Scheme.LineColor{7,1} = [0.74,0.65,0.49]; %Tan Opaque
        Scheme.LineColor{7,2} = [0.74,0.65,0.49,0.40]; %Tan Translucent        
        Scheme.LineColor{8,1}  = [1.00,0.78,0.16]; %Yellow Opaque
        Scheme.LineColor{8,2}  = [1.00,0.78,0.16,0.40]; %Yellow Translucent
        Scheme.LineColor{9,1}  = [0.06,1.00,1.00]; %Cyan Opaque
        Scheme.LineColor{9,2}  = [0.06,1.00,1.00,0.40]; %Cyan Translucent
        Scheme.LineColor{10,1} = [0.00,0.00,1.00]; %Blue Opaque
        Scheme.LineColor{10,2} = [0.00,0.00,1.00,0.40]; %Blue Translucent
        Scheme.LineColor{11,1} = [0.69,0.97,0.58]; %Pale Green Opaque
        Scheme.LineColor{11,2} = [0.69,0.97,0.58,0.40]; %Pale Green Opaque
        Scheme.LineColor{12,1} = [0.98,0.43,0.59]; %Pink Opaque
        Scheme.LineColor{12,2} = [0.98,0.43,0.59,0.40]; %Pink Translucent
        Scheme.LineColor{13,1} = [0.65,0.91,0.98]; %Pale Blue Opaque
        Scheme.LineColor{13,2} = [0.65,0.91,0.98,0.40]; %Pale Blue Translucent
        Scheme.LineColor{14,1}  = [1.00,0.07,0.65]; %Magenta Opaque
        Scheme.LineColor{14,2}  = [1.00,0.07,0.65,0.40]; %Magenta Translucent        
        Scheme.LineColor{15,1} = [0.15,0.53,0.00]; %Dark Green Opaque
        Scheme.LineColor{15,2} = [0.15,0.53,0.00,0.40]; %Dark Green Translucent

        %Repeat same stuff to not run out of lines at 15
        Scheme.LineColor{16,1}  = [0.07,0.62,1.00]; %Blue Opaque
        Scheme.LineColor{16,2}  = [0.07,0.62,1.00,0.40]; %Blue Translucent
        Scheme.LineColor{17,1}  = [1.00,0.00,0.00]; %Red Opaque
        Scheme.LineColor{17,2}  = [1.00,0.00,0.00,0.40]; %Red Translucent
        Scheme.LineColor{18,1}  = [1.00,0.60,0.16]; %Orange Opaque
        Scheme.LineColor{18,2}  = [1.00,0.60,0.16,0.40]; %Orange Translucent
        Scheme.LineColor{19,1}  = [0.72,0.27,1.00]; %Purple Opaque
        Scheme.LineColor{19,2}  = [0.72,0.27,1.00,0.40]; %Purple Translucent
        Scheme.LineColor{20,1}  = [0.00,1.00,0.53]; %Green Opaque
        Scheme.LineColor{20,2}  = [0.00,1.00,0.53,0.40]; %Green Translucent
        Scheme.LineColor{21,1}  = [0.74,0.08,0.20,]; %Burgundy Opaque
        Scheme.LineColor{21,2}  = [0.74,0.08,0.20,0.40]; %Burgundy Translucent        
        Scheme.LineColor{22,1} = [0.74,0.65,0.49]; %Tan Opaque
        Scheme.LineColor{22,2} = [0.74,0.65,0.49,0.40]; %Tan Translucent         
        Scheme.LineColor{23,1}  = [1.00,0.78,0.16]; %Yellow Opaque
        Scheme.LineColor{23,2}  = [1.00,0.78,0.16,0.40]; %Yellow Translucent
        Scheme.LineColor{24,1}  = [0.06,1.00,1.00]; %Cyan Opaque
        Scheme.LineColor{24,2}  = [0.06,1.00,1.00,0.40]; %Cyan Translucent
        Scheme.LineColor{25,1} = [0.00,0.00,1.00]; %Blue Opaque
        Scheme.LineColor{25,2} = [0.00,0.00,1.00,0.40]; %Blue Translucent
        Scheme.LineColor{26,1} = [0.69,0.97,0.58]; %Pale Green Opaque
        Scheme.LineColor{26,2} = [0.69,0.97,0.58,0.40]; %Pale Green Opaque
        Scheme.LineColor{27,1} = [0.98,0.43,0.59]; %Pink Opaque
        Scheme.LineColor{27,2} = [0.98,0.43,0.59,0.40]; %Pink Translucent
        Scheme.LineColor{28,1} = [0.65,0.91,0.98]; %Pale Blue Opaque
        Scheme.LineColor{28,2} = [0.65,0.91,0.98,0.40]; %Pale Blue Translucent
        Scheme.LineColor{29,1}  = [1.00,0.07,0.65]; %Magenta Opaque
        Scheme.LineColor{29,2}  = [1.00,0.07,0.65,0.40]; %Magenta Translucent         
        Scheme.LineColor{30,1} = [0.15,0.53,0.00]; %Dark Green Opaque
        Scheme.LineColor{30,2} = [0.15,0.53,0.00,0.40]; %Dark Green Translucent        
		
		Scheme.LineStyle{1} = '-';
		Scheme.LineStyle{2} = '--';
		Scheme.LineStyle{3} = '-.';
		Scheme.LineStyle{4} = ':';
		Scheme.LineStyle{5} = '-';
		Scheme.LineStyle{6} = '--';
		Scheme.LineStyle{7} = '-.';
		Scheme.LineStyle{8} = ':';
		Scheme.LineStyle{9} = '-';
		Scheme.LineStyle{10} = '--';
		Scheme.LineStyle{11} = '-.';
		Scheme.LineStyle{12} = ':';

		Scheme.Marker{1} = 'o';
		Scheme.Marker{2} = '*';
		Scheme.Marker{3} = 'd';
		Scheme.Marker{4} = 'x';
		Scheme.Marker{5} = '*';
		Scheme.Marker{6} = 'd';
		Scheme.Marker{7} = 'x';
		Scheme.Marker{8} = 'o';
		Scheme.Marker{9} = 'd';
		Scheme.Marker{10} = 'x';
		Scheme.Marker{11} = 'o';
		Scheme.Marker{12} = '*';
        
        %Create custom colormap
        CMTemp1 = parula(64);
        CMTemp1 = CMTemp1(1:end-12,:);
        CMTemp2 = jet(64);
        CMTemp2 = CMTemp2(45:64,:);
        
        Scheme.CMap = [CMTemp1; CMTemp2];
        Scheme.CMap(52,3) = 0.1;
        
        clear CMTemp1 CMTemp2     
end
end

function [FitData] = ModelOptimizer(RunData, FitData, OutputSet, d)
%Solver function for MF 6.1 Models

SA = RunData.Chnl.SA.*(pi./180).*-1; %-1 to ISO
SR = RunData.Chnl.SR;
P  = RunData.Chnl.P.*1000;
IA = RunData.Chnl.IA.*(pi./180);
FZ = RunData.Chnl.FZ.*-1; %-1 to ISO
V  = RunData.Chnl.V.*(0.27778);

FX = RunData.Chnl.FX;
MX = RunData.Chnl.MX;
FY = RunData.Chnl.FY.*-1; %-1 to ISO
MZ = RunData.Chnl.MZ.*-1; %-1 to ISO

switch OutputSet
    %-----------------------------------------------------------------------------------------------------%         
    case 'FX PURE'
        OptIC = FitData.Coeffs.FX;
%             OptFcn = 


    %-----------------------------------------------------------------------------------------------------%   
    case 'FX COMB'


        
        
    %-----------------------------------------------------------------------------------------------------%           
    case 'FY PURE'


        
        
    %-----------------------------------------------------------------------------------------------------%           
    case 'FY COMB'


        
    %-----------------------------------------------------------------------------------------------------%           
    case 'MX PURE'
        global MXobj
        MXobj = MF52_MX(FitData);
        OptIC = [];           
        OptIC(1) = MXobj.QSX1;
        OptIC(2) = MXobj.QSX2;            
        OptIC(3) = MXobj.QSX3;  
        XData(:,1) = FZ;
        XData(:,2) = IA;            
        XData(:,3) = FY;            
        YData = MX; 

        OptFcn = @MX_PURE_Wrapper;
    %-----------------------------------------------------------------------------------------------------%           
    case 'MZ PURE'

        
        
    %-----------------------------------------------------------------------------------------------------%   
    case 'MZ COMB'


        
    %-----------------------------------------------------------------------------------------------------%           
end

Options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',6000,'PlotFcn','optimplotresnorm');
OptCoeffs = lsqcurvefit(OptFcn,OptIC,XData,YData,[],[],Options);

% waitbar(0.1 + 0.45.*(i./SWEEPS_NUM),d,'FY Fitting');

function MX_PURE = MX_PURE_Wrapper(OptIC,XData)
    MXobj.QSX1 = OptIC(1); 
    MXobj.QSX2 = OptIC(2);            
    MXobj.QSX3 = OptIC(3); 
    FZ = XData(:,1);
    IA = XData(:,2);
    FY = XData(:,3);  
    
    MX_PURE = Calculate_Mx(MXobj,FZ,IA,FY);
end

end

