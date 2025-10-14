classdef DNU < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        TabGroup                 matlab.ui.container.TabGroup
        MenuTab                  matlab.ui.container.Tab
        DimensionlessNumbersEdupackLabel  matlab.ui.control.Label
        Image                    matlab.ui.control.Image
        Panel_2                  matlab.ui.container.Panel
        GridLayout2              matlab.ui.container.GridLayout
        FroudeFrWIPButton        matlab.ui.control.Button
        StokesStButton           matlab.ui.control.Button
        WeberWeWIPButton         matlab.ui.control.Button
        DeborahDeWIPButton       matlab.ui.control.Button
        OhnesorgeOhWIPButton     matlab.ui.control.Button
        PcletPeWIPButton         matlab.ui.control.Button
        CapillaryCaButton        matlab.ui.control.Button
        MachMaButton             matlab.ui.control.Button
        Button                   matlab.ui.control.Button
        ReynoldsReButton         matlab.ui.control.Button
        ReTab                    matlab.ui.container.Tab
        ReadMoreButton           matlab.ui.control.Button
        Simulation2Button        matlab.ui.control.Button
        ReAppText                matlab.ui.control.Label
        ReAppTitle               matlab.ui.control.Label
        ReEqnText                matlab.ui.control.Label
        ReEqnLatex               matlab.ui.control.Label
        ReEqnTitle               matlab.ui.control.Label
        ReInterpretText          matlab.ui.control.Label
        ReInterpretationTitle    matlab.ui.control.Label
        ReMenuButton             matlab.ui.control.Button
        Simulation1Button        matlab.ui.control.Button
        ReTheoryText             matlab.ui.control.Label
        ReTheoryTitle            matlab.ui.control.Label
        ReSimTab                 matlab.ui.container.Tab
        Label_4                  matlab.ui.control.Label
        ReSliderLabel            matlab.ui.control.Label
        ReSlider                 matlab.ui.control.Slider
        DropDown                 matlab.ui.control.DropDown
        BackButton               matlab.ui.control.Button
        MenuButton_2             matlab.ui.control.Button
        RunButton                matlab.ui.control.Button
        UIAxes_2                 matlab.ui.control.UIAxes
        UIAxes2                  matlab.ui.control.UIAxes
        ReSim2Tab                matlab.ui.container.Tab
        Label_3                  matlab.ui.control.Label
        StepHeightSlider         matlab.ui.control.Slider
        StepHeightSliderLabel    matlab.ui.control.Label
        StepWidthKnob_2          matlab.ui.control.Knob
        StepWidthKnob_2Label     matlab.ui.control.Label
        DropDown_2               matlab.ui.control.DropDown
        BackButton_4             matlab.ui.control.Button
        MenuButton_7             matlab.ui.control.Button
        ReSlider_2Label          matlab.ui.control.Label
        ReSlider_2               matlab.ui.control.Slider
        RunButton_3              matlab.ui.control.Button
        UIAxes_3                 matlab.ui.control.UIAxes
        UIAxes2_2                matlab.ui.control.UIAxes
        ReSim3Tab                matlab.ui.container.Tab
        Hyperlink                matlab.ui.control.Hyperlink
        MenuButton_10            matlab.ui.control.Button
        BackButton_7             matlab.ui.control.Button
        Re10000Label             matlab.ui.control.Label
        Re1000Label              matlab.ui.control.Label
        Image2                   matlab.ui.control.Image
        WhilethisLabel           matlab.ui.control.Label
        StoTab                   matlab.ui.container.Tab
        StEqnText                matlab.ui.control.Label
        StEqnLatex               matlab.ui.control.Label
        StEquationTitle          matlab.ui.control.Label
        StInterpretText          matlab.ui.control.Label
        StInterpretTitle         matlab.ui.control.Label
        StAppTitle               matlab.ui.control.Label
        StAppText                matlab.ui.control.Label
        StTheorytext             matlab.ui.control.Label
        StTheoryTitle            matlab.ui.control.Label
        SimulationButton         matlab.ui.control.Button
        StMenuButton             matlab.ui.control.Button
        StoSimTab                matlab.ui.container.Tab
        Label_2                  matlab.ui.control.Label
        StokesNumberSlider       matlab.ui.control.Slider
        StokesNumberSliderLabel  matlab.ui.control.Label
        RunButton_2              matlab.ui.control.Button
        BackButton_2             matlab.ui.control.Button
        MenuButton_4             matlab.ui.control.Button
        UIAxes                   matlab.ui.control.UIAxes
        MaTab                    matlab.ui.container.Tab
        ReAppText_2              matlab.ui.control.Label
        ReAppTitle_2             matlab.ui.control.Label
        ReEqnText_2              matlab.ui.control.Label
        ReEqnLatex_2             matlab.ui.control.Label
        ReEqnTitle_2             matlab.ui.control.Label
        ReInterpretText_2        matlab.ui.control.Label
        ReInterpretationTitle_2  matlab.ui.control.Label
        ReTheoryText_2           matlab.ui.control.Label
        ReTheoryTitle_2          matlab.ui.control.Label
        MenuButton_5             matlab.ui.control.Button
        SimulationButton_3       matlab.ui.control.Button
        MaSimTab                 matlab.ui.container.Tab
        ShockwavessthsthLabel    matlab.ui.control.Label
        RunButton_4              matlab.ui.control.Button
        MSlider                  matlab.ui.control.Slider
        MSliderLabel             matlab.ui.control.Label
        BackButton_3             matlab.ui.control.Button
        MenuButton_6             matlab.ui.control.Button
        UIAxes3                  matlab.ui.control.UIAxes
        CaTab                    matlab.ui.container.Tab
        ReEqnText_3              matlab.ui.control.Label
        ReAppText_3              matlab.ui.control.Label
        ReAppTitle_3             matlab.ui.control.Label
        ReEqnLatex_3             matlab.ui.control.Label
        ReEqnTitle_3             matlab.ui.control.Label
        ReInterpretText_3        matlab.ui.control.Label
        ReInterpretationTitle_3  matlab.ui.control.Label
        ReTheoryText_3           matlab.ui.control.Label
        ReTheoryTitle_3          matlab.ui.control.Label
        SimulationButton_4       matlab.ui.control.Button
        ReMenuButton_2           matlab.ui.control.Button
        CaSimTab_2               matlab.ui.container.Tab
        GridLayout3              matlab.ui.container.GridLayout
        Run                      matlab.ui.control.Button
        Panel_3                  matlab.ui.container.Panel
        ReInterpretText_4        matlab.ui.control.Label
        CaLabel                  matlab.ui.control.Label
        Ca                       matlab.ui.control.Slider
        MenuButton_9             matlab.ui.control.Button
        BackButton_6             matlab.ui.control.Button
        CaAxes                   matlab.ui.control.UIAxes
    end

      properties (Access = public)
        SimData % resim
        data
    end
    methods (Access = private)
        
        function results = func(app)
            
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: BackButton, BackButton_4, BackButton_7, 
        % ...and 1 other component
        function ReButtonPushed(app, event)
            app.TabGroup.SelectedTab=app.ReTab;
        end

        % Callback function
        function SliderValueChanging(app, event)
            changingValue = event.Value;
        end

        % Button pushed function: RunButton
        function NewReSim(app, event)
        app.SimData = ReSimLid('compute', app.ReSlider.Value);
ReSimLid('plot', app.SimData, app.DropDown.Value, app.UIAxes_2, app.UIAxes2);
        end

        % Button pushed function: BackButton_2, StokesStButton
        function StoButtonPushed(app, event)
            app.TabGroup.SelectedTab=app.StoTab;
        end

        % Button pushed function: MachMaButton
        function MaButtonPushed(app, event)
            app.TabGroup.SelectedTab=app.MaTab;
        end

        % Button pushed function: Simulation1Button
        function Simulation1ButtonPushed(app, event)
            app.TabGroup.SelectedTab=app.ReSimTab;
        end

        % Button pushed function: MenuButton_10, MenuButton_2, 
        % ...and 8 other components
        function MenuButtonPushed(app, event)
            app.TabGroup.SelectedTab=app.MenuTab;
            cla(app.UIAxes)
            cla(app.UIAxes2)
            cla(app.UIAxes_2)
            cla(app.UIAxes_3)
            cla(app.UIAxes2_2)
            cla(app.UIAxes3)
            cla(app.CaAxes)
        end

        % Callback function
        function RunButton_2Pushed(app, event)
            
        end

        % Button pushed function: RunButton_2
        function RunStokesSimButtonPushed(app, event)
            app.UIAxes.Visible = true;
            StSim(app.UIAxes, app.StokesNumberSlider.Value);
        end

        % Button pushed function: SimulationButton
        function SimulationButtonPushed(app, event)
            app.TabGroup.SelectedTab=app.StoSimTab;
        end

        % Value changed function: StokesNumberSlider
        function StokesNumberSliderValueChanged(app, event)
          
        end

        % Value changed function: DropDown
        function DropDownValueChanged(app, event)
        ReSimLid('plot', app.SimData, app.DropDown.Value, app.UIAxes_2, app.UIAxes2);


        end

        % Button pushed function: RunButton_3
        function RunButton_3Pushed(app, event)
% Compute
app.data = ReSimStep('compute', app.ReSlider_2.Value, ...
                            app.StepWidthKnob_2.Value/100, ...
                            app.StepHeightSlider.Value/100);

% Plot into app axes. Provide axes from the app to avoid use of gca.
ReSimStep('plot', app.data,app.DropDown_2.Value, app.UIAxes_3, app.UIAxes2_2);


        end

        % Value changed function: DropDown_2
        function DropDown_2ValueChanged(app, event)
          ReSimStep('plot', app.data,app.DropDown_2.Value, app.UIAxes_3, app.UIAxes2_2)
            
        end

        % Callback function
        function ReButtonPushed2(app, event)
            app.TabGroup.SelectedTab=app.ReSimTab;
        end

        % Button pushed function: BackButton_3
        function MaBack(app, event)
            app.TabGroup.SelectedTab=app.MaTab;
        end

        % Button pushed function: Simulation2Button
        function ToReSim2(app, event)
            app.TabGroup.SelectedTab=app.ReSim2Tab;
        end

        % Button pushed function: CapillaryCaButton
        function CaButtonPushed(app, event)
            app.TabGroup.SelectedTab=app.CaTab;
        end

        % Callback function
        function ReSimulationButton_2Pushed(app, event)
            
        end

        % Button pushed function: SimulationButton_4
        function CaSimPressed(app, event)
            app.TabGroup.SelectedTab=app.CaSimTab_2;
        end

        % Button pushed function: BackButton_6
        function CaBack(app, event)
            app.TabGroup.SelectedTab=app.CaTab;
        end

        % Image clicked function: Image
        function VisitURL(app, event)
             url = 'https://pubs.aip.org/aip/pof/article/37/9/097173/3364481';  
             web(url, '-browser');          
        end

        % Button pushed function: Run
        function CaRunPushed(app, event)

           CaSim(app.CaAxes,max(app.Ca.Value, 1e-8));            
        end

        % Callback function
        function mach(app, event)
           
        end

        % Button pushed function: RunButton_4
        function RunButton_4Pushed(app, event)
              MSim(app.MSlider.Value, app.UIAxes3)
        end

        % Button pushed function: SimulationButton_3
        function MaSim(app, event)
            app.TabGroup.SelectedTab=app.MaSimTab;
        end

        % Button pushed function: ReadMoreButton
        function ReadMore(app, event)
            app.TabGroup.SelectedTab=app.ReSim3Tab;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 641 452];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Resize = 'off';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.AutoResizeChildren = 'off';
            app.TabGroup.Position = [1 0 640 481];

            % Create MenuTab
            app.MenuTab = uitab(app.TabGroup);
            app.MenuTab.AutoResizeChildren = 'off';
            app.MenuTab.Title = 'Menu';

            % Create Panel_2
            app.Panel_2 = uipanel(app.MenuTab);
            app.Panel_2.AutoResizeChildren = 'off';
            app.Panel_2.Position = [2 3 638 313];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Panel_2);
            app.GridLayout2.ColumnWidth = {'1x', '1x', '1x'};
            app.GridLayout2.RowHeight = {'1x', '1x', '1x'};

            % Create ReynoldsReButton
            app.ReynoldsReButton = uibutton(app.GridLayout2, 'push');
            app.ReynoldsReButton.ButtonPushedFcn = createCallbackFcn(app, @ReButtonPushed, true);
            app.ReynoldsReButton.Layout.Row = 1;
            app.ReynoldsReButton.Layout.Column = 1;
            app.ReynoldsReButton.Text = 'Reynolds (Re)';

            % Create Button
            app.Button = uibutton(app.GridLayout2, 'push');
            app.Button.Layout.Row = 1;
            app.Button.Layout.Column = 2;
            app.Button.Text = '';

            % Create MachMaButton
            app.MachMaButton = uibutton(app.GridLayout2, 'push');
            app.MachMaButton.ButtonPushedFcn = createCallbackFcn(app, @MaButtonPushed, true);
            app.MachMaButton.Layout.Row = 1;
            app.MachMaButton.Layout.Column = 3;
            app.MachMaButton.Text = 'Mach (Ma)';

            % Create CapillaryCaButton
            app.CapillaryCaButton = uibutton(app.GridLayout2, 'push');
            app.CapillaryCaButton.ButtonPushedFcn = createCallbackFcn(app, @CaButtonPushed, true);
            app.CapillaryCaButton.Layout.Row = 2;
            app.CapillaryCaButton.Layout.Column = 1;
            app.CapillaryCaButton.Text = 'Capillary (Ca)';

            % Create PcletPeWIPButton
            app.PcletPeWIPButton = uibutton(app.GridLayout2, 'push');
            app.PcletPeWIPButton.Layout.Row = 2;
            app.PcletPeWIPButton.Layout.Column = 2;
            app.PcletPeWIPButton.Text = 'Péclet (Pe) WIP';

            % Create OhnesorgeOhWIPButton
            app.OhnesorgeOhWIPButton = uibutton(app.GridLayout2, 'push');
            app.OhnesorgeOhWIPButton.Layout.Row = 2;
            app.OhnesorgeOhWIPButton.Layout.Column = 3;
            app.OhnesorgeOhWIPButton.Text = 'Ohnesorge (Oh) WIP';

            % Create DeborahDeWIPButton
            app.DeborahDeWIPButton = uibutton(app.GridLayout2, 'push');
            app.DeborahDeWIPButton.Layout.Row = 3;
            app.DeborahDeWIPButton.Layout.Column = 1;
            app.DeborahDeWIPButton.Text = 'Deborah (De) WIP';

            % Create WeberWeWIPButton
            app.WeberWeWIPButton = uibutton(app.GridLayout2, 'push');
            app.WeberWeWIPButton.Layout.Row = 3;
            app.WeberWeWIPButton.Layout.Column = 2;
            app.WeberWeWIPButton.Text = 'Weber (We) WIP';

            % Create StokesStButton
            app.StokesStButton = uibutton(app.GridLayout2, 'push');
            app.StokesStButton.ButtonPushedFcn = createCallbackFcn(app, @StoButtonPushed, true);
            app.StokesStButton.Layout.Row = 1;
            app.StokesStButton.Layout.Column = 2;
            app.StokesStButton.Text = 'Stokes (St)';

            % Create FroudeFrWIPButton
            app.FroudeFrWIPButton = uibutton(app.GridLayout2, 'push');
            app.FroudeFrWIPButton.Layout.Row = 3;
            app.FroudeFrWIPButton.Layout.Column = 3;
            app.FroudeFrWIPButton.Text = 'Froude (Fr) WIP';

            % Create Image
            app.Image = uiimage(app.MenuTab);
            app.Image.ImageClickedFcn = createCallbackFcn(app, @VisitURL, true);
            app.Image.Position = [5 307 631 144];
            app.Image.ImageSource = fullfile(pathToMLAPP, 'MechBan.png');

            % Create DimensionlessNumbersEdupackLabel
            app.DimensionlessNumbersEdupackLabel = uilabel(app.MenuTab);
            app.DimensionlessNumbersEdupackLabel.BackgroundColor = [0 0 0];
            app.DimensionlessNumbersEdupackLabel.WordWrap = 'on';
            app.DimensionlessNumbersEdupackLabel.FontName = 'Arial';
            app.DimensionlessNumbersEdupackLabel.FontSize = 20;
            app.DimensionlessNumbersEdupackLabel.FontWeight = 'bold';
            app.DimensionlessNumbersEdupackLabel.FontColor = [1 1 1];
            app.DimensionlessNumbersEdupackLabel.Position = [44 329 191 51];
            app.DimensionlessNumbersEdupackLabel.Text = 'Dimensionless Numbers Edupack';

            % Create ReTab
            app.ReTab = uitab(app.TabGroup);
            app.ReTab.AutoResizeChildren = 'off';
            app.ReTab.Title = 'Re';

            % Create ReTheoryTitle
            app.ReTheoryTitle = uilabel(app.ReTab);
            app.ReTheoryTitle.FontWeight = 'bold';
            app.ReTheoryTitle.Position = [26 408 49 22];
            app.ReTheoryTitle.Text = 'Theory';

            % Create ReTheoryText
            app.ReTheoryText = uilabel(app.ReTab);
            app.ReTheoryText.Position = [26 387 582 22];
            app.ReTheoryText.Text = 'The Reynolds number (Re), compares the inertial forces to viscous forces in a fluid flow.';

            % Create Simulation1Button
            app.Simulation1Button = uibutton(app.ReTab, 'push');
            app.Simulation1Button.ButtonPushedFcn = createCallbackFcn(app, @Simulation1ButtonPushed, true);
            app.Simulation1Button.WordWrap = 'on';
            app.Simulation1Button.Position = [188 20 107 23];
            app.Simulation1Button.Text = 'Simulation 1';

            % Create ReMenuButton
            app.ReMenuButton = uibutton(app.ReTab, 'push');
            app.ReMenuButton.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.ReMenuButton.Position = [20 20 100 23];
            app.ReMenuButton.Text = 'Menu';

            % Create ReInterpretationTitle
            app.ReInterpretationTitle = uilabel(app.ReTab);
            app.ReInterpretationTitle.FontWeight = 'bold';
            app.ReInterpretationTitle.Position = [26 366 82 22];
            app.ReInterpretationTitle.Text = 'Interpretation';

            % Create ReInterpretText
            app.ReInterpretText = uilabel(app.ReTab);
            app.ReInterpretText.VerticalAlignment = 'top';
            app.ReInterpretText.WordWrap = 'on';
            app.ReInterpretText.Position = [26 337 582 30];
            app.ReInterpretText.Text = {'A low Re signifies a laminar (smooth) flow whereas a high Re signifies a turbulent (chaotic) flow dominated by eddies which can enhance mixing.'; ''};

            % Create ReEqnTitle
            app.ReEqnTitle = uilabel(app.ReTab);
            app.ReEqnTitle.FontWeight = 'bold';
            app.ReEqnTitle.Position = [26 304 56 22];
            app.ReEqnTitle.Text = 'Equation';

            % Create ReEqnLatex
            app.ReEqnLatex = uilabel(app.ReTab);
            app.ReEqnLatex.HorizontalAlignment = 'center';
            app.ReEqnLatex.WordWrap = 'on';
            app.ReEqnLatex.Interpreter = 'latex';
            app.ReEqnLatex.Position = [29 256 551 48];
            app.ReEqnLatex.Text = '$$Re = \frac{\rho U L}{\mu}$$';

            % Create ReEqnText
            app.ReEqnText = uilabel(app.ReTab);
            app.ReEqnText.VerticalAlignment = 'top';
            app.ReEqnText.WordWrap = 'on';
            app.ReEqnText.Position = [26 227 569 30];
            app.ReEqnText.Text = 'where ρ is the fluid density, U is the characteristic velocity, L is the characteristic length, and μ is the dynamic viscosity of the fluid.';

            % Create ReAppTitle
            app.ReAppTitle = uilabel(app.ReTab);
            app.ReAppTitle.FontWeight = 'bold';
            app.ReAppTitle.Position = [26 194 70 22];
            app.ReAppTitle.Text = 'Application';

            % Create ReAppText
            app.ReAppText = uilabel(app.ReTab);
            app.ReAppText.VerticalAlignment = 'top';
            app.ReAppText.WordWrap = 'on';
            app.ReAppText.Position = [26 89 582 103];
            app.ReAppText.Text = 'The Reynolds number is widely used in industry to predict and control fluid flow behavior in pipes, ducts, and around objects. It helps engineers determine whether flow is laminar or turbulent, influencing design choices for pumps, pipelines, heat exchangers, and aerodynamic surfaces. In chemical and process industries, it guides mixing efficiency, flow metering, and reactor design. In aerospace and automotive applications, it informs drag prediction and flow separation control. By comparing inertial to viscous forces, the Reynolds number enables scalable testing and accurate modeling of real-world fluid systems, ensuring optimal performance and energy efficiency across engineering applications.';

            % Create Simulation2Button
            app.Simulation2Button = uibutton(app.ReTab, 'push');
            app.Simulation2Button.ButtonPushedFcn = createCallbackFcn(app, @ToReSim2, true);
            app.Simulation2Button.Position = [357 20 100 22];
            app.Simulation2Button.Text = 'Simulation 2';

            % Create ReadMoreButton
            app.ReadMoreButton = uibutton(app.ReTab, 'push');
            app.ReadMoreButton.ButtonPushedFcn = createCallbackFcn(app, @ReadMore, true);
            app.ReadMoreButton.Position = [514 20 107 23];
            app.ReadMoreButton.Text = 'Read More...';

            % Create ReSimTab
            app.ReSimTab = uitab(app.TabGroup);
            app.ReSimTab.AutoResizeChildren = 'off';
            app.ReSimTab.Title = 'ReSim';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.ReSimTab);
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes2.Position = [302 174 338 261];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.ReSimTab);
            xlabel(app.UIAxes_2, 'X')
            ylabel(app.UIAxes_2, 'Y')
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_2.Position = [0 174 352 261];

            % Create RunButton
            app.RunButton = uibutton(app.ReSimTab, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @NewReSim, true);
            app.RunButton.VerticalAlignment = 'top';
            app.RunButton.Position = [430 75 100 23];
            app.RunButton.Text = 'Run';

            % Create MenuButton_2
            app.MenuButton_2 = uibutton(app.ReSimTab, 'push');
            app.MenuButton_2.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton_2.Position = [514 20 107 23];
            app.MenuButton_2.Text = 'Menu';

            % Create BackButton
            app.BackButton = uibutton(app.ReSimTab, 'push');
            app.BackButton.ButtonPushedFcn = createCallbackFcn(app, @ReButtonPushed, true);
            app.BackButton.Position = [20 20 100 23];
            app.BackButton.Text = 'Back';

            % Create DropDown
            app.DropDown = uidropdown(app.ReSimTab);
            app.DropDown.Items = {'Pressure', 'Horizontal Velocity', 'Vertical Velocity', 'Absolute Velocity'};
            app.DropDown.ValueChangedFcn = createCallbackFcn(app, @DropDownValueChanged, true);
            app.DropDown.Position = [178 431 100 22];
            app.DropDown.Value = 'Pressure';

            % Create ReSlider
            app.ReSlider = uislider(app.ReSimTab);
            app.ReSlider.Limits = [10 400];
            app.ReSlider.Position = [417 147 150 3];
            app.ReSlider.Value = 10;

            % Create ReSliderLabel
            app.ReSliderLabel = uilabel(app.ReSimTab);
            app.ReSliderLabel.HorizontalAlignment = 'right';
            app.ReSliderLabel.Position = [372 137 25 22];
            app.ReSliderLabel.Text = 'Re';

            % Create Label_4
            app.Label_4 = uilabel(app.ReSimTab);
            app.Label_4.WordWrap = 'on';
            app.Label_4.Position = [41 85 254 74];
            app.Label_4.Text = 'This simulation shows a lid-driven cavity. The container initially contains still fluid. The top lid then begins to move at a constant velocity rightwards. At higher Re, you can see the formation of a secondary vortex.';

            % Create ReSim2Tab
            app.ReSim2Tab = uitab(app.TabGroup);
            app.ReSim2Tab.AutoResizeChildren = 'off';
            app.ReSim2Tab.Title = 'ReSim2';

            % Create UIAxes2_2
            app.UIAxes2_2 = uiaxes(app.ReSim2Tab);
            xlabel(app.UIAxes2_2, 'X')
            ylabel(app.UIAxes2_2, 'Y')
            zlabel(app.UIAxes2_2, 'Z')
            app.UIAxes2_2.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes2_2.Position = [329 174 338 261];

            % Create UIAxes_3
            app.UIAxes_3 = uiaxes(app.ReSim2Tab);
            xlabel(app.UIAxes_3, 'X')
            ylabel(app.UIAxes_3, 'Y')
            zlabel(app.UIAxes_3, 'Z')
            app.UIAxes_3.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes_3.Position = [0 174 352 261];

            % Create RunButton_3
            app.RunButton_3 = uibutton(app.ReSim2Tab, 'push');
            app.RunButton_3.ButtonPushedFcn = createCallbackFcn(app, @RunButton_3Pushed, true);
            app.RunButton_3.Position = [398 63 100 23];
            app.RunButton_3.Text = 'Run';

            % Create ReSlider_2
            app.ReSlider_2 = uislider(app.ReSim2Tab);
            app.ReSlider_2.Limits = [10 100];
            app.ReSlider_2.Position = [392 127 105 3];
            app.ReSlider_2.Value = 10;

            % Create ReSlider_2Label
            app.ReSlider_2Label = uilabel(app.ReSim2Tab);
            app.ReSlider_2Label.HorizontalAlignment = 'right';
            app.ReSlider_2Label.Position = [427 143 25 22];
            app.ReSlider_2Label.Text = 'Re';

            % Create MenuButton_7
            app.MenuButton_7 = uibutton(app.ReSim2Tab, 'push');
            app.MenuButton_7.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton_7.Position = [514 20 107 23];
            app.MenuButton_7.Text = 'Menu';

            % Create BackButton_4
            app.BackButton_4 = uibutton(app.ReSim2Tab, 'push');
            app.BackButton_4.ButtonPushedFcn = createCallbackFcn(app, @ReButtonPushed, true);
            app.BackButton_4.Position = [20 20 100 23];
            app.BackButton_4.Text = 'Back';

            % Create DropDown_2
            app.DropDown_2 = uidropdown(app.ReSim2Tab);
            app.DropDown_2.Items = {'Pressure', 'Horizontal Velocity', 'Vertical Velocity', 'Absolute Velocity'};
            app.DropDown_2.ValueChangedFcn = createCallbackFcn(app, @DropDown_2ValueChanged, true);
            app.DropDown_2.Position = [178 429 100 22];
            app.DropDown_2.Value = 'Pressure';

            % Create StepWidthKnob_2Label
            app.StepWidthKnob_2Label = uilabel(app.ReSim2Tab);
            app.StepWidthKnob_2Label.HorizontalAlignment = 'center';
            app.StepWidthKnob_2Label.Position = [523 54 86 22];
            app.StepWidthKnob_2Label.Text = 'Step Width (%)';

            % Create StepWidthKnob_2
            app.StepWidthKnob_2 = uiknob(app.ReSim2Tab, 'continuous');
            app.StepWidthKnob_2.Limits = [5 75];
            app.StepWidthKnob_2.Position = [534 97 60 60];
            app.StepWidthKnob_2.Value = 5;

            % Create StepHeightSliderLabel
            app.StepHeightSliderLabel = uilabel(app.ReSim2Tab);
            app.StepHeightSliderLabel.HorizontalAlignment = 'right';
            app.StepHeightSliderLabel.Position = [303 33 90 22];
            app.StepHeightSliderLabel.Text = 'Step Height (%)';

            % Create StepHeightSlider
            app.StepHeightSlider = uislider(app.ReSim2Tab);
            app.StepHeightSlider.Limits = [5 75];
            app.StepHeightSlider.MajorTicks = [5 15 25 35 45 55 65 75];
            app.StepHeightSlider.Orientation = 'vertical';
            app.StepHeightSlider.Position = [339 65 3 108];
            app.StepHeightSlider.Value = 5;

            % Create Label_3
            app.Label_3 = uilabel(app.ReSim2Tab);
            app.Label_3.WordWrap = 'on';
            app.Label_3.Position = [41 85 254 74];
            app.Label_3.Text = 'This simulation shows the freestream flow past a backwards facing step. This is can be interpreted as a flow past an obstacle. Higher Re will lead to a larger wake behind the step. The wake contains chaotic vortices.';

            % Create ReSim3Tab
            app.ReSim3Tab = uitab(app.TabGroup);
            app.ReSim3Tab.AutoResizeChildren = 'off';
            app.ReSim3Tab.Title = 'ReSim3';

            % Create WhilethisLabel
            app.WhilethisLabel = uilabel(app.ReSim3Tab);
            app.WhilethisLabel.VerticalAlignment = 'top';
            app.WhilethisLabel.WordWrap = 'on';
            app.WhilethisLabel.Position = [60 293 514 111];
            app.WhilethisLabel.Text = 'Re can also give an approximation of the boundary layer thickness as it is can be approximated as 1/√Re. This implies that higher Re will lead to thinner boundary layers. While this result was formulated for rigid boundary layers, free surface boundary layers can also follow this trend. This is exemplified by the figure below showing the free surface boundary layers of solitary waves (by showing the vorticity contour plot) as studied by ';

            % Create Image2
            app.Image2 = uiimage(app.ReSim3Tab);
            app.Image2.Position = [48 -40 538 461];
            app.Image2.ImageSource = fullfile(pathToMLAPP, 'e5compoundvort2.jpg');

            % Create Re1000Label
            app.Re1000Label = uilabel(app.ReSim3Tab);
            app.Re1000Label.FontSize = 20;
            app.Re1000Label.Interpreter = 'latex';
            app.Re1000Label.Position = [97 230 105 27];
            app.Re1000Label.Text = '$Re=1000$';

            % Create Re10000Label
            app.Re10000Label = uilabel(app.ReSim3Tab);
            app.Re10000Label.FontSize = 20;
            app.Re10000Label.Interpreter = 'latex';
            app.Re10000Label.Position = [95 113 116 27];
            app.Re10000Label.Text = '$Re=10000$';

            % Create BackButton_7
            app.BackButton_7 = uibutton(app.ReSim3Tab, 'push');
            app.BackButton_7.ButtonPushedFcn = createCallbackFcn(app, @ReButtonPushed, true);
            app.BackButton_7.Position = [20 20 100 23];
            app.BackButton_7.Text = 'Back';

            % Create MenuButton_10
            app.MenuButton_10 = uibutton(app.ReSim3Tab, 'push');
            app.MenuButton_10.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton_10.Position = [514 20 107 23];
            app.MenuButton_10.Text = 'Menu';

            % Create Hyperlink
            app.Hyperlink = uihyperlink(app.ReSim3Tab);
            app.Hyperlink.URL = 'https://doi.org/10.1063/5.0287404';
            app.Hyperlink.Position = [315 327 141 22];
            app.Hyperlink.Text = 'Lau and Klettner (2025).';

            % Create StoTab
            app.StoTab = uitab(app.TabGroup);
            app.StoTab.AutoResizeChildren = 'off';
            app.StoTab.Title = 'Sto';

            % Create StMenuButton
            app.StMenuButton = uibutton(app.StoTab, 'push');
            app.StMenuButton.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.StMenuButton.Position = [20 20 100 23];
            app.StMenuButton.Text = 'Menu';

            % Create SimulationButton
            app.SimulationButton = uibutton(app.StoTab, 'push');
            app.SimulationButton.ButtonPushedFcn = createCallbackFcn(app, @SimulationButtonPushed, true);
            app.SimulationButton.Position = [514 20 107 23];
            app.SimulationButton.Text = 'Simulation';

            % Create StTheoryTitle
            app.StTheoryTitle = uilabel(app.StoTab);
            app.StTheoryTitle.FontWeight = 'bold';
            app.StTheoryTitle.Position = [26 408 45 22];
            app.StTheoryTitle.Text = 'Theory';

            % Create StTheorytext
            app.StTheorytext = uilabel(app.StoTab);
            app.StTheorytext.VerticalAlignment = 'top';
            app.StTheorytext.WordWrap = 'on';
            app.StTheorytext.Position = [26 374 585 33];
            app.StTheorytext.Text = 'The Stokes number (St) compares particle inertia to fluid drag forces which describes how much particles will be influenced by the surrounding flow.';

            % Create StAppText
            app.StAppText = uilabel(app.StoTab);
            app.StAppText.WordWrap = 'on';
            app.StAppText.Position = [26 99 580 68];
            app.StAppText.Text = 'The Stokes number predicts how particles behave in a fluid flow, indicating whether they follow the streamlines or move independently. It’s vital in industries like inkjet printing, aerosol design, combustion, and filtration, where particle motion, separation efficiency, and droplet dynamics affect performance, safety, and environmental control.';

            % Create StAppTitle
            app.StAppTitle = uilabel(app.StoTab);
            app.StAppTitle.FontWeight = 'bold';
            app.StAppTitle.Position = [26 167 70 22];
            app.StAppTitle.Text = {'Application'; ''};

            % Create StInterpretTitle
            app.StInterpretTitle = uilabel(app.StoTab);
            app.StInterpretTitle.FontWeight = 'bold';
            app.StInterpretTitle.Position = [26 353 82 22];
            app.StInterpretTitle.Text = 'Interpretation';

            % Create StInterpretText
            app.StInterpretText = uilabel(app.StoTab);
            app.StInterpretText.WordWrap = 'on';
            app.StInterpretText.Position = [26 310 582 44];
            app.StInterpretText.Text = {'Low values of St means droplets closely follow fluid streamlines, but as St increases they tend to ignore the streamlines more.'; ''; ''};

            % Create StEquationTitle
            app.StEquationTitle = uilabel(app.StoTab);
            app.StEquationTitle.FontWeight = 'bold';
            app.StEquationTitle.Position = [26 293 56 22];
            app.StEquationTitle.Text = 'Equation';

            % Create StEqnLatex
            app.StEqnLatex = uilabel(app.StoTab);
            app.StEqnLatex.HorizontalAlignment = 'center';
            app.StEqnLatex.WordWrap = 'on';
            app.StEqnLatex.Interpreter = 'latex';
            app.StEqnLatex.Position = [29 237 551 48];
            app.StEqnLatex.Text = '$$St = \frac{t}{T}$$';

            % Create StEqnText
            app.StEqnText = uilabel(app.StoTab);
            app.StEqnText.VerticalAlignment = 'top';
            app.StEqnText.WordWrap = 'on';
            app.StEqnText.Position = [26 201 592 30];
            app.StEqnText.Text = 'where t is the characteristic time of the particle motion and T is the characteristic time of the fluid flow. The characteristic time may be understood as the time needed to complete a characteristic motion.';

            % Create StoSimTab
            app.StoSimTab = uitab(app.TabGroup);
            app.StoSimTab.AutoResizeChildren = 'off';
            app.StoSimTab.Title = 'StoSim';

            % Create UIAxes
            app.UIAxes = uiaxes(app.StoSimTab);
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Color = [0 0 0];
            app.UIAxes.Box = 'on';
            app.UIAxes.Visible = 'off';
            app.UIAxes.Position = [76 137 496 304];

            % Create MenuButton_4
            app.MenuButton_4 = uibutton(app.StoSimTab, 'push');
            app.MenuButton_4.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton_4.Position = [514 20 107 23];
            app.MenuButton_4.Text = 'Menu';

            % Create BackButton_2
            app.BackButton_2 = uibutton(app.StoSimTab, 'push');
            app.BackButton_2.ButtonPushedFcn = createCallbackFcn(app, @StoButtonPushed, true);
            app.BackButton_2.Position = [20 20 100 23];
            app.BackButton_2.Text = 'Back';

            % Create RunButton_2
            app.RunButton_2 = uibutton(app.StoSimTab, 'push');
            app.RunButton_2.ButtonPushedFcn = createCallbackFcn(app, @RunStokesSimButtonPushed, true);
            app.RunButton_2.Position = [480 103 100 23];
            app.RunButton_2.Text = 'Run';

            % Create StokesNumberSliderLabel
            app.StokesNumberSliderLabel = uilabel(app.StoSimTab);
            app.StokesNumberSliderLabel.HorizontalAlignment = 'right';
            app.StokesNumberSliderLabel.Position = [354 78 88 22];
            app.StokesNumberSliderLabel.Text = 'Stokes Number';

            % Create StokesNumberSlider
            app.StokesNumberSlider = uislider(app.StoSimTab);
            app.StokesNumberSlider.Limits = [0.1 10];
            app.StokesNumberSlider.MajorTicks = [0.1 2 4 6 8 10];
            app.StokesNumberSlider.ValueChangedFcn = createCallbackFcn(app, @StokesNumberSliderValueChanged, true);
            app.StokesNumberSlider.Position = [463 87 150 3];
            app.StokesNumberSlider.Value = 0.1;

            % Create Label_2
            app.Label_2 = uilabel(app.StoSimTab);
            app.Label_2.WordWrap = 'on';
            app.Label_2.Position = [26 52 267 74];
            app.Label_2.Text = 'The particles have an initial horizontal velocity and enter a vortex. You will observe that at low St, the particles have weak inertia and are easily redirected by the vortex. At high St, the particles can persist on their original trajectory.';

            % Create MaTab
            app.MaTab = uitab(app.TabGroup);
            app.MaTab.AutoResizeChildren = 'off';
            app.MaTab.Title = 'Ma';

            % Create SimulationButton_3
            app.SimulationButton_3 = uibutton(app.MaTab, 'push');
            app.SimulationButton_3.ButtonPushedFcn = createCallbackFcn(app, @MaSim, true);
            app.SimulationButton_3.Position = [514 20 107 23];
            app.SimulationButton_3.Text = 'Simulation';

            % Create MenuButton_5
            app.MenuButton_5 = uibutton(app.MaTab, 'push');
            app.MenuButton_5.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton_5.Position = [20 20 100 23];
            app.MenuButton_5.Text = 'Menu';

            % Create ReTheoryTitle_2
            app.ReTheoryTitle_2 = uilabel(app.MaTab);
            app.ReTheoryTitle_2.FontWeight = 'bold';
            app.ReTheoryTitle_2.Position = [26 408 49 22];
            app.ReTheoryTitle_2.Text = 'Theory';

            % Create ReTheoryText_2
            app.ReTheoryText_2 = uilabel(app.MaTab);
            app.ReTheoryText_2.Position = [26 387 582 22];
            app.ReTheoryText_2.Text = 'The Mach number (Ma) compares the fluid velocity to speed of sound in the fluid flow.';

            % Create ReInterpretationTitle_2
            app.ReInterpretationTitle_2 = uilabel(app.MaTab);
            app.ReInterpretationTitle_2.FontWeight = 'bold';
            app.ReInterpretationTitle_2.Position = [26 366 82 22];
            app.ReInterpretationTitle_2.Text = 'Interpretation';

            % Create ReInterpretText_2
            app.ReInterpretText_2 = uilabel(app.MaTab);
            app.ReInterpretText_2.VerticalAlignment = 'top';
            app.ReInterpretText_2.WordWrap = 'on';
            app.ReInterpretText_2.Position = [26 293 582 74];
            app.ReInterpretText_2.Text = {'Values of Ma < 1 means that the flow is subsonic (lower than the speed of sound) where compressibility effects are small.'; 'At a value of Ma = 1, then the flow is sonic (equal to the speed of sound).'; 'Values of Ma > 1 means that it is supersonic (higher than the speed of sound), where shockwaves form and compressibility effects dominate.'; ''};

            % Create ReEqnTitle_2
            app.ReEqnTitle_2 = uilabel(app.MaTab);
            app.ReEqnTitle_2.FontWeight = 'bold';
            app.ReEqnTitle_2.Position = [26 263 56 22];
            app.ReEqnTitle_2.Text = 'Equation';

            % Create ReEqnLatex_2
            app.ReEqnLatex_2 = uilabel(app.MaTab);
            app.ReEqnLatex_2.HorizontalAlignment = 'center';
            app.ReEqnLatex_2.WordWrap = 'on';
            app.ReEqnLatex_2.Interpreter = 'latex';
            app.ReEqnLatex_2.Position = [26 222 551 48];
            app.ReEqnLatex_2.Text = '$$Ma = \frac{|u|}{c}$$';

            % Create ReEqnText_2
            app.ReEqnText_2 = uilabel(app.MaTab);
            app.ReEqnText_2.VerticalAlignment = 'top';
            app.ReEqnText_2.WordWrap = 'on';
            app.ReEqnText_2.Position = [26 201 569 22];
            app.ReEqnText_2.Text = 'where u is the speed of the fluid and c is the speed of sound.';

            % Create ReAppTitle_2
            app.ReAppTitle_2 = uilabel(app.MaTab);
            app.ReAppTitle_2.FontWeight = 'bold';
            app.ReAppTitle_2.Position = [26 167 70 22];
            app.ReAppTitle_2.Text = 'Application';

            % Create ReAppText_2
            app.ReAppText_2 = uilabel(app.MaTab);
            app.ReAppText_2.VerticalAlignment = 'top';
            app.ReAppText_2.WordWrap = 'on';
            app.ReAppText_2.Position = [26 109 582 59];
            app.ReAppText_2.Text = 'The Mach number is crucial in aerospace and automotive industries for designing aircraft, rockets, and high-speed vehicles, influencing aerodynamics, shock wave formation, and structural design under varying flow regimes from subsonic to hypersonic. Geometries that have ideal lift and drag characteristics in the subsonic regime may fare poorly in the supersonic regime due to the formation of shockwaves.';

            % Create MaSimTab
            app.MaSimTab = uitab(app.TabGroup);
            app.MaSimTab.AutoResizeChildren = 'off';
            app.MaSimTab.Title = 'MaSim';

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.MaSimTab);
            xlabel(app.UIAxes3, 'X')
            ylabel(app.UIAxes3, 'Y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes3.Position = [200 89 418 332];

            % Create MenuButton_6
            app.MenuButton_6 = uibutton(app.MaSimTab, 'push');
            app.MenuButton_6.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton_6.Position = [514 20 107 23];
            app.MenuButton_6.Text = 'Menu';

            % Create BackButton_3
            app.BackButton_3 = uibutton(app.MaSimTab, 'push');
            app.BackButton_3.ButtonPushedFcn = createCallbackFcn(app, @MaBack, true);
            app.BackButton_3.Position = [20 20 100 23];
            app.BackButton_3.Text = 'Back';

            % Create MSliderLabel
            app.MSliderLabel = uilabel(app.MaSimTab);
            app.MSliderLabel.HorizontalAlignment = 'right';
            app.MSliderLabel.Position = [215 46 25 22];
            app.MSliderLabel.Text = 'M';

            % Create MSlider
            app.MSlider = uislider(app.MaSimTab);
            app.MSlider.Limits = [0.5 5];
            app.MSlider.MajorTicks = [0.5 1 2 3 4 5];
            app.MSlider.MajorTickLabels = {'0.5', '1', '2', '3', '4', '5'};
            app.MSlider.Position = [261 55 150 3];
            app.MSlider.Value = 0.5;

            % Create RunButton_4
            app.RunButton_4 = uibutton(app.MaSimTab, 'push');
            app.RunButton_4.ButtonPushedFcn = createCallbackFcn(app, @RunButton_4Pushed, true);
            app.RunButton_4.Position = [428 36 68 22];
            app.RunButton_4.Text = 'Run';

            % Create ShockwavessthsthLabel
            app.ShockwavessthsthLabel = uilabel(app.MaSimTab);
            app.ShockwavessthsthLabel.WordWrap = 'on';
            app.ShockwavessthsthLabel.Position = [26 213 196 192];
            app.ShockwavessthsthLabel.Text = 'This simulation shows the freestream flow over a flat plate. At subsonic regimes (M<1), the Blasius solution for flow over a flat plate is recovered. However, at supersonic regimes, shock waves form which can be observed as drastic and discontinuous changes in the Mach number. While M varies smoothly in subsonic flows, M varies discontinuously through shock waves in supersonic flows.';

            % Create CaTab
            app.CaTab = uitab(app.TabGroup);
            app.CaTab.AutoResizeChildren = 'off';
            app.CaTab.Title = 'Ca';

            % Create ReMenuButton_2
            app.ReMenuButton_2 = uibutton(app.CaTab, 'push');
            app.ReMenuButton_2.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.ReMenuButton_2.Position = [20 20 100 23];
            app.ReMenuButton_2.Text = 'Menu';

            % Create SimulationButton_4
            app.SimulationButton_4 = uibutton(app.CaTab, 'push');
            app.SimulationButton_4.ButtonPushedFcn = createCallbackFcn(app, @CaSimPressed, true);
            app.SimulationButton_4.WordWrap = 'on';
            app.SimulationButton_4.Position = [514 20 107 23];
            app.SimulationButton_4.Text = 'Simulation';

            % Create ReTheoryTitle_3
            app.ReTheoryTitle_3 = uilabel(app.CaTab);
            app.ReTheoryTitle_3.FontWeight = 'bold';
            app.ReTheoryTitle_3.Position = [26 408 49 22];
            app.ReTheoryTitle_3.Text = 'Theory';

            % Create ReTheoryText_3
            app.ReTheoryText_3 = uilabel(app.CaTab);
            app.ReTheoryText_3.WordWrap = 'on';
            app.ReTheoryText_3.Position = [26 379 582 30];
            app.ReTheoryText_3.Text = 'The Capillary Number (Ca), compares the viscous forces to surface tension forces, which determines how easily a fluid can deform/spread due to viscous effects overcoming surface tension.';

            % Create ReInterpretationTitle_3
            app.ReInterpretationTitle_3 = uilabel(app.CaTab);
            app.ReInterpretationTitle_3.FontWeight = 'bold';
            app.ReInterpretationTitle_3.Position = [26 358 82 22];
            app.ReInterpretationTitle_3.Text = 'Interpretation';

            % Create ReInterpretText_3
            app.ReInterpretText_3 = uilabel(app.CaTab);
            app.ReInterpretText_3.VerticalAlignment = 'top';
            app.ReInterpretText_3.WordWrap = 'on';
            app.ReInterpretText_3.Position = [26 297 582 59];
            app.ReInterpretText_3.Text = {'Low values of Ca means that droplets retain their shape with little deformation.'; 'As Ca increases, droplets experience deformation but still remain intact.'; 'Eventually, at a high enough Ca (depending on the fluid), surface tension within droplets can no longer maintain its shape, causing spreading or film formation.'; ''};

            % Create ReEqnTitle_3
            app.ReEqnTitle_3 = uilabel(app.CaTab);
            app.ReEqnTitle_3.FontWeight = 'bold';
            app.ReEqnTitle_3.Position = [26 272 56 22];
            app.ReEqnTitle_3.Text = 'Equation';

            % Create ReEqnLatex_3
            app.ReEqnLatex_3 = uilabel(app.CaTab);
            app.ReEqnLatex_3.HorizontalAlignment = 'center';
            app.ReEqnLatex_3.WordWrap = 'on';
            app.ReEqnLatex_3.Interpreter = 'latex';
            app.ReEqnLatex_3.Position = [29 218 551 48];
            app.ReEqnLatex_3.Text = '$$Ca = \frac{\mu U}{\sigma}$$';

            % Create ReAppTitle_3
            app.ReAppTitle_3 = uilabel(app.CaTab);
            app.ReAppTitle_3.FontWeight = 'bold';
            app.ReAppTitle_3.Position = [28 137 70 22];
            app.ReAppTitle_3.Text = 'Application';

            % Create ReAppText_3
            app.ReAppText_3 = uilabel(app.CaTab);
            app.ReAppText_3.VerticalAlignment = 'top';
            app.ReAppText_3.WordWrap = 'on';
            app.ReAppText_3.Position = [28 65 582 59];
            app.ReAppText_3.Text = 'The Capillary number quantifies the balance between viscous and surface tension forces in fluids. It’s crucial in industries like oil recovery, inkjet printing, coating, and microfluidics, where controlling droplet formation, wetting, and flow in small channels determines process efficiency, product quality, and precise fluid manipulation';

            % Create ReEqnText_3
            app.ReEqnText_3 = uilabel(app.CaTab);
            app.ReEqnText_3.VerticalAlignment = 'top';
            app.ReEqnText_3.WordWrap = 'on';
            app.ReEqnText_3.Position = [20 188 569 35];
            app.ReEqnText_3.Text = 'where μ is the dynamic viscosity of the fluid, U is the droplet velocity, and σ is the surface tension of the droplet.';

            % Create CaSimTab_2
            app.CaSimTab_2 = uitab(app.TabGroup);
            app.CaSimTab_2.AutoResizeChildren = 'off';
            app.CaSimTab_2.Title = 'CaSim';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.CaSimTab_2);
            app.GridLayout3.ColumnWidth = {'1x', '1x', '1x'};
            app.GridLayout3.RowHeight = {'1x', 23};

            % Create CaAxes
            app.CaAxes = uiaxes(app.GridLayout3);
            xlabel(app.CaAxes, 'X')
            ylabel(app.CaAxes, 'Y')
            zlabel(app.CaAxes, 'Z')
            app.CaAxes.Layout.Row = 1;
            app.CaAxes.Layout.Column = [2 3];

            % Create BackButton_6
            app.BackButton_6 = uibutton(app.GridLayout3, 'push');
            app.BackButton_6.ButtonPushedFcn = createCallbackFcn(app, @CaBack, true);
            app.BackButton_6.Layout.Row = 2;
            app.BackButton_6.Layout.Column = 1;
            app.BackButton_6.Text = 'Back';

            % Create MenuButton_9
            app.MenuButton_9 = uibutton(app.GridLayout3, 'push');
            app.MenuButton_9.ButtonPushedFcn = createCallbackFcn(app, @MenuButtonPushed, true);
            app.MenuButton_9.Layout.Row = 2;
            app.MenuButton_9.Layout.Column = 3;
            app.MenuButton_9.Text = 'Menu';

            % Create Panel_3
            app.Panel_3 = uipanel(app.GridLayout3);
            app.Panel_3.AutoResizeChildren = 'off';
            app.Panel_3.BorderType = 'none';
            app.Panel_3.Layout.Row = 1;
            app.Panel_3.Layout.Column = 1;

            % Create Ca
            app.Ca = uislider(app.Panel_3);
            app.Ca.Limits = [0 10000];
            app.Ca.Orientation = 'vertical';
            app.Ca.Position = [95 18 10 248];

            % Create CaLabel
            app.CaLabel = uilabel(app.Panel_3);
            app.CaLabel.HorizontalAlignment = 'right';
            app.CaLabel.FontSize = 18;
            app.CaLabel.Position = [37 167 44 24];
            app.CaLabel.Text = 'Ca =';

            % Create ReInterpretText_4
            app.ReInterpretText_4 = uilabel(app.Panel_3);
            app.ReInterpretText_4.VerticalAlignment = 'top';
            app.ReInterpretText_4.WordWrap = 'on';
            app.ReInterpretText_4.Position = [8 294 188 89];
            app.ReInterpretText_4.Text = 'This simulation shows a droplet in a shearing flow. The top moves leftwards but the bottom moves righwards. At low Ca, the droplet is barely affected. At high Ca, the droplet is deformed.';

            % Create Run
            app.Run = uibutton(app.GridLayout3, 'push');
            app.Run.ButtonPushedFcn = createCallbackFcn(app, @CaRunPushed, true);
            app.Run.Layout.Row = 2;
            app.Run.Layout.Column = 2;
            app.Run.Text = 'Run';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DNU

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end