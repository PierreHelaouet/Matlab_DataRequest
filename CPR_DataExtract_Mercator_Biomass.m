%% Program extraction CPR data for Mercatpr
% Contact: Alexandre Mignot <amignot@mercator-ocean.fr>
% Description: mesozooplanbkton biomass per longhurst areas
% le 01/02/2023
% Pierre helaouet
%% Load database
load Database_2022_ENV_06Dec22
%% Load shapefile
S_Longhurst = shaperead('D:\Matlab\Data\Longhurst areas\Longhurst_world_v4_2010.shp');
%% Attribute code to database using shapefile
% Note: test were conducted to evaluate potential ways to use 'inpolygon'
% effiently. Results show that parallel calculation is the way to go.
% See 'Matlab_InpolygonTest.xlsx'

% Transfer data to GPU
xq = gpuArray(Database_2022.SpaceTime.Longitude);
yq = gpuArray(Database_2022.SpaceTime.Latitude);

for ii = 1:length(S_Longhurst)
    % Display current areas
    disp(['Area: ' num2str(ii) ' / 54 :' S_Longhurst(ii).ProvDescr])
    % Set selected polygon in GPU
    clear xv yv
    xv = gpuArray(S_Longhurst(ii).X);
    yv = gpuArray(S_Longhurst(ii).Y);

    % Run calculation on GPU
    clear IN
    IN = inpolygon(xq, yq, xv, yv);

    % Brings results to Workspace
    IN_Local = gather(IN);
    if sum(IN_Local) > 0
        S_Longhurst(ii).Is_CPR = 1;
        S_Longhurst(ii).Idx_CPR = IN_Local;
    else
        S_Longhurst(ii).Is_CPR = 0;
        S_Longhurst(ii).Idx_CPR = [];
    end
end
clearvars -except Database_2022 S_Longhurst
%% Save 
cd('D:\Matlab\Data\matrices analyses\Mercator')
save('CPRtoMercator_Biomass_1', 'Database_2022', 'S_Longhurst');
%% Extract Idx to trim shapefile by removing areas without data
V_Idx_S = logical(extractfield(S_Longhurst, 'Is_CPR'));
S_Longhurst(~V_Idx_S) = [];
clearvars -except Database_2022 S_Longhurst
%% Build CPR subset: step 1: add sample info and spacetime
for ii = 1:length(S_Longhurst)
    S_CPR(ii, 1).Sample = Database_2022.Sample(S_Longhurst(ii).Idx_CPR, :);
    S_CPR(ii, 1).SpaceTime = Database_2022.SpaceTime(S_Longhurst(ii).Idx_CPR, :);
end
%% Map
figure
axm = axesm('MapProjection','Lambertstd');
for ii =  1:length(S_Longhurst)
    scatterm(S_CPR(ii).SpaceTime.Latitude, S_CPR(ii).SpaceTime.Longitude, 10, rgb('DarkRed'));
    hold on
    h_patch(ii) = patchm(S_Longhurst(ii).Y, S_Longhurst(ii).X, 'k');
    h_patch(ii).FaceColor ='none';
    h_patch(ii).EdgeColor = 'k';
end
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5]);
%% Trim areas using sampling effort
V_YearMonth = [...
    reshape(repmat(1946:2022, 12, 1), length(1946:2022)*12, 1), ...
    repmat((1:12)', length(1946:2022), 1) ...
    ];

for ii = 1:length(S_CPR)
    T_stats_E = grpstats(S_CPR(ii).SpaceTime,{'Year','Month'});
    S_CPR(ii).SpaceTime
end


[~, Idx_Year] = ismember(T_stats_E.Year, V_YearMonth(:, 1));
V_RegYear = zeros(length(V_YearMonth), 1);
V_RegYear(Idx_Year) = T_stats_E.Year;

[~, Idx_Month] = ismember(T_stats_E.Month, V_YearMonth(:, 2));
V_RegMonth = zeros(length(V_YearMonth), 1);
V_RegMonth(Idx_Month) = T_stats_E.Month;




for ii = 1:length(S_CPR)
    V_Count(ii,1) = sum(S_Longhurst(ii).Idx_CPR);
end
%% Extract data with size
T_Traits = readtable('D:\Science\ICPR\iCPR_Database\iCPR_Functional\Zoo_Traits.xlsx');
T_Size = T_Traits(:, [2,10:12]);
%% Match size with taxa list




    S_CPR(ii).Data.Data_EyeC = Database_2022.Data.Data_EyeC(S_Longhurst(ii).Idx_CPR, :);
    S_CPR(ii).Data.Data_Trav = Database_2022.Data.Data_Trav(S_Longhurst(ii).Idx_CPR, :);
    S_CPR(ii).List.List_EyeC = Database_2022.List.List_EyeC;
    S_CPR(ii).List.List_Trav = Database_2022.List.List_Trav;
    S_CPR(ii).Taxo.Taxo_EyeC = Database_2022.Taxo.Taxo_EyeC;
    S_CPR(ii).taxo.Taxo_Trav = Database_2022.Taxo.Taxo_Trav;
