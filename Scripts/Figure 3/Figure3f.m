%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3f -- Using video data to predict projections onto CDramp on 
% single trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc,close all  
%% Set paths

% Base path for where code lives
basepth = 'C:\Code';

% add paths
utilspth = [basepth '\Uninstructed-Movements\For sharing'];
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'Utils')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % (1) all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % (2) right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % (3) left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % (4) error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % (5) error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};              % (6) no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};              % (7) no left, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % (8) all hits, no stim, aw off


params.tmin = -2.5;     % min time (in s) relative to alignEvent 
params.tmax = 2.5;      % max time (in s) relative to alignEvent
params.dt = 1/100;      % size of time bin

% smooth with causal gaussian kernel
params.smooth = 15;

% Sorted unit qualities to use (can specify whether you only want to use
% single units or multi units)
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% Kinematic features that you want to load
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
params.bctype = 'reflect'; % options are : reflect  zeropad  none
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';

meta = [];

% Scripts for loading data from each animal  
meta = loadJEB6_ALMVideo(meta,datapth);     
meta = loadJEB7_ALMVideo(meta,datapth);              
meta = loadEKH1_ALMVideo(meta,datapth);     
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);     
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);  
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating coding dimensions----')
cond2use = [2 3];   % Conditions to use to calculate CDchoice (sometimes it may say CDlate, these are the same thing) and CDaction
                    % right hits, left hits (with reference to to PARAMS.CONDITION)
rampcond = 8;       % Condition to use to calculate CDramp (with reference to PARAMS.CONDITION)
cond2proj = 2:7;    % Conditions that you want to project onto the CDs
                    % right hits, left hits, right miss, left miss, right no, left no (corresponding to PARAMS.CONDITION)
cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
regr = getCodingDimensions_2afc(obj,params,cond2use,rampcond, cond2proj);   % Calculate the coding dimensions

disp('----Projecting single trials onto CDRamp----')
cd = 'ramping';     % CD that you want to project single trials onto 
smooth = 60;        % How much you want to smooth the single trial projections by
regr = getSingleTrialProjs(regr,obj,cd,smooth);     % Will add a field called 'SingleProj' to 'regr' variable
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
clearvars -except obj meta kin regr params me
%% Set parameters for CDramp decoding from kinematic features%%

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)

% Set how many time bins before and after current time bin you want to use for decoding
par.pre=6;          % number of time bins prior to output used for decoding
par.post=0;         % number of time bins after output used for decoding
par.dt = params(1).dt;  % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using
par.post_s = par.post .* params(1).dt;

% Set parameters for what time points within the trial you want to use for decoding
trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent)+0.4);
start = find(obj(1).time>trialstart,1,'first');
go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time>go,1,'first');
par.timerange = start:stop;

% cross val folds
par.nFolds = 4;

% data sets
par.train = 0.6;        % fraction of trials to be used for the training data 
par.test = 1 - par.train;

% kinematic features to use to decode
par.feats = {'motion','nos','jaw'};     % Any kinematic features that contain those strings will be used
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});

% what trials you want to include in the decoding analysis
par.cond2use = [2 3]; % R and L hits, DR trials only (not WC, ~early) -- with reference to params.condition

% model parameters
par.regularize = 1;                 % if 0, linear regression. if 1, ridge regression
par.lambdas = logspace(-3,3,20);    % regularization parameters to search over
%% Do the decoding (predict CDramp on single trials from kinematic features)

close all

for sessix = 1:numel(meta)                          % For each session...
    disp(['Decoding for session ' ' ' num2str(sessix) ' / ' num2str(numel(meta))])
    
    % Organize the kinematic and neural data into the proper format
    [X,Y,par] = preparePredictorsRegressors_v2(par, sessix, kin, regr,params);

    % Search over regularization parameters
    for ilambda = 1:numel(par.lambdas)            % For every lambda value in your range...
        lambda = par.lambdas(ilambda);
        cvmdl{ilambda} = fitrlinear(X.train,Y.train,'Learner','svm','KFold',par.nFolds,'Regularization','ridge','Lambda',lambda);   % Fit a cross validated model with this lambda
        loss(ilambda) = cvmdl{ilambda}.kfoldLoss;       % Get the loss for the cvmdl with this lambda
    end
    [~,bestmodel] = min(loss);                          % Find the best lambda
    mdl = cvmdl{bestmodel};                             % we now have the best lambda, and then trained cvmodel with that lambda,
    % can predict test data using best cv mdl
    for i = 1:par.nFolds
        mdl_loss(i) = mdl.Trained{i}.loss(X.train,Y.train);
    end
    [~,bestmodel] = min(mdl_loss);
    testmdl = mdl.Trained{bestmodel}; % we now have the best lambda, and the trained model with that lambda,
    
    pred = testmdl.predict(X.test);   % use left-out test data for the prediction
    
    % Save the predictor coefficients from the testmdl
    loadings(:,sessix) = testmdl.Beta;  % Average the coefficients for each predictor term across folds; save these for each time point

    y = reshape(Y.test,Y.size.test(1),Y.size.test(2));      % original input (neural) data (standardized)
    yhat = reshape(pred,Y.size.test(1),Y.size.test(2));     % predicted input (neural) data

    % find which trials were Rhit and which were Lhit in test set
    Rhit_trials = ismember(par.trials.test,par.trials.Rhit); 
    Lhit_trials = ismember(par.trials.test,par.trials.Lhit); 

    % neural data, split into trial type
    trueVals.Rhit{sessix} = y(:,Rhit_trials);
    trueVals.Lhit{sessix} = y(:,Lhit_trials);

    % predicted neural data, split into ground truth
    modelpred.Rhit{sessix} = yhat(:,Rhit_trials);
    modelpred.Lhit{sessix} = yhat(:,Lhit_trials);
end

disp('---FINISHED DECODING FOR ALL SESSIONS---')
clearvars -except datapth kin me meta obj params regr nSessions exsess modelpred trueVals par loadings
%% Baseline subtract CDramp

% Times that you want to use as the baseline (trialstart to the sample tone onset)
trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

cond2use = {'Rhit','Lhit'};
for sessix = 1:length(meta)             % For every session...
    for c = 1:length(cond2use)          % For both right and left hit trials
        cond = cond2use{c};             % Get the current condition
        curr = modelpred.(cond){sessix};                % Get the predicted neural data from all trials of this condition from this session 
        presampcurr = curr(start:stop,:);               % Get the presample values (time x trials)
        presampcurr = mean(presampcurr,1,'omitnan');    % Take the average across the presample period (1 x trials)
        presampcurr = mean(presampcurr,'omitnan');      % Take the average across trials (1 x 1)

        modelpred.(cond){sessix} = modelpred.(cond){sessix}-presampcurr;    % Subtract that baseline value from all values

        curr = trueVals.(cond){sessix};                 % Do the same thing for ground truth neural data
        presampcurr = curr(start:stop,:);
        presampcurr = mean(presampcurr,1,'omitnan');
        presampcurr = mean(presampcurr,'omitnan');

        trueVals.(cond){sessix} = trueVals.(cond){sessix}-presampcurr;
    end
end
%% Plot example of motion energy and projections onto CDramp, averaged across right and left trials
exsess = 19;
colors = getColors();

nSessions = numel(meta);
cond2plot = [2 3];          % Which conditions you want to include in the plot (will be grouped together after) 
MEix = find(strcmp(kin(1).featLeg,'motion_energy'));    % Get the index within 'kin' that corresponds to 'motion_energy'

% time parameters
times.trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
times.startix = find(obj(1).time>times.trialstart,1,'first');
times.samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
times.stopix = find(obj(1).time<times.samp,1,'last');

% plotting parameters
cols = {colors.rhit,colors.lhit};
alph = 0.2;
sm = 21;
for sessix = 1:numel(meta)                  % For every session...
    condME = [];
    condRamp = [];
    for c = 1:length(cond2plot)             % For each condition...
        cond = cond2plot(c);                
        condtrix = params(sessix).trialid{cond};            % Trials that belong to that condition
        temp = squeeze(kin(sessix).dat(:,condtrix,MEix));   % Get motion energy from those trials (time x nTrials in condition)
        condME = [condME,temp];                             % Concatenate motion energy traces across conditions (time x [nTrials in cond1 + nTrials in cond2])

         if c==1
            dir = 'Rhit';
        else
            dir = 'Lhit';
        end
        temp = trueVals.(dir){sessix};                      % Concatenate CDramp projections across trials from both conditions
        condRamp = [condRamp,temp];
    end

    % Baseline subtract ME (baseline is 1st percentile)
    pctME = prctile(condME,1);
    condME = condME-pctME;
    
    nTrials = size(condME,2);
    % Plot average motion energy
    ax1 = subplot(2,1,1);
    toplot = mean(mySmooth(condME,31),2,'omitnan');
    err = 1.96*(std(mySmooth(condME,31),0,2,'omitnan')./sqrt(nTrials));
    ax = gca;
    shadedErrorBar(obj(sessix).time,mySmooth(toplot,sm),err,{'Color',[0.2 0.2 0.2],'LineWidth',2},alph,ax); hold on;
    % Plot average projection onto CDramp
    ax2 = subplot(2,1,2);
    toplot = mean(mySmooth(condRamp(par.timerange,:),21),2,'omitnan');
    err = 1.96*(std(mySmooth(condRamp(par.timerange,:),21),0,2,'omitnan')./sqrt(nTrials));
    ax = gca;
    shadedErrorBar(obj(sessix).time(par.timerange),toplot,err,{'Color',[0.2 0.2 0.2],'LineWidth',2},alph,ax); hold on;
end
xlabel(ax1,'Time from go cue (s)')
ylabel(ax1,'Motion energy (a.u.)')
xline(ax1,0,'k--','LineWidth',1)
xline(ax1,-0.9,'k--','LineWidth',1)
xline(ax1,-2.2,'k--','LineWidth',1)
xlim(ax1,[-2.3 0])
set(ax1,'TickDir','out'); 

xlabel(ax2,'Time from go cue (s)')
ylabel(ax2,'Projection onto CDramp (a.u.)')
xline(ax2,0,'k--','LineWidth',1)
xline(ax2,-0.9,'k--','LineWidth',1)
xline(ax2,-2.2,'k--','LineWidth',1)
xlim(ax2,[-2.3 0])
set(ax2,'TickDir','out'); 
%% Make heatmaps for a single session showing projections onto CDramp across trials and predicted CDRamp

sm = 15;

load('C:\Code\Uninstructed-Movements\LeftRightDiverging_Colormap.mat')

% Times that you want to use to sort CDramp
del = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>del,1,'first');
resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
stop = find(obj(1).time<resp,1,'last');

cond2plot = {'Lhit','Rhit'};
for sessix = 1:length(meta)                                                                   % For each session...
    figure();
    cnt = 0;
    tempTrue = [];
    tempPred = [];
    % Combine the true values for CDramp and the model predicted
    % values across conditions
    % tempTrue = (time x [num left trials + num right trials])
    for c = 1:length(cond2plot)                                                 % For left and right trials...
        cond = cond2plot{c};
        currTrue = trueVals.(cond){sessix};                                     % Get the true single trial CDramp projections for that condition and session
        tempTrue = [tempTrue,currTrue];
        currPred = modelpred.(cond){sessix};                                    % Get the model predicted single trial CDramp projections
        tempPred = [tempPred,currPred];
    end
    [~,sortix] = sort(mean(tempTrue(start:stop,:),1,'omitnan'),'descend');      % Sort the true projections by average magnitude during the delay period
    True2plot = mySmooth(tempTrue(:,sortix),31);
    Pred2plot = tempPred(:,sortix);                                             % Sort the model predictions in the same order

    nTrials = size(True2plot,2);                                                % Total number of trials that are being plotted
    % Heatmap of true neural data 
    ax1 = subplot(1,2,1);                                                       % Plot true CDramp data on left subplot
    imagesc(obj(sessix).time(par.timerange),1:nTrials,True2plot'); hold on     

    % Heatmap of predicted neural data
    ax2 = subplot(1,2,2);
    imagesc(obj(sessix).time(par.timerange),1:nTrials,mySmooth(Pred2plot,21)'); hold on     
    title(ax1,'CDramping - recorded neural data')
    colorbar(ax1)
    colormap(LeftRightDiverging_Colormap)
    clim(ax1,[-3 3])
    xlabel(ax1,'Time from go cue (s)')
    xline(ax1,0,'k--','LineWidth',1)
    xline(ax1,-0.9,'k--','LineWidth',1)
    xline(ax1,-2.2,'k--','LineWidth',1)
    xlim(ax1,[-2.5 0])
    set(gca,'TickDir','out');

    title(ax2,'Model prediction')
    xlabel(ax2,'Time from go cue (s)')
    xline(ax2,0,'k--','LineWidth',1)
    xline(ax2,-0.9,'k--','LineWidth',1)
    xline(ax2,-2.2,'k--','LineWidth',1)
    xlim(ax2,[-2.5 0])
    set(gca,'TickDir','out');

    colorbar(ax2)
    colormap(LeftRightDiverging_Colormap)
    clim(ax2,[-2.5 2.5])

    sgtitle(['Example session:  ' meta(sessix).anm ' ' meta(sessix).date])
end
%% Example plots by session for relating predicted and true CDTrialType
delR2_ALL = [];

for sessix = 1:length(meta)
    %%% Plot a scatter plot for a single session of true CDramp and predicted CDramp for each trial
    %%% Each dot = an average value of CDramp during the delay period
    figure();
    subplot(1,2,2)
    tempR2 = Scatter_ModelPred_TrueCDTrialType(trueVals, modelpred, sessix, start, stop,meta,'no');
    % Save R2 value for that session
    delR2_ALL = [delR2_ALL, tempR2];

    % Calculate averages and standard deviation for true CD and predicted CD  this session
    [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix);

    colors = getColors();
    alph  = 0.2;

    %%% Plot an example session of CDramp prediction vs true value
    subplot(1,2,1)
    plotExampleCDRamp_Pred(colors, obj, par, meta, avgCD, stdCD, sessix, trueVals,alph, tempR2,'no');
end
%% Plot bar plot to show average R2 values across sessions
colors = getColors();
delR2_ALL = abs(delR2_ALL);

nSessions =length(meta);

markerSize = 60;
figure();
% Bar plot = mean R2 value across all sessions
b = bar(mean(delR2_ALL),'FaceColor',colors.afc); hold on;               
ix2plot = 1:nSessions;
ix2plot(exsess) = [];
% Plot the R2 values as a dot for each session
scatter(ones(nSessions-1,1),delR2_ALL(ix2plot),markerSize,'filled','o','MarkerFaceColor',...,
    'k','XJitter','randn','XJitterWidth',0.25); hold on;
% Plot the R2 value as a different colored dot for the example session
scatter(1,delR2_ALL(exsess),markerSize,'o','MarkerEdgeColor',colors.afc)
% Plot error bar
errorbar(b.XEndPoints,mean(delR2_ALL,'omitnan'),std(delR2_ALL,'omitnan'),'LineStyle','none','Color','k','LineWidth',1)
ylim([0 0.9])
set(gca,'TickDir','out');
ax = gca;
ax.FontSize = 16;
title(['CDramp: Ex session = ' meta(exsess).anm meta(exsess).date])
%% FUNCTIONS
function [avgCD,stdCD] = getAvgStd(trueVals,modelpred,sessix)
temp = [trueVals.Rhit{sessix},trueVals.Lhit{sessix}];
avgCD.allhits.true = mean(mySmooth(temp,81),2,'omitnan');
stdCD.allhits.true = mean(mySmooth(temp,81),2,'omitnan');
avgCD.allhits.nTrials = size(temp,2);

avgCD.Rhit.true = mean(mySmooth(trueVals.Rhit{sessix},81),2,'omitnan');      % Get average true CDlate for R and L hits for this session
avgCD.Lhit.true = mean(mySmooth(trueVals.Lhit{sessix},81),2,'omitnan');
stdCD.Rhit.true = std(mySmooth(trueVals.Rhit{sessix},81),0,2,'omitnan');     % Get standard deviation of true CDlate for R and L hits
stdCD.Lhit.true = std(mySmooth(trueVals.Lhit{sessix},81),0,2,'omitnan');

modelpred.Rhit{sessix} = fillmissing(modelpred.Rhit{sessix},"nearest");
modelpred.Lhit{sessix} = fillmissing(modelpred.Lhit{sessix},"nearest");
infix = find(isinf(modelpred.Rhit{sessix})); modelpred.Rhit{sessix}(infix) = 0;
infix = find(isinf(modelpred.Lhit{sessix})); modelpred.Lhit{sessix}(infix) = 0;

temp = [modelpred.Rhit{sessix},modelpred.Lhit{sessix}];
avgCD.allhits.pred = mean(mySmooth(temp,81),2,'omitnan');
stdCD.allhits.pred = mean(mySmooth(temp,81),2,'omitnan');

avgCD.Rhit.pred = mean(mySmooth(modelpred.Rhit{sessix},81),2,'omitnan');     % Get average predicted CDlate for R and L hits for this session
avgCD.Lhit.pred = mean(mySmooth(modelpred.Lhit{sessix},81),2,'omitnan');
stdCD.Rhit.pred = std(mySmooth(modelpred.Rhit{sessix},81),0,2,'omitnan');    % Get stdev of predicted CDlate for R and L hits for this session
stdCD.Lhit.pred = std(mySmooth(modelpred.Lhit{sessix},81),0,2,'omitnan');
end