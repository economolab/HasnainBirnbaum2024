function tonguefeat = getRawTongueKin(meta,feat2use,cond2use,params,times,nLicks,lickDur,kin,featix,condfns)
for sessix = 1:length(meta)                                                 % For every session...
    for feat = 1:length(feat2use)                                           % For tongue length and angle...
        metric = feat2use{feat};
        for c = 1:length(cond2use)                                          % For each condition...
            cond = cond2use(c);
            condtrix = params(sessix).trialid{cond};                        % Get the trials for this condition
            condfeat = squeeze(kin(sessix).dat(:,condtrix,featix(feat)));   % Get the kinematic data for these trials and current feature
            respfeat = condfeat(times.goix:times.respix,:);                             % Take only desired post-goCue time points

            maxDur = 20;
            nTrials = size(respfeat,2);                                     % number of trials in this condition
            nTimePts = size(respfeat,1);                                    % number of time-points
            temp = cell(nLicks,nTrials);                                    % Pre-allocate (nLicks x nTrials) for this condition
            for trix = 1:nTrials                                            % For each trial in this condition...
                licks.start = [];
                licks.stop = [];
                % Find indices corresponding to lick start and lick end
                for pt = 4:nTimePts-4                                       % For every time point...
                    % If the tongue is visible for this time point and the subsequent two time points (but not the previous time
                    % point, take this point to be the start of a lick)
                    check1 = ~isnan(respfeat(pt,trix));                     
                    check2 = isnan(respfeat(pt-1,trix));
                    check3 = ~isnan(respfeat(pt+1,trix));
                    check4 = ~isnan(respfeat(pt+2,trix));
                    allcheck = check1+check2+check3+check4;
                    if allcheck==4
                        licks.start = [licks.start,pt];
                    end

                    % If the tongue is visible for this time point but not the subsequent two time points (but is visible the previous time
                    % point, take this point to be the end of a lick)
                    check1 = ~isnan(respfeat(pt,trix));                     
                    check2 = ~isnan(respfeat(pt-1,trix));
                    check3 = isnan(respfeat(pt+1,trix));
                    check4 = isnan(respfeat(pt+2,trix));
                    allcheck = check1+check2+check3+check4;
                    if allcheck==4
                        licks.stop = [licks.stop,pt];
                    end
                end

                for lix = 1:nLicks                                                  % For each lick that you want to look at
                    if lix<size(licks.start,2)&&lix<size(licks.stop,2)              % If there are this many licks in the current trial...
                        currstart = licks.start(lix);
                        currstop = licks.stop(lix);
                        if currstart<currstop
                            currtongue = respfeat(currstart:currstop,trix);             % Take the tongue length/angle for this current lick and store it
                            currdur = size(currtongue,1);                               % Get the duration in frames of each individual lick
%                             avgDur(lix,trix) = currdur;
%                             tongueinterp = interp1(1:length(currtongue),currtongue, linspace(1,length(currtongue),lickDur));
                            temp{lix,trix} = currtongue;
                        end
                    end
                end
            end
            tonguefeat(sessix).(metric).(condfns{c}) = temp;
        end
    end
end