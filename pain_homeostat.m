function S = pain_homeostat(wakeDurs,sleepDurs,wakeOnset,sleepOnset,t, swu, sso)

%time constants
ti = 15.78*60;
td = 3.37*60;

%asymptotes
UA = 0.4125;
LA = -0.5088;

%functions
swake = @(t, te, swu) (swu - UA)*exp(-(t - te)/ti) + UA;
ssleep = @(t, te, sso) (sso - LA)*exp(-(t - te)/td) + LA;

last_index = 2;
te = 0;
cuml = 0;

% if first wake time is larger than first sleep time, then it started in wake 
if wakeOnset(1) > sleepOnset(1)
    S = swu;
    if length(wakeDurs) == length(sleepDurs)
        d = length(wakeDurs);
    else
        d = length(wakeDurs) - 1;
    end
    for k = 1:d
        % wake exponential
        t_index = find(t>=(wakeDurs(k)+cuml),1); 
        t_period = t(last_index:t_index);
        func_eval = swake(t_period,te,swu); 
        S = [S; func_eval]; 
        sso = func_eval(end); 
        cuml = cuml + wakeDurs(k); 
        last_index = t_index + 1;
        te = t(last_index);

        % sleep exponential
        t_index = find(t>=(sleepDurs(k)+cuml),1);
        t_period = t(last_index:t_index);
        func_eval = ssleep(t_period,te,sso); 
        S = [S; func_eval]; 
        swu = func_eval(end); 
        if t(t_index) < t(end)
            cuml = cuml + sleepDurs(k);
            last_index = t_index + 1;
            te = t(last_index);
        end
    end
    % if it ends in wake, but started in wake, it will need extra wake exponential but not sleep
   if d<max(length(wakeDurs),length(sleepDurs))
        t_index = find(abs(t-(wakeDurs(end)+cuml))<0.001,1);
        t_period = t(last_index:t_index);
        func_eval = swake(t_period,te,swu); 
        S = [S; func_eval];
   end
else % it started in sleep
     S = sso;
    if length(wakeDurs) == length(sleepDurs)
        d = length(wakeDurs);
    else 
        d = length(sleepDurs) - 1;
    end
    for k = 1:d
        % sleep exponential
        t_index = find(t>=(sleepDurs(k)+cuml),1);
        t_period = t(last_index:t_index);
        func_eval = ssleep(t_period,te,sso); 
        S = [S; func_eval]; 
        swu = func_eval(end); 
        cuml = cuml + wakeDurs(k); 
        last_index = t_index + 1;
        te = t(last_index);

         % wake exponential
        t_index = find(t>=(wakeDurs(k)+cuml),1); 
        t_period = t(last_index:t_index);
        func_eval = swake(t_period,te,swu); 
        S = [S; func_eval]; 
        sso = func_eval(end); 
        if t(t_index) < t(end)
            cuml = cuml + wakeDurs(k);
            last_index = t_index + 1;
            te = t(last_index);
        end
    end
    % if it ends in sleep and started in sleep, it will need sleep wake exponential but not sleep
   if d<max(length(wakeDurs),length(sleepDurs))
        t_index = find(abs(t-(sleepDurs(end)+cuml))<0.001,1);
        t_period = t(last_index:t_index);
        func_eval = ssleep(t_period,te,sso); 
        S = [S; func_eval];
   end

end

end
