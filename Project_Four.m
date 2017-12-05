function Prject_Four

%Question 1
fr = 100;
tSim = 1000; %in ms
dt = 1; %in ms
nTrials = 100; % number of simulations
[spikeMat, tVec, r, spikes] = poissonSpikeGen(fr, tSim, nTrials, dt);

figure();
plot(r);
title('Average Spiking Rates');
xlabel('Time (in ms )');
ylabel('Spike rate ( in #spikes per second )');
xlim([0,tSim/dt]);
ylim([40,50]);

figure();
plot(spikes);
title('Spike Train');
xlabel('Time ( in ms )');
ylabel('Spike');
xlim([0,tSim/dt]);
ylim([0,1]);

spikeMat = double(spikeMat);
spikecount = sum(spikeMat);
PSTH = spikecount/nTrials*1000;


figure();
plot(PSTH);
title('PSTH');
xlabel('Time ( in ms )');
ylabel('Spike rate ( in #spikes per second )');
xlim([0,tSim/dt]);

trialcount = [10, 20, 40, 80, 160, 320, 640, 1280];
l = size(trialcount);
m = l(2);

for j = 1:m
    nTrials = trialcount(j);
    [spikeMat, tVec, r, spikes] = poissonSpikeGen(fr, tSim, nTrials, dt);
    spikeMat = double(spikeMat);
    spikecount = sum(spikeMat);
    PST = spikecount/nTrials*1000;
    error = PST - r;
 
    rmse(j) = rms(error);
end

figure(4);
plot(trialcount,rmse);
title('RMSE vs # of trials');
xlabel('# of trials');
ylabel('RMSE');


%Stimulus S and D responses
t=1;
for i = 1:15
    for j = 1:50
        if i ~= 8
            rS_ODD(t) = 10;
            rD_ODD(t) = 2.5;
        else
            rS_ODD(t) = 2.5;
            rD_ODD(t) = 10;
        end
        t = t + 1;
    end
    for j = 1:250
        rS_ODD(t) = 0.5;
        rD_ODD(t) = 0.5;
        t = t + 1;
    end
end

t = t - 1;

P_S_ODD = rS_ODD.*0.001;
P_D_ODD = rD_ODD.*0.001;

for i = 1:t;
    SpikeS_ODD(i) = rand < P_S_ODD(i);
end

for i = 1:t;
    SpikeD_ODD(i) = rand < P_D_ODD(i);
end

for  i = 1:t;
    if SpikeS_ODD(i) == 1
        M_S_ODD(i) = 1;
        i = i + 1;
        M_S_ODD(i) = 1;
    else
        M_S_ODD(i) = 0;
    end
end

for  i = 1:t
    if SpikeD_ODD(i) == 1
        M_D_ODD(i) = 1;
        i = i + 1;
        M_D_ODD(i) = 1;
    else
        M_D_ODD(i) = 0;
    end
end

for loop = 1:1;
    %initializing the values
    timesteps = t;
    tsyn = 10;
    B = 5;
    tref = 2;
    ws = 0.2.*ones(1,timesteps);
    wd = 0.2.*ones(1,timesteps);
    wsp = 0.11.*ones(1,timesteps);
    wsl4 = 0.02.*ones(1,timesteps);
    wdl4 = 0.02.*ones(1,timesteps); %all time constants are in milli seconds
    S.xr = 1;
    S.xi = 0;
    S.xe = 0;
    S.tre = 0.9;
    S.tei = 10;
    S.tir = 5000;
    S.M = M_S_ODD;
    
    
    S.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            S.g(i) = S.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                S.g(i) = S.g(i-1)*exp(-1/tsyn);
            end
        end
        if S.M(i) == 1
            flag = 1;
        end
    end
    D.xr = 1;
    D.xi = 0;
    D.xe = 0;
    D.tre = 0.9;
    D.tei = 10;
    D.tir = 5000;
    D.M = M_D_ODD;
    D.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            D.g(i) = D.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                D.g(i) = D.g(i-1)*exp(-1/tsyn);
            end
        end     
        if D.M(i)==1
            flag=1;
        end
    end
    SP.xr = 1;
    SP.xi = 0;
    SP.xe = 0;
    SP.tre = 0.9;
    SP.tei = 27;
    SP.tir = 5000;
    SP.M = zeros(1,timesteps);
    SP.g = zeros(1,timesteps);
    countervsp = 0;
    countervl4 = 0;
    Vsp = zeros(1,timesteps+1);
    Vl4 = zeros(1,timesteps);
    Vl4M = zeros(1,timesteps);
    flag=0;%For SP.g update
    
    tm= linspace(1,4500,4500);
    options =odeset('MaxStep',5);
    [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    figure();
    plot(T_S, S.x(:,1),'r-');
    S.xr = interp1(T_S,S.x(:,1),tm,'spline');
    title('S.xr');
    
    figure();
    plot(T_S, S.x(:,2),'r-');
    S.xe = interp1(T_S,S.x(:,2),tm,'spline');
    title('S.xe');
    
    figure();
    plot(T_S,S.x(:,3),'r-');
    S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    title('S.xi');
    
    [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
    figure();
    plot(T_S,D.x(:,1),'r-');
    D.xr = interp1(T_S,D.x(:,1),tm,'spline');
    title('D.xr');
    
    figure();
    plot(T_S,D.x(:,2),'r-');
    D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    title('D.xe');
    
    figure();
    plot(T_S,D.x(:,3),'r-');
    D.xi = interp1(T_S,D.x(:,3),tm,'spline');
    title('D.xi');
    
    for t = 1:timesteps
        if flag == 1
            SP.g(t) = SP.g(t-1) + 1;
            flag = 0;
        else
            if t ~= 1
                SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
            end
        end
        if SP.M(t) == 1
            flag = 1;
        end
        %Calculate the changes and assign them to a temporary variable
        if countervsp == 0
            Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
        else
            if SP.M(t) == 1
                Vsp(t+1) = -B;
            else
                Vsp(t+1) = Vsp(t)*exp(-1/tref);
            end
        countervsp = countervsp - 1;
        end
        
        if Vsp(t+1) > 0.05
            SP.M(t+1) = 1;
            countervsp = 20;
        end
    end
    [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    figure(11);
    plot(T_S, SP.x(:,1),'r-');
    SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');
    title('SP.xr');
    
    figure(12);
    plot(T_S, SP.x(:,2),'r-');
    SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');
    title('SP.xe');
    
    figure(13);
    plot(T_S,SP.x(:,3),'r-');
    SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
    title('SP.xi');
    
    for t = 1:timesteps
        if countervl4==0
            Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
        else
            if Vl4M(t)==1
                Vl4(t+1)=-B;
            else
                Vl4(t+1)=Vl4(t)*exp(-1/tref);
            end
            countervl4=countervl4-1;
        end
        if Vl4(t+1)>0.05
            Vl4M(t+1)=1;
            countervl4=20;
        end
    end
    
end

figure();
subplot(4,1,1);
plot(M_S_ODD)
title('Stimulus to S');
xlabel('Time ( in ms )');
ylabel('Spike at S');

subplot(4,1,2);
plot(M_D_ODD);
title('Stimulus to D');
xlabel('Time ( in ms )');
ylabel('Spike at D');

subplot(4,1,3);
plot(Vsp);
title('SP neuron');
xlabel('Time ( in ms )');
ylabel('Vsp');

subplot(4,1,4);
plot(Vl4);
title('L4 neuron');
xlabel('Time ( in ms )');
ylabel('Vl4');

















%Question 2

%Stimulus S and D responses
t=1;
for i = 1:15
    for j = 1:50
        if i ~= 8
            rS_ODD(t) = 10;
            rD_ODD(t) = 2.5;
        else
            rS_ODD(t) = 2.5;
            rD_ODD(t) = 10;
        end
        t = t + 1;
    end
    for j = 1:250
        rS_ODD(t) = 0.5;
        rD_ODD(t) = 0.5;
        t = t + 1;
    end
end

t = t - 1;

P_S_ODD = rS_ODD.*0.001;
P_D_ODD = rD_ODD.*0.001;


SPM = [];
L4M = [];
for loop = 1:50;
    for i = 1:t;    
        SpikeS_ODD(i) = rand < P_S_ODD(i);
    end
    for i = 1:t;    
        SpikeD_ODD(i) = rand < P_D_ODD(i);
    end   
    for  i = 1:t;    
        if SpikeS_ODD(i) == 1
        
            M_S_ODD(i) = 1;        
            i = i + 1;        
            M_S_ODD(i) = 1;  
        else            
            M_S_ODD(i) = 0;    
        end        
    end       
    for  i = 1:t    
        if SpikeD_ODD(i) == 1        
            M_D_ODD(i) = 1;
            i = i + 1;
            M_D_ODD(i) = 1;    
        else            
            M_D_ODD(i) = 0;    
        end        
    end
    %initializing the values
    timesteps = t;
    tsyn = 10;
    B = 5;
    tref = 2;
    ws = 0.2.*ones(1,timesteps);
    wd = 0.2.*ones(1,timesteps);
    wsp = 0.11.*ones(1,timesteps);
    wsl4 = 0.02.*ones(1,timesteps);
    wdl4 = 0.02.*ones(1,timesteps); %all time constants are in milli seconds
    S.xr = 1;
    S.xi = 0;
    S.xe = 0;
    S.tre = 0.9;
    S.tei = 10;
    S.tir = 5000;
    S.M = M_S_ODD;
    
    
    S.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            S.g(i) = S.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                S.g(i) = S.g(i-1)*exp(-1/tsyn);
            end
        end
        if S.M(i) == 1
            flag = 1;
        end
    end
    D.xr = 1;
    D.xi = 0;
    D.xe = 0;
    D.tre = 0.9;
    D.tei = 10;
    D.tir = 5000;
    D.M = M_D_ODD;
    D.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            D.g(i) = D.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                D.g(i) = D.g(i-1)*exp(-1/tsyn);
            end
        end     
        if D.M(i)==1
            flag=1;
        end
    end
    SP.xr = 1;
    SP.xi = 0;
    SP.xe = 0;
    SP.tre = 0.9;
    SP.tei = 27;
    SP.tir = 5000;
    SP.M = zeros(1,timesteps);
    SP.g = zeros(1,timesteps);
    countervsp = 0;
    countervl4 = 0;
    Vsp = zeros(1,timesteps+1);
    Vl4 = zeros(1,timesteps);
    Vl4M = zeros(1,timesteps);
    flag=0;%For SP.g update
    
    tm= linspace(1,4500,4500);
    options =odeset('MaxStep',5);
    [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);   
    S.xr = interp1(T_S,S.x(:,1),tm,'spline');

    S.xe = interp1(T_S,S.x(:,2),tm,'spline');

    S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    
    [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
    D.xr = interp1(T_S,D.x(:,1),tm,'spline');
    
    D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    
    D.xi = interp1(T_S,D.x(:,3),tm,'spline');

    
    for t = 1:timesteps
        if flag == 1
            SP.g(t) = SP.g(t-1) + 1;
            flag = 0;
        else
            if t ~= 1
                SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
            end
        end
        if SP.M(t) == 1
            flag = 1;
        end
        %Calculate the changes and assign them to a temporary variable
        if countervsp == 0
            Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
        else
            if SP.M(t) == 1
                Vsp(t+1) = -B;
            else
                Vsp(t+1) = Vsp(t)*exp(-1/tref);
            end
        countervsp = countervsp - 1;
        end
        
        if Vsp(t+1) > 0.05
            SP.M(t+1) = 1;
            countervsp = 20;
        end
    end
    [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');


    SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');

    SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
    
    for t = 1:timesteps
        if countervl4==0
            Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
        else
            if Vl4M(t)==1
                Vl4(t+1)=-B;
            else
                Vl4(t+1)=Vl4(t)*exp(-1/tref);
            end
            countervl4=countervl4-1;
        end
        if Vl4(t+1)>0.05
            Vl4M(t+1)=1;
            countervl4=20;
        end
    end
    SPM = [SPM; SP.M];
    L4M = [L4M;Vl4M];
    
end

PSTHsp=sum(SPM)./(loop);
PSTHl4=sum(L4M)./(loop);
psthsp=[];
psthl4=[];
time=[];
length = t;
for i=1:length/10
    sum1=0;
    sum2=0;
    for j=1:10
        sum1=sum1+PSTHsp((i-1)*10+j);
        sum2=sum2+PSTHl4((i-1)*10+j);
    end
    psthsp=[psthsp,sum1];
    psthl4=[psthl4,sum2];
    time(i)=(i-1)*10;
end
psthsp=psthsp.*100;
psthl4=psthl4.*100;

figure();
subplot(2,1,1);
plot(psthsp);
title('PSTH SP neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 5000ms )');

subplot(2,1,2);
plot(psthl4);
title('PSTH L4 neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 5000ms )');









%Question 3.a

%Stimulus S and D responses
t=1;
for i = 1:15
    for j = 1:50
        if i ~= 8
            rS_ODD(t) = 10;
            rD_ODD(t) = 2.5;
        else
            rS_ODD(t) = 2.5;
            rD_ODD(t) = 10;
        end
        t = t + 1;
    end
    for j = 1:250
        rS_ODD(t) = 0.5;
        rD_ODD(t) = 0.5;
        t = t + 1;
    end
end

t = t - 1;

P_S_ODD = rS_ODD.*0.001;
P_D_ODD = rD_ODD.*0.001;


SPM = [];
L4M = [];
for loop = 1:50;
    for i = 1:t;    
        SpikeS_ODD(i) = rand < P_S_ODD(i);
    end
    for i = 1:t;    
        SpikeD_ODD(i) = rand < P_D_ODD(i);
    end   
    for  i = 1:t;    
        if SpikeS_ODD(i) == 1
        
            M_S_ODD(i) = 1;        
            i = i + 1;        
            M_S_ODD(i) = 1;  
        else            
            M_S_ODD(i) = 0;    
        end        
    end       
    for  i = 1:t    
        if SpikeD_ODD(i) == 1        
            M_D_ODD(i) = 1;
            i = i + 1;
            M_D_ODD(i) = 1;    
        else            
            M_D_ODD(i) = 0;    
        end        
    end
    %initializing the values
    timesteps = t;
    tsyn = 10;
    B = 5;
    tref = 2;
    ws = 0.2.*ones(1,timesteps);
    wd = 0.2.*ones(1,timesteps);
    wsp = 0.11.*ones(1,timesteps);
    wsl4 = 0.02.*ones(1,timesteps);
    wdl4 = 0.02.*ones(1,timesteps); %all time constants are in milli seconds
    S.xr = 1;
    S.xi = 0;
    S.xe = 0;
    S.tre = 0.9;
    S.tei = 10;
    S.tir = 1000;
    S.M = M_S_ODD;
    
    
    S.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            S.g(i) = S.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                S.g(i) = S.g(i-1)*exp(-1/tsyn);
            end
        end
        if S.M(i) == 1
            flag = 1;
        end
    end
    D.xr = 1;
    D.xi = 0;
    D.xe = 0;
    D.tre = 0.9;
    D.tei = 10;
    D.tir = 1000;
    D.M = M_D_ODD;
    D.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            D.g(i) = D.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                D.g(i) = D.g(i-1)*exp(-1/tsyn);
            end
        end     
        if D.M(i)==1
            flag=1;
        end
    end
    SP.xr = 1;
    SP.xi = 0;
    SP.xe = 0;
    SP.tre = 0.9;
    SP.tei = 27;
    SP.tir = 1000;
    SP.M = zeros(1,timesteps);
    SP.g = zeros(1,timesteps);
    countervsp = 0;
    countervl4 = 0;
    Vsp = zeros(1,timesteps+1);
    Vl4 = zeros(1,timesteps);
    Vl4M = zeros(1,timesteps);
    flag=0;%For SP.g update
    
    tm= linspace(1,4500,4500);
    options =odeset('MaxStep',5);
    [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);   
    S.xr = interp1(T_S,S.x(:,1),tm,'spline');

    S.xe = interp1(T_S,S.x(:,2),tm,'spline');

    S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    
    [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
    D.xr = interp1(T_S,D.x(:,1),tm,'spline');
    
    D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    
    D.xi = interp1(T_S,D.x(:,3),tm,'spline');

    
    for t = 1:timesteps
        if flag == 1
            SP.g(t) = SP.g(t-1) + 1;
            flag = 0;
        else
            if t ~= 1
                SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
            end
        end
        if SP.M(t) == 1
            flag = 1;
        end
        %Calculate the changes and assign them to a temporary variable
        if countervsp == 0
            Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
        else
            if SP.M(t) == 1
                Vsp(t+1) = -B;
            else
                Vsp(t+1) = Vsp(t)*exp(-1/tref);
            end
        countervsp = countervsp - 1;
        end
        
        if Vsp(t+1) > 0.05
            SP.M(t+1) = 1;
            countervsp = 20;
        end
    end
    [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');


    SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');

    SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
    
    for t = 1:timesteps
        if countervl4==0
            Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
        else
            if Vl4M(t)==1
                Vl4(t+1)=-B;
            else
                Vl4(t+1)=Vl4(t)*exp(-1/tref);
            end
            countervl4=countervl4-1;
        end
        if Vl4(t+1)>0.05
            Vl4M(t+1)=1;
            countervl4=20;
        end
    end
    SPM = [SPM; SP.M];
    L4M = [L4M;Vl4M];
    
end

PSTHsp=sum(SPM)./(loop);
PSTHl4=sum(L4M)./(loop);
psthsp=[];
psthl4=[];
time=[];
length = t;
for i=1:length/10
    sum1=0;
    sum2=0;
    for j=1:10
        sum1=sum1+PSTHsp((i-1)*10+j);
        sum2=sum2+PSTHl4((i-1)*10+j);
    end
    psthsp=[psthsp,sum1];
    psthl4=[psthl4,sum2];
    time(i)=(i-1)*10;
end
psthsp=psthsp.*100;
psthl4=psthl4.*100;

figure();
subplot(2,1,1);
plot(psthsp);
title('PSTH SP neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 1000ms )');

subplot(2,1,2);
plot(psthl4);
title('PSTH L4 neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 1000ms )');














%Question 3.b

%Stimulus S and D responses
t=1;
for i = 1:15
    for j = 1:50
        if i ~= 8
            rS_ODD(t) = 10;
            rD_ODD(t) = 2.5;
        else
            rS_ODD(t) = 2.5;
            rD_ODD(t) = 10;
        end
        t = t + 1;
    end
    for j = 1:250
        rS_ODD(t) = 0.5;
        rD_ODD(t) = 0.5;
        t = t + 1;
    end
end

t = t - 1;

P_S_ODD = rS_ODD.*0.001;
P_D_ODD = rD_ODD.*0.001;


SPM = [];
L4M = [];
for loop = 1:50;
    for i = 1:t;    
        SpikeS_ODD(i) = rand < P_S_ODD(i);
    end
    for i = 1:t;    
        SpikeD_ODD(i) = rand < P_D_ODD(i);
    end   
    for  i = 1:t;    
        if SpikeS_ODD(i) == 1
        
            M_S_ODD(i) = 1;        
            i = i + 1;        
            M_S_ODD(i) = 1;  
        else            
            M_S_ODD(i) = 0;    
        end        
    end       
    for  i = 1:t    
        if SpikeD_ODD(i) == 1        
            M_D_ODD(i) = 1;
            i = i + 1;
            M_D_ODD(i) = 1;    
        else            
            M_D_ODD(i) = 0;    
        end        
    end
    %initializing the values
    timesteps = t;
    tsyn = 10;
    B = 5;
    tref = 2;
    ws = 0.2.*ones(1,timesteps);
    wd = 0.2.*ones(1,timesteps);
    wsp = 0.11.*ones(1,timesteps);
    wsl4 = 0.02.*ones(1,timesteps);
    wdl4 = 0.02.*ones(1,timesteps); %all time constants are in milli seconds
    S.xr = 1;
    S.xi = 0;
    S.xe = 0;
    S.tre = 0.9;
    S.tei = 10;
    S.tir = 3000;
    S.M = M_S_ODD;
    
    
    S.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            S.g(i) = S.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                S.g(i) = S.g(i-1)*exp(-1/tsyn);
            end
        end
        if S.M(i) == 1
            flag = 1;
        end
    end
    D.xr = 1;
    D.xi = 0;
    D.xe = 0;
    D.tre = 0.9;
    D.tei = 10;
    D.tir = 3000;
    D.M = M_D_ODD;
    D.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            D.g(i) = D.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                D.g(i) = D.g(i-1)*exp(-1/tsyn);
            end
        end     
        if D.M(i)==1
            flag=1;
        end
    end
    SP.xr = 1;
    SP.xi = 0;
    SP.xe = 0;
    SP.tre = 0.9;
    SP.tei = 27;
    SP.tir = 3000;
    SP.M = zeros(1,timesteps);
    SP.g = zeros(1,timesteps);
    countervsp = 0;
    countervl4 = 0;
    Vsp = zeros(1,timesteps+1);
    Vl4 = zeros(1,timesteps);
    Vl4M = zeros(1,timesteps);
    flag=0;%For SP.g update
    
    tm= linspace(1,4500,4500);
    options =odeset('MaxStep',5);
    [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);   
    S.xr = interp1(T_S,S.x(:,1),tm,'spline');

    S.xe = interp1(T_S,S.x(:,2),tm,'spline');

    S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    
    [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
    D.xr = interp1(T_S,D.x(:,1),tm,'spline');
    
    D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    
    D.xi = interp1(T_S,D.x(:,3),tm,'spline');

    
    for t = 1:timesteps
        if flag == 1
            SP.g(t) = SP.g(t-1) + 1;
            flag = 0;
        else
            if t ~= 1
                SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
            end
        end
        if SP.M(t) == 1
            flag = 1;
        end
        %Calculate the changes and assign them to a temporary variable
        if countervsp == 0
            Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
        else
            if SP.M(t) == 1
                Vsp(t+1) = -B;
            else
                Vsp(t+1) = Vsp(t)*exp(-1/tref);
            end
        countervsp = countervsp - 1;
        end
        
        if Vsp(t+1) > 0.05
            SP.M(t+1) = 1;
            countervsp = 20;
        end
    end
    [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');


    SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');

    SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
    
    for t = 1:timesteps
        if countervl4==0
            Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
        else
            if Vl4M(t)==1
                Vl4(t+1)=-B;
            else
                Vl4(t+1)=Vl4(t)*exp(-1/tref);
            end
            countervl4=countervl4-1;
        end
        if Vl4(t+1)>0.05
            Vl4M(t+1)=1;
            countervl4=20;
        end
    end
    SPM = [SPM; SP.M];
    L4M = [L4M;Vl4M];
    
end

PSTHsp=sum(SPM)./(loop);
PSTHl4=sum(L4M)./(loop);
psthsp=[];
psthl4=[];
time=[];
length = t;
for i=1:length/10
    sum1=0;
    sum2=0;
    for j=1:10
        sum1=sum1+PSTHsp((i-1)*10+j);
        sum2=sum2+PSTHl4((i-1)*10+j);
    end
    psthsp=[psthsp,sum1];
    psthl4=[psthl4,sum2];
    time(i)=(i-1)*10;
end
psthsp=psthsp.*100;
psthl4=psthl4.*100;

figure();
subplot(2,1,1);
plot(psthsp);
title('PSTH SP neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 3000ms )');

subplot(2,1,2);
plot(psthl4);
title('PSTH L4 neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 3000ms)');









%Question 3.c

%Stimulus S and D responses
t=1;
for i = 1:15
    for j = 1:50
        if i ~= 8
            rS_ODD(t) = 10;
            rD_ODD(t) = 2.5;
        else
            rS_ODD(t) = 2.5;
            rD_ODD(t) = 10;
        end
        t = t + 1;
    end
    for j = 1:250
        rS_ODD(t) = 0.5;
        rD_ODD(t) = 0.5;
        t = t + 1;
    end
end

t = t - 1;

P_S_ODD = rS_ODD.*0.001;
P_D_ODD = rD_ODD.*0.001;


SPM = [];
L4M = [];
for loop = 1:50;
    for i = 1:t;    
        SpikeS_ODD(i) = rand < P_S_ODD(i);
    end
    for i = 1:t;    
        SpikeD_ODD(i) = rand < P_D_ODD(i);
    end   
    for  i = 1:t;    
        if SpikeS_ODD(i) == 1
        
            M_S_ODD(i) = 1;        
            i = i + 1;        
            M_S_ODD(i) = 1;  
        else            
            M_S_ODD(i) = 0;    
        end        
    end       
    for  i = 1:t    
        if SpikeD_ODD(i) == 1        
            M_D_ODD(i) = 1;
            i = i + 1;
            M_D_ODD(i) = 1;    
        else            
            M_D_ODD(i) = 0;    
        end        
    end
    %initializing the values
    timesteps = t;
    tsyn = 10;
    B = 5;
    tref = 2;
    ws = 0.2.*ones(1,timesteps);
    wd = 0.2.*ones(1,timesteps);
    wsp = 0.11.*ones(1,timesteps);
    wsl4 = 0.02.*ones(1,timesteps);
    wdl4 = 0.02.*ones(1,timesteps); %all time constants are in milli seconds
    S.xr = 1;
    S.xi = 0;
    S.xe = 0;
    S.tre = 0.9;
    S.tei = 10;
    S.tir = 10000;
    S.M = M_S_ODD;
    
    
    S.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            S.g(i) = S.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                S.g(i) = S.g(i-1)*exp(-1/tsyn);
            end
        end
        if S.M(i) == 1
            flag = 1;
        end
    end
    D.xr = 1;
    D.xi = 0;
    D.xe = 0;
    D.tre = 0.9;
    D.tei = 10;
    D.tir = 10000;
    D.M = M_D_ODD;
    D.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            D.g(i) = D.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                D.g(i) = D.g(i-1)*exp(-1/tsyn);
            end
        end     
        if D.M(i)==1
            flag=1;
        end
    end
    SP.xr = 1;
    SP.xi = 0;
    SP.xe = 0;
    SP.tre = 0.9;
    SP.tei = 27;
    SP.tir = 10000;
    SP.M = zeros(1,timesteps);
    SP.g = zeros(1,timesteps);
    countervsp = 0;
    countervl4 = 0;
    Vsp = zeros(1,timesteps+1);
    Vl4 = zeros(1,timesteps);
    Vl4M = zeros(1,timesteps);
    flag=0;%For SP.g update
    
    tm= linspace(1,4500,4500);
    options =odeset('MaxStep',5);
    [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);   
    S.xr = interp1(T_S,S.x(:,1),tm,'spline');

    S.xe = interp1(T_S,S.x(:,2),tm,'spline');

    S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    
    [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
    D.xr = interp1(T_S,D.x(:,1),tm,'spline');
    
    D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    
    D.xi = interp1(T_S,D.x(:,3),tm,'spline');

    
    for t = 1:timesteps
        if flag == 1
            SP.g(t) = SP.g(t-1) + 1;
            flag = 0;
        else
            if t ~= 1
                SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
            end
        end
        if SP.M(t) == 1
            flag = 1;
        end
        %Calculate the changes and assign them to a temporary variable
        if countervsp == 0
            Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
        else
            if SP.M(t) == 1
                Vsp(t+1) = -B;
            else
                Vsp(t+1) = Vsp(t)*exp(-1/tref);
            end
        countervsp = countervsp - 1;
        end
        
        if Vsp(t+1) > 0.05
            SP.M(t+1) = 1;
            countervsp = 20;
        end
    end
    [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');


    SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');

    SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
    
    for t = 1:timesteps
        if countervl4==0
            Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
        else
            if Vl4M(t)==1
                Vl4(t+1)=-B;
            else
                Vl4(t+1)=Vl4(t)*exp(-1/tref);
            end
            countervl4=countervl4-1;
        end
        if Vl4(t+1)>0.05
            Vl4M(t+1)=1;
            countervl4=20;
        end
    end
    SPM = [SPM; SP.M];
    L4M = [L4M;Vl4M];
    
end

PSTHsp=sum(SPM)./(loop);
PSTHl4=sum(L4M)./(loop);
psthsp=[];
psthl4=[];
time=[];
length = t;
for i=1:length/10
    sum1=0;
    sum2=0;
    for j=1:10
        sum1=sum1+PSTHsp((i-1)*10+j);
        sum2=sum2+PSTHl4((i-1)*10+j);
    end
    psthsp=[psthsp,sum1];
    psthl4=[psthl4,sum2];
    time(i)=(i-1)*10;
end
psthsp=psthsp.*100;
psthl4=psthl4.*100;

figure(15);
subplot(2,1,1);
plot(psthsp);
title('PSTH SP neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 10000ms)');

subplot(2,1,2);
plot(psthl4);
title('PSTH L4 neuron');
xlabel('Time/10 ( in ms )');
ylabel('Spike rate ( in #spikes per second at Tre = 10000ms)');





















%Question 4
%%Stimulus S and D responses
t=1;
diffstim=1200;
seq=[];
for i=1:diffstim
    base=rand();
    if(base<0.9)
        seq=[seq,1];%S=1, D=0
    else
        seq=[seq,0];
    end
end


for i=1:diffstim
    for j=1:50
        if seq(i)==1
            rS(t)=10;
            rD(t)=2.5;
        else
            rS(t)=2.5;
            rD(t)=10;
        end
        t=t+1;
    end
    for j=1:250
        rS(t)=0.5;
        rD(t)=0.5;
        t=t+1;
    end
end

n=t-1;
t = t- 1;

P_S = rS.*0.001;
P_D = rD.*0.001;

for i = 1:t;
    SpikeS(i) = rand < P_S(i);
end

for i = 1:t;
    SpikeD(i) = rand < P_D(i);
end


for  i = 1:t;
    if SpikeS(i) == 1
        M_S(i) = 1;
        i = i + 1;
        M_S(i) = 1;
    else
        M_S(i) = 0;
    end
end

for  i = 1:t;
    if SpikeD(i) == 1
        M_D(i) = 1;
        i = i + 1;
        M_D(i) = 1;
    else
        M_D(i) = 0;
    end
end

for loop = 1:1;
    %initializing the values
    timesteps = t;
    tsyn = 10;
    B = 5;
    tref = 2;
    ws = 0.2.*ones(1,timesteps);
    wd = 0.2.*ones(1,timesteps);
    wsp = 0.11.*ones(1,timesteps);
    wsl4 = 0.02.*ones(1,timesteps);
    wdl4 = 0.02.*ones(1,timesteps); %all time constants are in milli seconds
    tltp=13;
    tltd=20;
    altp=0.015;
    altd=0.021;
    postl4sp = -1;%post spike time for sp on L4 neuron
    postl4s = -1;%post spike time for s on L4 neuron
    postl4d = -1;%post spike time for d on L4 neuron
    pres = -1;%pre spike time of s on l4
    pred = -1;%pre spike time of d on l4
    presp =- 1;%pre spike time of sp on l4

    S.xr = 1;
    S.xi = 0;
    S.xe = 0;
    S.tre = 0.9;
    S.tei = 10;
    S.tir = 500;
    S.M = M_S;
    
    
    S.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            S.g(i) = S.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                S.g(i) = S.g(i-1)*exp(-1/tsyn);
            end
        end
        if S.M(i) == 1
            flag = 1;
        end
    end
    D.xr = 1;
    D.xi = 0;
    D.xe = 0;
    D.tre = 0.9;
    D.tei = 10;
    D.tir = 500;
    D.M = M_D;
    D.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            D.g(i) = D.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                D.g(i) = D.g(i-1)*exp(-1/tsyn);
            end
        end     
        if D.M(i)==1
            flag=1;
        end
    end
    SP.xr = 1;
    SP.xi = 0;
    SP.xe = 0;
    SP.tre = 0.9;
    SP.tei = 27;
    SP.tir = 500;
    SP.M = zeros(1,timesteps);
    SP.g = zeros(1,timesteps);
    countervsp = 0;
    countervl4 = 0;
    Vsp = zeros(1,timesteps+1);
    Vl4 = zeros(1,timesteps);
    Vl4M = zeros(1,timesteps);
    flag=0;%For SP.g update
    
    tm= linspace(1,t,t);
    options =odeset('MaxStep',5);
    [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);   
    S.xr = interp1(T_S,S.x(:,1),tm,'spline');

    S.xe = interp1(T_S,S.x(:,2),tm,'spline');

    S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    
    [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
    D.xr = interp1(T_S,D.x(:,1),tm,'spline');
   
    D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    
    D.xi = interp1(T_S,D.x(:,3),tm,'spline');
    
    
    for t = 1:timesteps
        if flag == 1
            SP.g(t) = SP.g(t-1) + 1;
            flag = 0;
        else
            if t ~= 1
                SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
            end
        end
        if SP.M(t) == 1
            flag = 1;
        end
        %Calculate the changes and assign them to a temporary variable
        if countervsp == 0
            Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
        else
            if SP.M(t) == 1
                Vsp(t+1) = -B;
            else
                Vsp(t+1) = Vsp(t)*exp(-1/tref);
            end
        countervsp = countervsp - 1;
        end
        
        if Vsp(t+1) > 0.05
            SP.M(t+1) = 1;
            countervsp = 20;
        end
    end
    [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');


    SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');

    SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
    
    for t = 1:timesteps
        if countervl4==0
            Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
        else
            if Vl4M(t)==1
                Vl4(t+1)=-B;
            else
                Vl4(t+1)=Vl4(t)*exp(-1/tref);
            end
            countervl4=countervl4-1;
        end
        if Vl4(t+1)>0.05
            Vl4M(t+1)=1;
            countervl4=20;
        end
    end
    for t=1:timesteps        %Updating weights for LTP/LTD
        if(Vl4M(t)==1)
            postl4sp = t;
            postl4s = t;
            postl4d = t;
        end
        if(SP.M(t)==1)
            presp=t;
        end
        if(S.M(t)==1)
            pres=t;
        end
        if(D.M(t)==1)
            pred=t;
        end
        if(t~=1)
            wsp(t) = wsp(t-1);
            wsl4(t)= wsl4(t-1);
            wdl4(t)= wdl4(t-1);
        end
        if(postl4sp~=-1)
            if(presp~=-1&&(Vl4M(t)==1||SP.M(t)==1))
                dt=postl4sp-presp;
                if(dt>0)
                    wsp(t)=wsp(t)*(1+altp*exp(-dt/tltp));
                else
                    wsp(t)=wsp(t)*(1-altd*exp(dt/tltd));
                end
                %postl4sp=-1;
                %presp=-1;
            end
        end      
        if(postl4s~=-1)
            if(pres~=-1&&(Vl4M(t)==1||S.M(t)==1))
                dt=postl4s-pres;
                if(dt>0)
                    wsl4(t)=wsl4(t)*(1+altp*exp(-dt/tltp));
                else
                    wsl4(t)=wsl4(t)*(1-altd*exp(dt/tltd));
                end
                %postl4s=-1;
                %pres=-1;
            end
        end
        if(postl4d~=-1)
            if(pred~=-1&&(Vl4M(t)==1||D.M(t)==1))
                dt=postl4d-pred;
                if(dt>0)
                    wdl4(t)=wdl4(t)*(1+altp*exp(-dt/tltp));
                else
                    wdl4(t)=wdl4(t)*(1-altd*exp(dt/tltd));
                end
                %postl4d=-1;
                %pred=-1;
            end
        end
    end
end

figure();
subplot(3,1,1);
plot(wsp);
title('SP-->L4');
xlabel('Time ( in ms )');
ylabel('w (sp->l4)');

subplot(3,1,2);
plot(wsl4);
title('S-->L4');
xlabel('Time ( in ms )');
ylabel('w (s->l4)');

subplot(3,1,3);
plot(wdl4);
title('D-->L4');
xlabel('Time ( in ms )');
ylabel('w (d-->l4)');

WSP=wsp;
WSL4=wsl4;
WDL4=wdl4;

t=1;
for i = 1:15
    for j = 1:50
        if i ~= 8
            rS_ODD(t) = 10;
            rD_ODD(t) = 2.5;
        else
            rS_ODD(t) = 2.5;
            rD_ODD(t) = 10;
        end
        t = t + 1;
    end
    for j = 1:250
        rS_ODD(t) = 0.5;
        rD_ODD(t) = 0.5;
        t = t + 1;
    end
end

n = t- 1;
t = t - 1;

P_S_ODD = rS_ODD.*0.001;
P_D_ODD = rD_ODD.*0.001;

size(wsl4)

L4M = [];
SPM = [];
length = n;

timeval=[10000,25000,35000,40000,52000];

for val=1:5
    for loop = 1:50
        for i = 1:t           
            SpikeS_ODD(i) = rand < P_S_ODD(i);            
        end        
        for i = 1:t  
            SpikeD_ODD(i) = rand < P_D_ODD(i);         
        end     
        for  i = 1:t      
            if SpikeS_ODD(i) == 1         
                M_S_ODD(i) = 1;                      
                i = i + 1;                      
                M_S_ODD(i) = 1;                   
            else              
                M_S_ODD(i) = 0;               
            end           
        end        
        for  i = 1:t        
            if SpikeD_ODD(i) == 1                    
                M_D_ODD(i) = 1;                        
                i = i + 1;                        
                M_D_ODD(i) = 1;                    
            else                
                M_D_ODD(i) = 0;
            end            
        end
        %initializing the values
        timesteps = t;
        tsyn = 10;
        B = 5;
        tref = 2;
        ws = 0.2.*ones(1,timesteps);
        wd = 0.2.*ones(1,timesteps);
        wsp=WSP(timeval(val)).*ones(1,timesteps);
        wsl4=WSL4(timeval(val))*ones(1,timesteps);
        wdl4=WDL4(timeval(val)).*ones(1,timesteps);;
        %all time constants are in milli seconds
        S.xr = 1;
        S.xi = 0;
        S.xe = 0;
        S.tre = 0.9;
        S.tei = 10;
        S.tir = 3000;
        S.M = M_S_ODD;            
        S.g = zeros(1,timesteps);
        flag = 0;
        for i = 1:timesteps
            if flag == 1                                
                S.g(i) = S.g(i-1) + 1;
                flag = 0;
            else
                 if i ~= 1
                     S.g(i) = S.g(i-1)*exp(-1/tsyn);
                 end
            end
            if S.M(i) == 1
                flag = 1;
            end
        end
        D.xr = 1;
        D.xi = 0;
        D.xe = 0;
        D.tre = 0.9;
        D.tei = 10;
        D.tir = 3000;
        D.M = M_D_ODD;
        D.g = zeros(1,timesteps);
        flag = 0;
        for i = 1:timesteps
            if flag == 1
                D.g(i) = D.g(i-1) + 1;
                flag = 0;
            else
                if i ~= 1
                    D.g(i) = D.g(i-1)*exp(-1/tsyn);
                end
            end     
            if D.M(i)==1
                flag=1;
            end
        end
        SP.xr = 1;
        SP.xi = 0;
        SP.xe = 0;
        SP.tre = 0.9;
        SP.tei = 27;
        SP.tir = 3000;
        SP.M = zeros(1,timesteps);
        SP.g = zeros(1,timesteps);
        countervsp = 0;
        countervl4 = 0;
        Vsp = zeros(1,timesteps+1);
        Vl4 = zeros(1,timesteps);
        Vl4M = zeros(1,timesteps);
        flag=0;%For SP.g update
    
        tm= linspace(1,t,t);
        options =odeset('MaxStep',5);
        [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);   
        S.xr = interp1(T_S,S.x(:,1),tm,'spline');

        S.xe = interp1(T_S,S.x(:,2),tm,'spline');

        S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    
        [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
        D.xr = interp1(T_S,D.x(:,1),tm,'spline');
    
        D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    
        D.xi = interp1(T_S,D.x(:,3),tm,'spline');

    
        for t = 1:timesteps
            if flag == 1
                SP.g(t) = SP.g(t-1) + 1;
                flag = 0;
            else
                 if t ~= 1
                 SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
                 end
            end
            if SP.M(t) == 1
                flag = 1;
            end
            %Calculate the changes and assign them to a temporary variable
            if countervsp == 0
                Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
            else
                if SP.M(t) == 1
                    Vsp(t+1) = -B;
                else
                    Vsp(t+1) = Vsp(t)*exp(-1/tref);
                end
            countervsp = countervsp - 1;
            end
        
            if Vsp(t+1) > 0.05
                SP.M(t+1) = 1;
                countervsp = 20;
            end
        end
        [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
        SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');


        SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');

        SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
     
        for t = 1:timesteps
            if countervl4==0
                Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
            else
                if Vl4M(t)==1
                    Vl4(t+1)=-B;
                else
                    Vl4(t+1)=Vl4(t)*exp(-1/tref);
                end
                countervl4=countervl4-1;
            end
            if Vl4(t+1)>0.05
                Vl4M(t+1)=1;
                countervl4=20;
            end
        end
        SPM = [SPM; SP.M];
        L4M = [L4M;Vl4M];   
    end
    PSTHsp=sum(SPM)./(loop);
    PSTHl4=sum(L4M)./(loop);
    psthsp=[];
    psthl4=[];
    length = t;
    time=[];
    for i=1:length/10
        sum1=0;
        sum2=0;
            for j=1:10
            sum1=sum1+PSTHsp((i-1)*10+j);
            sum2=sum2+PSTHl4((i-1)*10+j);
            end
            psthsp=[psthsp,sum1];  
            psthl4=[psthl4,sum2];
            time(i)=(i-1)*10;
    end

    psthsp=psthsp.*100;
    psthl4=psthl4.*100;
    figure();
    subplot(2,1,1);
    plot(psthsp);
    title('PSTH SP neuron');
    xlabel('Time/10 ( in ms )');
    ylabel('Spike rate ( in #spikes per second )');

    subplot(2,1,2);
    plot(psthl4);

    title('PSTH L4 neuron');

    xlabel('Time/10 ( in ms )');

    ylabel('Spike rate ( in #spikes per second )');
end












%Question 5

%%Stimulus S and D responses
t=1;
diffstim=1200;
seq=[];
for i=1:diffstim
    base=rand();
    if(base<0.5)
        seq=[seq,1];%S=1, D=0
    else
        seq=[seq,0];
    end
end


for i=1:diffstim
    for j=1:50
        if seq(i)==1
            rS(t)=10;
            rD(t)=2.5;
        else
            rS(t)=2.5;
            rD(t)=10;
        end
        t=t+1;
    end
    for j=1:250
        rS(t)=0.5;
        rD(t)=0.5;
        t=t+1;
    end
end

n=t-1;
t = t- 1;

P_S = rS.*0.001;
P_D = rD.*0.001;

for i = 1:t;
    SpikeS(i) = rand < P_S(i);
end

for i = 1:t;
    SpikeD(i) = rand < P_D(i);
end


for  i = 1:t;
    if SpikeS(i) == 1
        M_S(i) = 1;
        i = i + 1;
        M_S(i) = 1;
    else
        M_S(i) = 0;
    end
end

for  i = 1:t;
    if SpikeD(i) == 1
        M_D(i) = 1;
        i = i + 1;
        M_D(i) = 1;
    else
        M_D(i) = 0;
    end
end

for loop = 1:1;
    %initializing the values
    timesteps = t;
    tsyn = 10;
    B = 5;
    tref = 2;
    ws = 0.2.*ones(1,timesteps);
    wd = 0.2.*ones(1,timesteps);
    wsp = 0.11.*ones(1,timesteps);
    wsl4 = 0.02.*ones(1,timesteps);
    wdl4 = 0.02.*ones(1,timesteps); %all time constants are in milli seconds
    tltp=13;
    tltd=20;
    altp=0.015;
    altd=0.021;
    postl4sp = -1;%post spike time for sp on L4 neuron
    postl4s = -1;%post spike time for s on L4 neuron
    postl4d = -1;%post spike time for d on L4 neuron
    pres = -1;%pre spike time of s on l4
    pred = -1;%pre spike time of d on l4
    presp =- 1;%pre spike time of sp on l4

    S.xr = 1;
    S.xi = 0;
    S.xe = 0;
    S.tre = 0.9;
    S.tei = 10;
    S.tir = 100;
    S.M = M_S;
    
    
    S.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            S.g(i) = S.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                S.g(i) = S.g(i-1)*exp(-1/tsyn);
            end
        end
        if S.M(i) == 1
            flag = 1;
        end
    end
    D.xr = 1;
    D.xi = 0;
    D.xe = 0;
    D.tre = 0.9;
    D.tei = 10;
    D.tir = 100;
    D.M = M_D;
    D.g = zeros(1,timesteps);
    flag = 0;
    for i = 1:timesteps
        if flag == 1
            D.g(i) = D.g(i-1) + 1;
            flag = 0;
        else
            if i ~= 1
                D.g(i) = D.g(i-1)*exp(-1/tsyn);
            end
        end     
        if D.M(i)==1
            flag=1;
        end
    end
    SP.xr = 1;
    SP.xi = 0;
    SP.xe = 0;
    SP.tre = 0.9;
    SP.tei = 27;
    SP.tir = 100;
    SP.M = zeros(1,timesteps);
    SP.g = zeros(1,timesteps);
    countervsp = 0;
    countervl4 = 0;
    Vsp = zeros(1,timesteps+1);
    Vl4 = zeros(1,timesteps);
    Vl4M = zeros(1,timesteps);
    flag=0;%For SP.g update
    
    tm= linspace(1,t,t);
    options =odeset('MaxStep',5);
    [T_S, S.x] = ode45(@(t,y) myode(t, y, double(S.M), double(S.tre), double(S.tir), double(S.tei), tm),[1 4500],[0.05 0.15 0.8],options);   
    S.xr = interp1(T_S,S.x(:,1),tm,'spline');

    S.xe = interp1(T_S,S.x(:,2),tm,'spline');

    S.xi = interp1(T_S,S.x(:,3),tm,'spline');
    
    [T_S, D.x] = ode45(@(t,y) myode(t, y, double(D.M), double(D.tre), double(D.tir), double(D.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    
    D.xr = interp1(T_S,D.x(:,1),tm,'spline');
    
    D.xe = interp1(T_S,D.x(:,2),tm,'spline');
    
    D.xi = interp1(T_S,D.x(:,3),tm,'spline');
    
    
    for t = 1:timesteps
        if flag == 1
            SP.g(t) = SP.g(t-1) + 1;
            flag = 0;
        else
            if t ~= 1
                SP.g(t) = SP.g(t-1)*exp(-1/tsyn);
            end
        end
        if SP.M(t) == 1
            flag = 1;
        end
        %Calculate the changes and assign them to a temporary variable
        if countervsp == 0
            Vsp(t+1) = S.g(t)*ws(t)*(double(S.xe(t)))+D.g(t)*wd(t)*(double(D.xe(t)))-0.1*Vsp(t);
        else
            if SP.M(t) == 1
                Vsp(t+1) = -B;
            else
                Vsp(t+1) = Vsp(t)*exp(-1/tref);
            end
        countervsp = countervsp - 1;
        end
        
        if Vsp(t+1) > 0.05
            SP.M(t+1) = 1;
            countervsp = 20;
        end
    end
    [T_S, SP.x] = ode45(@(t,y) myode(t, y, double(SP.M), double(SP.tre), double(SP.tir), double(SP.tei), tm),[1 4500],[0.05 0.15 0.8],options);
    SP.xr = interp1(T_S,SP.x(:,1),tm,'spline');


    SP.xe = interp1(T_S,SP.x(:,2),tm,'spline');

    SP.xi = interp1(T_S,SP.x(:,3),tm,'spline');
    
    for t = 1:timesteps
        if countervl4==0
            Vl4(t+1)=S.g(t)*wsl4(t)*S.xe(t)+D.g(t)*wdl4(t)*D.xe(t)+SP.g(t)*wsp(t)*SP.xe(t)-0.1*Vl4(t);
        else
            if Vl4M(t)==1
                Vl4(t+1)=-B;
            else
                Vl4(t+1)=Vl4(t)*exp(-1/tref);
            end
            countervl4=countervl4-1;
        end
        if Vl4(t+1)>0.05
            Vl4M(t+1)=1;
            countervl4=20;
        end
    end
    for t=1:timesteps        %Updating weights for LTP/LTD
        if(Vl4M(t)==1)
            postl4sp = t;
            postl4s = t;
            postl4d = t;
        end
        if(SP.M(t)==1)
            presp=t;
        end
        if(S.M(t)==1)
            pres=t;
        end
        if(D.M(t)==1)
            pred=t;
        end
        if(t~=1)
            wsp(t) = wsp(t-1);
            wsl4(t)= wsl4(t-1);
            wdl4(t)= wdl4(t-1);
        end
        if(postl4sp~=-1)
            if(presp~=-1&&(Vl4M(t)==1||SP.M(t)==1))
                dt=postl4sp-presp;
                if(dt>0)
                    wsp(t)=wsp(t)*(1+altp*exp(-dt/tltp));
                else
                    wsp(t)=wsp(t)*(1-altd*exp(dt/tltd));
                end
                %postl4sp=-1;
                %presp=-1;
            end
        end      
        if(postl4s~=-1)
            if(pres~=-1&&(Vl4M(t)==1||S.M(t)==1))
                dt=postl4s-pres;
                if(dt>0)
                    wsl4(t)=wsl4(t)*(1+altp*exp(-dt/tltp));
                else
                    wsl4(t)=wsl4(t)*(1-altd*exp(dt/tltd));
                end
                %postl4s=-1;
                %pres=-1;
            end
        end
        if(postl4d~=-1)
            if(pred~=-1&&(Vl4M(t)==1||D.M(t)==1))
                dt=postl4d-pred;
                if(dt>0)
                    wdl4(t)=wdl4(t)*(1+altp*exp(-dt/tltp));
                else
                    wdl4(t)=wdl4(t)*(1-altd*exp(dt/tltd));
                end
                %postl4d=-1;
                %pred=-1;
            end
        end
    end
end

figure();
subplot(3,1,1);
plot(wsp);
title('SP-->L4');
xlabel('Time ( in ms )');
ylabel('w (sp->l4)');

subplot(3,1,2);
plot(wsl4);
title('S-->L4');
xlabel('Time ( in ms )');
ylabel('w (s->l4)');

subplot(3,1,3);
plot(wdl4);
title('D-->L4');
xlabel('Time ( in ms )');
ylabel('w (d-->l4)');















function [spikeMat, tVec, r, spikes] = poissonSpikeGen(fr, tSim, nTrials, dt)
nBins = floor(tSim/dt);
r = ceil(rand(1, nBins).*10+40); %Average rate function
p = r.*0.001; % Probability of spikes in each bin

spikes = rand(1, nBins) < p; % Spike with given probability
spikes = double(spikes);
rnd = rand(nTrials, nBins); %Spike with given probability for more than one trials
for i = 1:nTrials
    spikeMat(i,:) = rnd(i,:) < p;
end
spikeMat = double(spikeMat);
tVec = 0:tSim;
end

function dydt = myode(t, y, M, tre, tir, tei, tm)
    dydt = zeros(3,1);
    f = interp1(tm, M, t); % Interpolate the data set (ft,f) at time t
    dydt(1)= -f*y(1)/tre + y(3)/tir;
    dydt(2)= f*y(1)/tre - y(2)/tei;
    dydt(3)= y(2)/tei - y(3)/tir;
end

end

