%HOMEWORK 3  -  Ghiotti Maddalena s281604, Sambo Aldo s282334

clc
clear all
global DueDate n s sinit p m

%% DATI

m = 5; %numero macchine
n = 100; %numero job
numSwitches = 50; %numero scambi nell'algoritmo Rand
numMovings = 50; %numero spostamenti nell'algoritmo Rand
numMultistart = 20; %numero sequenze di partenza per gli algoritmi multistart

%(codice fornito dal professore)
p = randi([20,100], m, n); %tempi di processo
s = randi([12,24], n, n); %tempi di setup
for i = 2:m
    s(:,:,i) = randi([12,24], n, n);
end
smin = min(s); %tempo minimo di setup: utile per generare DueDate sensate
sinit = randi([12,24], m,n); % tempi di setup al tempo 0 (per inizializzare la prima volta la macchina)
for i = 1:n
    s(i,i) = 0; %diagonale = 0
end
smedia = 18;
pmedia = 60;
l = 0; %l e' un coefficiente utile a generare delle DueDate sensate
for i = 1:n
    l = l + p(i) +smin(i);
end
l = l/m;
tau = 0.4; %parametro per le DueDate
R = 0.8; %parametro per le DueDate
edx = l*((1-tau+R)/2);
DueDate = randi(floor(edx), n, 1);
Job = zeros(n,1);
for i = 1:n
    Job(i) = i;
end

%Stampa dei dati fondamentali modificabili
fprintf("Numero di job: %.0f\n", n);
fprintf("Numero di macchine: %.0f\n", m);
fprintf("Numero di scambi nell'algoritmo Rand: %.0f\n", numSwitches);
fprintf("Numero di spostamenti nell'algoritmo Rand: %.0f\n", numMovings);
fprintf("Numero di sequenze di partenza per gli algoritmi multistart: %.0f\n", numMultistart);
fprintf("\n")
fprintf("\n")

%% EURISTICA COSTRUTTIVA

tic
[~, orderDueDate] = sort(DueDate);  %sequenza dei job con date di consegna ordinate.
[latenessMaxDueDate, jobLMaxDueDate, configurationDueDate, completingTimesDueDate] = orderSplitCompletingTimes(orderDueDate);
timeDueDate = toc;
disp("Euristica costruttiva: EDD (per min tempi di completamento)");
formatPrint(configurationDueDate, latenessMaxDueDate, timeDueDate)


%% EURISTICHE ITERATIVE

%Configurazione di partenza. 
startingConfiguration = randomConfiguration();
startingLateness = latenessMax(startingConfiguration);

%Stampo la configurazione di partenza.
    fprintf("CONFIGURAZIONE DI PARTENZA RANDOM");
    fprintf("\n")
    fprintf("La configurazione di partenza è:\n");
    for i = 1:m
        fprintf("Macchina %d:", i);
        for j = 1:n
            if startingConfiguration(i,j) ~= 0 
                fprintf("\t %d", startingConfiguration(i,j));
            end
        end
        fprintf("\n")
    end
    fprintf("\n")
    fprintf("La lateness max di partenza è: %d.\n", startingLateness)
    fprintf("\n")
    fprintf("\n")

%Euristica: LATENESS MASSIMA
tic
[lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(startingConfiguration);
timeLatenessMax = toc;
disp("Euristica iterativa: LATENESS MASSIMA");
percentageImprovementLatenessMax = formatPrintI(configurationLatenessMax, lMaxLatenessMax, startingLateness, timeLatenessMax);

%Euristica: SETUP MASSIMO
tic
[lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(startingConfiguration);
timeSetUpMax = toc;
disp("Euristica iterativa: SETUP MASSIMO");
percentageImprovementSetUpMax = formatPrintI(configurationSetUpMax, lMaxSetUpMax, startingLateness, timeSetUpMax);

%Euristica: RANDOM
tic
[lMaxRand, configurationRand] = algorithmRand(startingConfiguration, numSwitches, numMovings);
timeRand = toc;
disp("Euristica iterativa: RANDOM");
percentageImprovementRand = formatPrintI(configurationRand, lMaxRand, startingLateness, timeRand);


%% EURISTICHE ITERATIVE - MULTISTART

%Inizializzo matrici in cui inserire i dati ottenuti nelle iterazioni.
lMaxMatrix = zeros(3, numMultistart);
configurationMatrix = zeros(m, n, numMultistart, 3);

for t = 1:numMultistart
    %Configurazione di partenza. 
    startingConfiguration = randomConfiguration();
    startingLateness = latenessMax(startingConfiguration);
    
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(startingConfiguration);
    lMaxMatrix(1,t) = lMaxLatenessMax;
    configurationMatrix(:,:,t,1) = configurationLatenessMax;

    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(startingConfiguration);
    lMaxMatrix(2,t) = lMaxSetUpMax;
    configurationMatrix(:,:,t,2) = configurationSetUpMax;

    %Euristica: RANDOM
    [lMaxRand, configurationRand] = algorithmRand(startingConfiguration, numSwitches, numMovings);
    lMaxMatrix(3,t) = lMaxRand;
    configurationMatrix(:,:,t,3) = configurationRand;
end
[lMax, indexLMax] = min(lMaxMatrix,[], 2);   %minimo delle lateness massime ottenute e indici delle rispettive sequenze.

disp("Euristica iterativa multistart: LATENESS MASSIMA");
formatPrint(configurationMatrix(:,:,indexLMax(1), 1), lMax(1), timeLatenessMax*numMultistart);
disp("Euristica iterativa multistart: SETUP MASSIMO");
formatPrint(configurationMatrix(:,:,indexLMax(2), 2), lMax(2), timeSetUpMax*numMultistart);
disp("Euristica iterativa multistart: RANDOM");
formatPrint(configurationMatrix(:,:,indexLMax(3), 3), lMax(3), timeRand*numMultistart);


%% EURISTICHE MISTE

%Configurazione di partenza: EDD.
startingConfiguration = configurationDueDate;
startingLateness = latenessMaxDueDate;

%% EDD + Euristiche iterative

%Euristica: LATENESS MASSIMA
tic
[lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(startingConfiguration);
timeMixLMax = timeDueDate + toc;
disp("Euristica mista: EDD + LATENESS MASSIMA");
formatPrint(configurationLatenessMax, lMaxLatenessMax, timeMixLMax);

%Euristica: SETUP MASSIMO
tic
[lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(startingConfiguration);
timeSetUpMax = timeDueDate + toc;
disp("Euristica mista: EDD + SETUP MASSIMO");
formatPrint(configurationSetUpMax, lMaxSetUpMax, timeSetUpMax);

%Euristica: RANDOM
tic
[lMaxRand, configurationRand] = algorithmRand(startingConfiguration, numSwitches, numMovings);
timeRand = timeDueDate + toc;
disp("Euristica mista: EDD + RANDOM");
formatPrint(configurationRand, lMaxRand, timeRand);

%% EDD + Setup massimo + Lateness massima

tic
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(startingConfiguration);
%Euristica: LATENESS MASSIMA
[lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);

improving =  lMaxSetUpMax - lMaxLatenessMax;

while improving > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(configurationLatenessMax);
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);

    improving =  lMaxSetUpMax - lMaxLatenessMax;
end
timeSetupLateness = timeDueDate + toc;
disp("Euristica mista: EDD + SETUP MASSIMO + LATENESS MASSIMA");
formatPrint(configurationLatenessMax, lMaxLatenessMax, timeSetupLateness);

%% EDD + Lateness massima + Setup massimo

tic
%Euristica: LATENESS MASSIMA
[lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(startingConfiguration);
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(configurationLatenessMax);

improving =  lMaxLatenessMax - lMaxSetUpMax;

while improving > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(configurationLatenessMax);

    improving =  lMaxLatenessMax - lMaxSetUpMax;
end
timeLatenessSetup = timeDueDate + toc;
disp("Euristica mista: EDD + LATENESS MASSIMA + SETUP MASSIMO");
formatPrint(configurationSetUpMax, lMaxSetUpMax, timeLatenessSetup);

%% EDD + Setup massimo + Lateness massima + Random

tic
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(startingConfiguration);
%Euristica: LATENESS MASSIMA
[~, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);
%Euristica: RANDOM
[lMaxRand, configurationRand] = algorithmRand(configurationLatenessMax, numSwitches, numMovings);

improving =  lMaxSetUpMax - lMaxRand;

while improving > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(configurationRand);
    %Euristica: LATENESS MASSIMA
    [~, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);    
    %Euristica: RANDOM
    [lMaxRand, configurationRand] = algorithmRand(configurationLatenessMax, numSwitches, numMovings);

    improving =  lMaxSetUpMax - lMaxRand;
end
timeSetupLatenessRand = timeDueDate + toc;
disp("Euristica mista: EDD + SETUP MASSIMO + LATENESS MASSIMA + RANDOM");
formatPrint(configurationRand, lMaxRand, timeSetupLatenessRand);

%% EDD + (Setup massimo + Lateness massima) + Random

tic
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(startingConfiguration);
%Euristica: LATENESS MASSIMA
[lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);

improving =  lMaxSetUpMax - lMaxLatenessMax;

while improving > 0   %itero i primi due algoritmi fino a quando la lateness massima smette di diminuire.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(configurationLatenessMax);
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);

    improving =  lMaxSetUpMax - lMaxLatenessMax;
end
%Euristica: RANDOM
[lMaxRand, configurationRand] = algorithmRand(configurationLatenessMax, numSwitches, numMovings);

improvingExt = lMaxLatenessMax - lMaxRand;

while improvingExt > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(configurationRand);
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);

    improving =  lMaxSetUpMax - lMaxLatenessMax;

    while improving > 0   %itero i primi due algoritmi fino a quando la lateness massima smette di diminuire.
        %Euristica: SETUP MASSIMO
        [lMaxSetUpMax, configurationSetUpMax] = algorithmSetUpMax(configurationLatenessMax);
        %Euristica: LATENESS MASSIMA
        [lMaxLatenessMax, configurationLatenessMax] = algorithmLatenessMax(configurationSetUpMax);

        improving =  lMaxSetUpMax - lMaxLatenessMax;
    end
    %Euristica: RANDOM
    [lMaxRand, configurationRand] = algorithmRand(configurationLatenessMax, numSwitches, numMovings);

    improvingExt = lMaxLatenessMax - lMaxRand;
end
timeSetupLatenessRand2 = timeDueDate + toc;
disp("Euristica mista: EDD + (SETUP MASSIMO + LATENESS MASSIMA) + RANDOM");
formatPrint(configurationRand, lMaxRand, timeSetupLatenessRand2);








%% FUNZIONE: DA ORDINE UNICO A SEQUENZE SEPARATE PER TEMPI DI COMPLETAMENTO
% @input order, vettore contenente i job nell'ordine di partenza con cui inserirli.
% @output lMax, lateness massima della configurazione di arrivo.
% @output jobLMax, job con lateness massima. 
% @output configuration, configurazione di arrivo.
% @output completingTimes, vettore contenente i tempi di completamento dei job.
function [lMax, jobLMax, configuration, completingTimes] = orderSplitCompletingTimes(order)
    global DueDate n s sinit p m
    configuration = zeros(m, n);
    completingTimes = zeros(n, 1);
    %Informazioni sull'ultimo job aggiunto per macchina.
    lastAdded = zeros(m, 1);
    completingTimeLastAdded = zeros(m, 1);
    positionLastAdded = zeros(m,1);
    %Per ogni job calcolo il tempo di completamento che avrebbe in ogni macchina, seleziono il minore e inserisco il job nella sequenza.
    for q = 1:n
        currentCompletingTimes = zeros(m,1);
        for r = 1:m
            if lastAdded(r) == 0
                currentCompletingTimes(r) = sinit(r,order(q)) + p(r,order(q));
            else 
                currentCompletingTimes(r) = completingTimeLastAdded(r) + s(lastAdded(r),order(q),r) + p(r,order(q));
            end
        end
        [minCompleting, machineMinCompleting] = min(currentCompletingTimes);
        configuration(machineMinCompleting, positionLastAdded(machineMinCompleting)+1) = order(q);
        %Aggiornamento informazioni ultimo job aggiunto per macchina.
        lastAdded(machineMinCompleting) = order(q);
        completingTimeLastAdded(machineMinCompleting) = minCompleting;
        positionLastAdded(machineMinCompleting) = positionLastAdded(machineMinCompleting)+1;
        completingTimes(order(q)) = minCompleting;
    end
    lateness = completingTimes - DueDate;
    [lMax, jobLMax] = max(lateness);
end

%% FUNZIONE: SINGOLA ITERAZIONE LATENESS MASSIMA
% @input configuration, matrice rappresentante una certa configurazione di job.
% @output currentConfiguration, configurazione migliore (con minore lateness max) tra quelle ottenute con anticipazione del job con lateness max o scambio con job vicini temporalmente.
% @output currentLMax, lateness massima relativa a currentConfiguration.
function [currentConfiguration, currentLMax] = singleIterationLatenessMax(configuration)
    global n m
    %Trovo la posizione del job con la lateness massima.
    [~, jobLMax, ~, completingTimes] = latenessMax(configuration);
    [machineJobMax, positionJobMax] = find(configuration == jobLMax);
    %TROVO TUTTE LE PERMUTAZIONI DELL'INTORNO.
    permutationMatrix = [];
    %Anticipazione del job nella macchina stessa.
    if positionJobMax ~= 1
        for u = 1:positionJobMax-1
            permutation = configuration;   
            permutation(machineJobMax, u:positionJobMax) = [jobLMax, permutation(machineJobMax, u:positionJobMax-1)];
            permutationMatrix(:,:,u) = permutation;
        end
    end
    u = positionJobMax;
    for i = 1:m
        if i ~= machineJobMax
            j = 1;
            %Anticipazione del job spostandolo di macchina (fino alla posizione subito precedente l'ultimo job con tempo di completamento strettamente inferiore).
            while configuration(i,j) ~= 0 && completingTimes(configuration(i,j)) < completingTimes(jobLMax)
                permutation = configuration;
                permutation(i,j) = jobLMax;
                permutation(i,j+1:n) = configuration(i,j:n-1);
                permutation(machineJobMax, positionJobMax:n-1) = configuration(machineJobMax, positionJobMax+1:n);
                permutationMatrix(:,:,u) = permutation;
                u = u+1;
                j = j+1;
            end
            %Spostamento del job in un'altra macchina, in posizione subito precedente il primo job (della macchina di arrivo) con tempo
            %di completamento maggiore se questo esiste, altrimenti alla fine della sequenza (in sostituzione del primo 0).
            permutation = configuration;
            permutation(i,j) = jobLMax;
            permutation(i,j+1:n) = configuration(i,j:n-1);
            permutation(machineJobMax, positionJobMax:n) = [configuration(machineJobMax, positionJobMax+1:n),0];
            permutationMatrix(:,:,u) = permutation;
            u = u+1;
            %Scambio con l'ultimo job con tempo di completamento strettamente inferiore nella sequenza di un'altra macchina.
            if j ~= 1
                permutation = configuration;
                permutation(i,j-1) = jobLMax;
                permutation(machineJobMax, positionJobMax) = configuration(i,j-1);
                permutationMatrix(:,:,u) = permutation;
                u = u+1;
            end
        end
    end
    [currentLMax, ~, currentConfiguration] = latenessMax(permutationMatrix);
end 

%% FUNZIONE: ALGORITMO LATENESS MASSIMA
% @input startingConfiguration, configurazione di partenza da cui comincia l'esecuzione dell'algoritmo.
% @output lMax, lateness massima migliore (più piccola) raggiunta dall'algoritmo.
% @output configuration, configurazione migliore tra quelle testate, cioè con lateness massima lMax. 
function [lMax, configuration] = algorithmLatenessMax(startingConfiguration)
    lMax = latenessMax(startingConfiguration);
    configuration = startingConfiguration;
    [currentConfiguration, currentLMax] = singleIterationLatenessMax(startingConfiguration);  %eseguo la prima iterazione cominciando dalla configurazione di partenza. 
    while currentLMax < lMax   %continuo ad iterare solo se la lateness massima continua ad essere migliore.
        configuration = currentConfiguration;
        lMax = currentLMax;
        [currentConfiguration, currentLMax] = singleIterationLatenessMax(configuration);
    end
end

%% FUNZIONE: SINGOLA ITERAZIONE SETUP MASSIMO
% @input configuration, matrice rappresentante una certa configurazione di job.
% @input completingTimes, vettore contenente i tempi di completamento dei job.
% @output currentConfiguration, configurazione migliore (con minore lateness max) tra quelle ottenute con scambi e spostamenti del job avente setup max con job vicini temporalmente.
% @output currentLMax, lateness massima relativa a currentConfiguration.
function [currentConfiguration, currentLMax, currentCompletingTimes] = singleIterationSetUpMax(configuration, completingTimes)
    global n s sinit m
    %Trovo la posizione del job con setup massimo.
    setUpVector = zeros(n,1);
    for i = 1:m
        j = 1;
        while configuration(i,j) ~= 0            
            if j == 1
                setUpVector(configuration(i,1)) = sinit(i, configuration(i,1));
            else 
                setUpVector(configuration(i,j)) = s(configuration(i,j-1), configuration(i,j),i);
            end
            j = j+1;
        end
    end
    [~, jobMaxSetUp] = max(setUpVector);
    [machineJobMax, positionJobMax] = find(configuration == jobMaxSetUp);
    %TROVO TUTTE LE PERMUTAZIONI DELL'INTORNO.
    permutationMatrix = [];
    %Permutazione dei due job precedenti il setup max.
    count = 1;
    if positionJobMax ~= 1 && positionJobMax ~= 2
        permutation1 = configuration;   
        permutation1(machineJobMax, positionJobMax-2:positionJobMax-1) = [configuration(machineJobMax, positionJobMax-1), configuration(machineJobMax, positionJobMax-2)];
        permutationMatrix(:,:,count) = permutation1;
        count = count+1;
    end
    %Permutazione dei due job adiacenti al setup max.
    if positionJobMax ~= 1
        permutation2 = configuration;   
        permutation2(machineJobMax, positionJobMax-1:positionJobMax) = [jobMaxSetUp, configuration(machineJobMax, positionJobMax-1)];
        permutationMatrix(:,:,count) = permutation2;
        count = count+1;
    end
    %Permutazione dei due job seguenti il setup max.
    if positionJobMax ~= n && configuration(machineJobMax, positionJobMax+1) ~= 0 
        permutation3 = configuration;   
        permutation3(machineJobMax, positionJobMax:positionJobMax+1) = [configuration(machineJobMax, positionJobMax+1), jobMaxSetUp];
        permutationMatrix(:,:,count) = permutation3;
        count = count+1;
    end 
    for i = 1:m
        if i ~= machineJobMax
            j = 1;
            while configuration(i,j) ~= 0 && completingTimes(configuration(i,j)) < completingTimes(jobMaxSetUp)
                j = j+1;
            end
            %Spostamento del job in un'altra macchina, in posizione subito precedente il primo job (della macchina di arrivo) con tempo
            %di completamento maggiore se questo esiste, altrimenti alla fine della sequenza (in sostituzione del primo 0).
            permutation6 = configuration;
            permutation6(i,j) = jobMaxSetUp;
            if j ~= n
                permutation6(i,j+1:n) = configuration(i,j:n-1);
            end
            permutation6(machineJobMax, positionJobMax:n) = [configuration(machineJobMax, positionJobMax+1:n),0];
            permutationMatrix(:,:,count) = permutation6;
            count = count+1;
            if j ~= 1
                %Spostamento del job in posizione subito precedente l'ultimo job con tempo di completamento strettamente inferiore nella sequenza di un'altra macchina.
                permutation5 = configuration;
                permutation5(i,j-1) = jobMaxSetUp;
                permutation5(i,j:n) = configuration(i,j-1:n-1);
                permutation5(machineJobMax, positionJobMax:n-1) = configuration(machineJobMax, positionJobMax+1:n);
                permutationMatrix(:,:,count) = permutation5;
                count = count+1;
                %Scambio con l'ultimo job con tempo di completamento strettamente inferiore nella sequenza di un'altra macchina.
                permutation4 = configuration;
                permutation4(i,j-1) = jobMaxSetUp;
                permutation4(machineJobMax, positionJobMax) = configuration(i,j-1);
                permutationMatrix(:,:,count) = permutation4;
                count = count+1;
            end
        end
    end
    [currentLMax, ~, currentConfiguration, currentCompletingTimes] = latenessMax(permutationMatrix);
end

%% FUNZIONE: ALGORITMO SETUP MASSIMO
% @input startingConfiguration, configurazione di partenza da cui comincia l'esecuzione dell'algoritmo.
% @output lMax, lateness massima migliore (più piccola) raggiunta dall'algoritmo.
% @output configuration, configurazione migliore tra quelle testate, cioè con lateness massima lMax. 
function [lMax, configuration] = algorithmSetUpMax(startingConfiguration)
    [lMax,~,~, startingCompletingTimes] = latenessMax(startingConfiguration);
    configuration = startingConfiguration;
    [currentConfiguration, currentLMax, currentCompletingTimes] = singleIterationSetUpMax(startingConfiguration, startingCompletingTimes);  %eseguo la prima iterazione cominciando dalla configurazione di partenza. 
    while currentLMax < lMax   %continuo ad iterare solo se la lateness massima continua ad essere migliore.
        configuration = currentConfiguration;
        lMax = currentLMax;
        [currentConfiguration, currentLMax, currentCompletingTimes] = singleIterationSetUpMax(configuration, currentCompletingTimes);
    end
end

%% FUNZIONE: SINGOLA ITERAZIONE RANDOM
% @input configuration, matrice rappresentante una certa configurazione di job.
% @input numSwitches, numero di sequenze vicine ottenute con scambi casuali. 
% @input numMovings, numero di sequenze vicine ottenute con spostamenti casuali.
% @output currentConfiguration, configurazione migliore (con minore lateness max) tra quelle ottenute con numSwitches scambi e numMovings spostamenti random.
% @output currentLMax, lateness massima relativa a currentConfiguration.
function [currentConfiguration, currentLMax] = singleIterationRand(configuration, numSwitches, numMovings)
    global n m
    permutationMatrix = [];
    
    %SCAMBIO due job casuali della configurazione.
    for i = 1:numSwitches
        a = randi([1,n]);
        b = randi([1,n-1]);   %evito di generare a=b
        if b >= a
            b = b+1;
        end
        [aMachine, aPosition] = find(configuration == a);
        [bMachine, bPosition] = find(configuration == b);
        switchConfiguration = configuration;
        switchConfiguration(aMachine, aPosition) = b;
        switchConfiguration(bMachine, bPosition) = a;
        permutationMatrix(:,:,i) = switchConfiguration;
    end
    
    %SPOSTO un job casuale della configurazione in un'altra posizione casuale.
    for i = numSwitches+1:numSwitches+numMovings
        a = randi([1,n]);   %posizione job da scambiare.
        b = randi([1,n+m-2]);   %posizione dove inserire il job.
        moveConfiguration = configuration;
        %Trovo le "coordinate" di a nella configurazione.
        count = 0;
        aMachine = 1;
        while count + nnz(configuration(aMachine,:)) < a
            count = count + nnz(configuration(aMachine,:));
            aMachine = aMachine+1;
        end
        aPosition = a-count;
        %Evito di tenere il job in a fermo nella stessa posizione.
        if b >= a+aMachine-1
            b = b+2;
        end
        %Trovo le "coordinate" di b nella configurazione.
        count = 0;
        bMachine = 1;
        while count+nnz(configuration(bMachine,:))+1 < b
            count = count+nnz(configuration(bMachine,:))+1;
            bMachine = bMachine+1;
        end
        bPosition = b-count;
        %Effettuo gli spostamenti.
        if aMachine == bMachine   %spostamento di un job nella macchina stessa.
            if bPosition > aPosition
                moveConfiguration(aMachine,aPosition:bPosition-1) = [moveConfiguration(aMachine,aPosition+1:bPosition-1),configuration(aMachine,aPosition)];
            else
                moveConfiguration(bMachine,bPosition:aPosition)=[configuration(aMachine,aPosition),configuration(bMachine,bPosition:aPosition-1)];
            end
        else   %spostamento di un job in un'altra macchina.
            moveConfiguration(bMachine,bPosition:n) = [configuration(aMachine,aPosition),moveConfiguration(bMachine,bPosition:n-1)];
            moveConfiguration(aMachine,aPosition:n) = [moveConfiguration(aMachine,aPosition+1:n),0];
        end
        permutationMatrix(:,:,i) = moveConfiguration;
    end
    [currentLMax, ~, currentConfiguration] = latenessMax(permutationMatrix);
end

%% FUNZIONE: ALGORITMO RANDOM
% @input startingConfiguration, configurazione di partenza da cui comincia l'esecuzione dell'algoritmo.
% @input numSwitches, numero di sequenze vicine ottenute con scambi casuali. 
% @input numMovings, numero di sequenze vicine ottenute con spostamenti casuali.
% @output lMax, lateness massima migliore (più piccola) raggiunta dall'algoritmo.
% @output configuration, configurazione migliore tra quelle testate, cioè con lateness massima lMax.
function [lMax, configuration] = algorithmRand(startingConfiguration, numSwitches, numMovings)
    lMax = latenessMax(startingConfiguration);
    configuration = startingConfiguration;
    [currentConfiguration, currentLMax] = singleIterationRand(startingConfiguration, numSwitches, numMovings);  %eseguo la prima iterazione cominciando dalla configurazione di partenza. 
    while currentLMax < lMax   %continuo ad iterare solo se la lateness massima continua ad essere migliore.
        configuration = currentConfiguration;
        lMax = currentLMax;
        [currentConfiguration, currentLMax] = singleIterationRand(configuration, numSwitches, numMovings);
    end
end

%% FUNZIONE: LATENESS MASSIMA
% @input Q, matrice tridimensionale con configurazioni di job di cui voglio calcolare la migliore lateness max.
% @output lMax, lateness massima migliore tra le lMax delle singole configurazioni.
% @output jobLMax, job con lateness massima. 
% @output configuration, configurazione di Q con lateness massima minima.
% @output completingTimes, vettore tempi di completamento.
function [lMax, jobLMax, configuration, completingTimes] = latenessMax(Q)
    global DueDate n s sinit p m
    vectorLMax = zeros(size(Q,3), 1);
    jobLMaxVector = zeros(size(Q,3), 1);
    for i = 1:size(Q,3)   %iterazione sulle configurazioni.
        completingTimes = zeros(n, 1); 
        for j = 1:m   %iterazione sulle macchine.
            k = 1;   %iterazione sulle posizioni dei job nelle singole sequenze.
            while Q(j,k,i) ~= 0
                if k == 1
                    completingTimes(Q(j,k,i)) = sinit(j,Q(j,k,i)) + p(j,Q(j,k,i));   %calcolo del tempo di completamento del primo job su una macchina.
                else
                    completingTimes(Q(j,k,i)) = completingTimes(Q(j,k-1,i)) + s(Q(j,k-1,i), Q(j,k,i),j) + p(j,Q(j,k,i));   %calcolo del tempo di completamento dei job successivi.
                end
                k = k + 1;
            end
        end
        L = completingTimes - DueDate;  %vettore delle lateness di ogni job.
        [vectorLMax(i), jobLMaxVector(i)] = max(L);
    end
    [lMax, indexConfiguration] = min(vectorLMax);
    configuration = Q(:,:,indexConfiguration);
    jobLMax = jobLMaxVector(indexConfiguration);
end

%% FUNZIONE: CONFIGURAZIONE RANDOM
% @output rs, matrice contenente configurazione casuale di n jobs.
function rs = randomConfiguration()
    global n m
    rs = zeros(m,n);
    t = randi([0,100], 1, n);
    [~,randomOrder] = sort(t);
    numPerMachine = floor(n/m);   %distribuisco i job tra le macchine in modo che il numero di quelli eseguiti su ogni macchina sia simile (in numero uguale sulle prime m-1, la m-esima macchina risente invece del resto della divisione)
    for i = 1:m-1
        rs(i,1:numPerMachine) = randomOrder((i-1)*numPerMachine+1:i*numPerMachine);
    end
    count = 1;
    for j = (m-1)*numPerMachine+1:n
        rs(m,count) = randomOrder(j);
        count = count+1;
    end    
end

%% FUNZIONE: STAMPA RISULTATI FORMATTATA (euristica iterativa)
% @input a, configurazione di arrivo.
% @input la, lateness max configurazione di arrivo.
% @input lp, lateness max configurazione di partenza.
% @input t, tempo di processamento dell'algoritmo.
% @output percentageImprovement, miglioramento percentuale della lateness max rispetto alla lateness max di partenza.
function percentageImprovement = formatPrintI(a, la, lp, t)
    global n m
    fprintf("La configurazione ottenuta è:\n");
    for i = 1:m
        fprintf("Macchina %d:", i);
        for j = 1:n
            if a(i,j) ~= 0 
                fprintf("\t %d", a(i,j));
            end
        end
        fprintf("\n")
    end
    fprintf("\n")
    fprintf("La lateness max ottenuta è: %d.\n", la)
    %Valuto il miglioramento sulla lateness massima.
    fprintf("Il miglioramento sulla lateness max è: %d.\n", la - lp)
    percentageImprovement = round((lp - la)*100/lp, 2);
    fprintf("Il miglioramento percentuale sulla lateness max è: %s%%.\n", num2str(percentageImprovement))
    %Stampo il tempo di esecuzione dell'algoritmo.
    fprintf("Il tempo di esecuzione dell'algoritmo è: %f secondi.\n", t)
    fprintf("\n")
    fprintf("\n")
end

%% FUNZIONE: STAMPA RISULTATI FORMATTATA
% @input s, configurazione.
% @input l, lateness max.
% @input t, tempo di esecuzione dell'algoritmo.
function formatPrint(s, l, t)
    global n m
    fprintf("La configurazione ottenuta è:\n");
    for i = 1:m
        fprintf("Macchina %d:", i);
        for j = 1:n
            if s(i,j) ~= 0 
                fprintf("\t %d", s(i,j));
            end
        end
        fprintf("\n")
    end
    fprintf("\n")
    fprintf("La lateness max ottenuta è: %d.\n", l)
    %Stampo il tempo di esecuzione dell'algoritmo.
    fprintf("Il tempo di esecuzione dell'algoritmo è: %f secondi.\n", t)
    fprintf("\n")
    fprintf("\n")
end