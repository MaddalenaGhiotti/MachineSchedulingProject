%HOMEWORK 2  -  Ghiotti Maddalena s281604, Sambo Aldo s282334

clc
clear all
global DueDate n s sinit p

%% DATI

m = 1; %numero macchine
n = 7; %numero job
numIterations = 4; %numero scambi nell'algoritmo Rand.
numMultistart = 10; %numero sequenze di partenza per gli algoritmi multistart.

M = 10000; %grande M
%(codice fornito dal professore)
p = randi([20,100], n, 1); %tempi di processo
s = randi([12,24], n, n); %tempi di setup
smin = min(s); %tempo minimo di setup: utile per generare DueDate sensate
sinit = randi([12,24], n,1); % tempi di setup al tempo 0 (per inizializzare la prima volta la macchina)
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
fprintf("Numero di scambi nell'algoritmo Rand: %.0f\n", numIterations);
fprintf("Numero di sequenze di partenza per gli algoritmi multistart: %.0f\n", numMultistart);
fprintf("\n")
fprintf("\n")


%% EURISTICHE COSTRUTTIVE

%EDD
tic
[~, sequenceDueDate] = sort(DueDate);  %sequenza dei job con date di consegna ordinate.
latenessDueDate = latenessMax(sequenceDueDate');
timeDueDate = toc;
disp("Euristica costruttiva: EDD");
formatPrint(sequenceDueDate, latenessDueDate, timeDueDate);

%SETUP
tic
sequenceSetUp = zeros(n, 1);
[~, jobMinSinit] = min(sinit);
sequenceSetUp(1) = jobMinSinit;
for q = 2:n
    [sortedRawTimes, sortedRawJobs] = sort(s(sequenceSetUp(q-1),:));
    count = 1;
    while ismember(sortedRawJobs(count), sequenceSetUp) == 1
        count = count + 1;
    end
    sequenceSetUp(q) = sortedRawJobs(count);
end
latenessSetUp = latenessMax(sequenceSetUp');
timeSetup = toc;
disp("Euristica costruttiva: SETUP");
formatPrint(sequenceSetUp, latenessSetUp, timeSetup);

%PROCESSING TIME (PROD)
tic
[~, sequenceProcess] = sort(p');   %sequenza dei job con processing time ordinati.
latenessProcess = latenessMax(sequenceProcess);
timeProcess = toc;
disp("Euristica costruttiva: PROCESSING TIME (PROD)");
formatPrint(sequenceProcess, latenessProcess, timeProcess);


%% EURISTICHE ITERATIVE 

%Sequenza di partenza: vettore RIGA. 
startingSequence = randomSequence();
startingLateness = latenessMax(startingSequence);
%Stampo la sequenza di partenza.
    fprintf("SEQUENZA DI PARTENZA RANDOM");
    fprintf("\n")
    fprintf("La sequenza di partenza è: ");
    for i = 1:n
        fprintf("\t %d", startingSequence(i));
    end
    fprintf("\n")
    fprintf("La lateness max di partenza è: %d.\n", startingLateness)
    fprintf("\n")
    fprintf("\n")

%Euristica: LATENESS MASSIMA
tic
[lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(startingSequence);
timeLatenessMax = toc;
disp("Euristica iterativa: LATENESS MASSIMA");
percentageImprovementLatenessMax = formatPrintI(sequenceLatenessMax, lMaxLatenessMax, startingLateness, timeLatenessMax);

%Euristica: SETUP MASSIMO
tic
[lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(startingSequence);
timeSetUpMax = toc;
disp("Euristica iterativa: SETUP MASSIMO");
percentageImprovementSetUpMax = formatPrintI(sequenceSetUpMax, lMaxSetUpMax, startingLateness, timeSetUpMax);

%Euristica: RANDOM
tic
[lMaxRand, sequenceRand] = algorithmRand(startingSequence, numIterations);
timeRand = toc;
disp("Euristica iterativa: RANDOM");
percentageImprovementRand = formatPrintI(sequenceRand, lMaxRand, startingLateness, timeRand);


%% EURISTICHE ITERATIVE - MULTISTART

%Inizializzo matrici in cui inserire i dati ottenuti nelle iterazioni.
lMaxMatrix = zeros(3, numMultistart);
sequenceMatrix = zeros(numMultistart, n, 3);

for t = 1:numMultistart
    %Sequenza di partenza: vettore RIGA. 
    startingSequence = randomSequence();
    startingLateness = latenessMax(startingSequence);
    
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(startingSequence);
    lMaxMatrix(1,t) = lMaxLatenessMax;
    sequenceMatrix(t, :, 1) = sequenceLatenessMax;

    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(startingSequence);
    lMaxMatrix(2,t) = lMaxSetUpMax;
    sequenceMatrix(t, :, 2) = sequenceSetUpMax;

    %Euristica: RANDOM
    [lMaxRand, sequenceRand] = algorithmRand(startingSequence, numIterations);
    lMaxMatrix(3,t) = lMaxRand;
    sequenceMatrix(t, :, 3) = sequenceRand;
end
[lMax, indexLMax] = min(lMaxMatrix,[], 2);   %minimo delle lateness massime ottenute e indici delle rispettive sequenze.

disp("Euristica iterativa multistart: LATENESS MASSIMA");
formatPrint(sequenceMatrix(indexLMax(1), :, 1), lMax(1), timeLatenessMax*numMultistart);
disp("Euristica iterativa multistart: SETUP MASSIMO");
formatPrint(sequenceMatrix(indexLMax(2), :, 2), lMax(2), timeSetUpMax*numMultistart);
disp("Euristica iterativa multistart: RANDOM");
formatPrint(sequenceMatrix(indexLMax(3), :, 3), lMax(3), timeRand*numMultistart);


%% EURISTICHE MISTE

%Sequenza di partenza: EDD.
startingSequence = sequenceDueDate';
startingLateness = latenessDueDate;

%% EDD + Euristiche iterative

%Euristica: LATENESS MASSIMA
tic
[lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(startingSequence);
timeMixLMax = timeDueDate + toc;
disp("Euristica mista: EDD + LATENESS MASSIMA");
formatPrint(sequenceLatenessMax, lMaxLatenessMax, timeMixLMax);

%Euristica: SETUP MASSIMO
tic
[lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(startingSequence);
timeSetUpMax = timeDueDate + toc;
disp("Euristica mista: EDD + SETUP MASSIMO");
formatPrint(sequenceSetUpMax, lMaxSetUpMax, timeSetUpMax);

%Euristica: RANDOM
tic
[lMaxRand, sequenceRand] = algorithmRand(startingSequence, numIterations);
timeRand = timeDueDate + toc;
disp("Euristica mista: EDD + RANDOM");
formatPrint(sequenceRand, lMaxRand, timeRand);

%% EDD + Setup massimo + Lateness massima

tic
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(startingSequence);

%Euristica: LATENESS MASSIMA
[lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);

improving =  lMaxSetUpMax - lMaxLatenessMax;

while improving > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(sequenceLatenessMax);
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);

    improving =  lMaxSetUpMax - lMaxLatenessMax;
end
timeSetupLateness = timeDueDate + toc;
disp("Euristica mista: EDD + SETUP MASSIMO + LATENESS MASSIMA");
formatPrint(sequenceLatenessMax, lMaxLatenessMax, timeSetupLateness);

%% EDD + Lateness massima + Setup massimo

tic
%Euristica: LATENESS MASSIMA
[lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(startingSequence);
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(sequenceLatenessMax);

improving =  lMaxLatenessMax - lMaxSetUpMax;

while improving > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(sequenceLatenessMax);

    improving =  lMaxLatenessMax - lMaxSetUpMax;
end
timeLatenessSetup = timeDueDate + toc;
disp("Euristica mista: EDD + LATENESS MASSIMA + SETUP MASSIMO");
formatPrint(sequenceSetUpMax, lMaxSetUpMax, timeLatenessSetup);

%% EDD + Setup massimo + Lateness massima + Random

tic
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(startingSequence);
%Euristica: LATENESS MASSIMA
[~, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);
%Euristica: RANDOM
[lMaxRand, sequenceRand] = algorithmRand(sequenceLatenessMax, numIterations);

improving =  lMaxSetUpMax - lMaxRand;

while improving > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(sequenceRand);
    %Euristica: LATENESS MASSIMA
    [~, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);  
    %Euristica: RANDOM
    [lMaxRand, sequenceRand] = algorithmRand(sequenceLatenessMax, numIterations);

    improving =  lMaxSetUpMax - lMaxRand;
end
timeSetupLatenessRand = timeDueDate + toc;
disp("Euristica mista: EDD + SETUP MASSIMO + LATENESS MASSIMA + RANDOM");
formatPrint(sequenceRand, lMaxRand, timeSetupLatenessRand);

%% EDD + (Setup massimo + Lateness massima) + Random

tic
%Euristica: SETUP MASSIMO
[lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(startingSequence);
%Euristica: LATENESS MASSIMA
[lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);

improving =  lMaxSetUpMax - lMaxLatenessMax;

while improving > 0   %itero i primi due algoritmi fino a quando la lateness massima smette di diminuire.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(sequenceLatenessMax);
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);

    improving =  lMaxSetUpMax - lMaxLatenessMax;
end
%Euristica: RANDOM
[lMaxRand, sequenceRand] = algorithmRand(sequenceLatenessMax, numIterations);

improvingExt = lMaxLatenessMax - lMaxRand;

while improvingExt > 0   %applico nuovamente gli algoritmi solo se la lateness massima è diminuita.
    %Euristica: SETUP MASSIMO
    [lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(sequenceRand);
    %Euristica: LATENESS MASSIMA
    [lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);

    improving =  lMaxSetUpMax - lMaxLatenessMax;

    while improving > 0   %itero i primi due algoritmi fino a quando la lateness massima smette di diminuire.
        %Euristica: SETUP MASSIMO
        [lMaxSetUpMax, sequenceSetUpMax] = algorithmSetUpMax(sequenceLatenessMax);
        %Euristica: LATENESS MASSIMA
        [lMaxLatenessMax, sequenceLatenessMax] = algorithmLatenessMax(sequenceSetUpMax);

        improving =  lMaxSetUpMax - lMaxLatenessMax;
    end
    %Euristica: RANDOM
    [lMaxRand, sequenceRand] = algorithmRand(sequenceLatenessMax, numIterations);

    improvingExt = lMaxLatenessMax - lMaxRand;
end
timeSetupLatenessRand2 = timeDueDate + toc;
disp("Euristica mista: EDD + (SETUP MASSIMO + LATENESS MASSIMA) + RANDOM");
formatPrint(sequenceRand, lMaxRand, timeSetupLatenessRand2);


%% MODELLO MILP

if n < 8
    fprintf("MODELLO MILP\n");
    tic
    setUpMatrix = [0 sinit'; zeros(n, 1) s];

    %PROBLEMA
    prob = optimproblem('ObjectiveSense', 'min');

    %VARIABLI
    latenessMaxMILP = optimvar('latenessMaxMILP', 1, 1);
    completingTime = optimvar('completingTime', 1, n, 'LowerBound', 0);
    sequenceJobLocation = optimvar('sequenceJobLocation', 1, n, 'LowerBound', 1, 'UpperBound', n, 'Type', 'integer');
    isImmediatelyPrevious = optimvar('isImmediatelyPrevious', n+1, n+1, 'LowerBound', 0, 'UpperBound', 1, 'Type', 'integer');
    isPrevious = optimvar('isPrevious', n, n, 'LowerBound', 0, 'UpperBound', 1, 'Type', 'integer');

    %FUNZIONE OBIETTIVO
    prob.Objective = latenessMaxMILP;

    %VINCOLI
    %vincolo sulla lateness max
    latenessMaxConstr = optimconstr(1, n);
    latenessMaxConstr = latenessMaxMILP*ones(1, n) >= completingTime - DueDate';
    prob.Constraints.latenessMaxConstr = latenessMaxConstr;

    %vincolo nodi in uscita
    outgoingJobConstr = optimconstr(n+1, 1);
    for i = 1:n + 1
        outgoingJobConstr(i) = sum(isImmediatelyPrevious(i, :)) == 1;
    end
    prob.Constraints.outgoingJobConstr = outgoingJobConstr;

    %vincolo nodi in entrata
    incomingJobConstr = optimconstr(n+1, 1);
    for i = 1:n + 1
        incomingJobConstr(i) = sum(isImmediatelyPrevious(:, i)) == 1;
    end
    prob.Constraints.incomingJobConstr = incomingJobConstr;

    %vincolo su tempi di completamento iniziale
    completingTimeInitConstr = optimconstr(1, n);
    for i=1:n
        completingTimeInitConstr(i) = p(i) + isImmediatelyPrevious(i,i+1)*setUpMatrix(1,i+1)<=completingTime(i);
    end
    prob.Constraints.completingTimeInitConstr = completingTimeInitConstr;

    %vincolo su matrice isImmediatelyPrevious
    nodeImmediatelyPreviousItselfConstr = optimconstr(1, n+1);
    for i = 1:n+1
        nodeImmediatelyPreviousItselfConstr(i) = isImmediatelyPrevious(i, i) == 0;
    end
    prob.Constraints.nodeImmediatelyPreviousItselfConstr = nodeImmediatelyPreviousItselfConstr;

    %vincolo in cui assegno a ciascun job il rispettivo posto nella sequenza (vincolo subtour)
    jobEtiquetteConstr = optimconstr(n, n);
    for i = 1:n
        for j = 1:n
            jobEtiquetteConstr(i, j) = sequenceJobLocation(i) - sequenceJobLocation(j) + n*isImmediatelyPrevious(i+1, j+1) <= n-1;
        end
    end
    prob.Constraints.jobEtiquetteConstr = jobEtiquetteConstr;

    %vincolo su matrice isPrevious
    previousNodeConstr = optimconstr(n, n);
    for i = 1:n
        for j = 1:n
            if j == i
                previousNodeConstr(i, j) = isPrevious(i, j) == 0;
            else
                previousNodeConstr(i, j) = isPrevious(i, j) == 1 - isPrevious(j, i);
            end
        end
    end
    prob.Constraints.previousNodeConstr = previousNodeConstr;

    %vincolo di linking tra isPrevious e sequenceJobLocation
    previousAndSequenceConstr = optimconstr(n, n);
    for i = 1:n
        for j = 1:n
            previousAndSequenceConstr(i, j) = sequenceJobLocation(j)  <= isPrevious(i, j)*M + sequenceJobLocation(i);
        end
    end
    prob.Constraints.previousAndSequenceConstr = previousAndSequenceConstr;

    %vincolo sui tempi di completamento
    completingTimeConstr = optimconstr(n, n);
    for i = 1:n
        for j = 1:n
            completingTimeConstr(i, j) = completingTime(i) + p(j) + sum(isImmediatelyPrevious(2:n+1, j+1).*setUpMatrix(2:n+1, j+1)) <= completingTime(j) + M*(1-isPrevious(i, j));
        end
    end
    prob.Constraints.completingTimeConstr = completingTimeConstr;

    %SOLUZIONE_PROBLEMA
    opts = optimoptions(@intlinprog, 'Display', 'off');
    [sol, value] = solve(prob, 'Options', opts);
    timeMILP = toc;
    formatPrint(sol.sequenceJobLocation, value, timeMILP);
end








%% FUNZIONE: SINGOLA ITERAZIONE LATENESS MASSIMA
% @input sequence, vettore rappresentante una certa sequenza di job.
% @output currentSequence, sequenza migliore (con minore lateness max) tra quelle ottenute con anticipazione del job con lateness max.
% @output currentLMax, lateness massima relativa a currentSequence.
function [currentSequence, currentLMax] = singleIterationLatenessMax(sequence)
    [lMaxSequence, indexLMax] = latenessMax(sequence);
    permutationMatrix = [];
    if indexLMax ~= 1
        for i = 1:indexLMax-1   
            permutation = sequence;   
            permutation(i:indexLMax) = [permutation(indexLMax), permutation(i:indexLMax-1)];
            permutationMatrix = [permutationMatrix; permutation];
        end
        [currentLMax, ~, currentSequence] = latenessMax(permutationMatrix);
    else
        currentLMax = lMaxSequence;
        currentSequence = sequence;
    end
end

%% FUNZIONE: ALGORITMO LATENESS MASSIMA
% @input startingSequence, sequenza di partenza da cui comincia l'esecuzione dell'algoritmo.
% @output lMax, lateness massima migliore (più piccola) raggiunta dall'algoritmo.
% @output sequence, sequenza migliore tra quelle testate, cioè con lateness massima lMax.
function [lMax, sequence] = algorithmLatenessMax(startingSequence)
    lMax = latenessMax(startingSequence);
    sequence = startingSequence;
    [currentSequence, currentLMax] = singleIterationLatenessMax(startingSequence);  %eseguo la prima iterazione cominciando dalla sequenza di partenza. 
    while currentLMax < lMax   %continuo ad iterare solo se la lateness massima continua ad essere migliore.
        sequence = currentSequence;
        lMax = currentLMax;
        [currentSequence, currentLMax] = singleIterationLatenessMax(sequence);
    end
end

%% FUNZIONE: SINGOLA ITERAZIONE SETUP MASSIMO
% @input sequence, vettore rappresentante una certa sequenza di job.
% @output currentSequence, sequenza migliore (con minore lateness max) tra quelle ottenute con permutazioni attorno a setup max.
% @output currentLMax, lateness massima relativa a currentSequence.
function [currentSequence, currentLMax] = singleIterationSetUpMax(sequence)
    global n s sinit
    %Trovo la posizione del job con il setup maggiore.
    setUpVector = zeros(n,1);
    for i = 1:n
        if i == 1
            setUpVector(i) = sinit(sequence(i));
        else 
            setUpVector(i) = s(sequence(i-1), sequence(i));
        end
    end
    [~, indexMaxSetUp] = max(setUpVector);
    %PERMUTAZIONI ATTORNO A SETUP MAX
    permutationMatrix = [];
    %Permutazione dei due job precedenti il setup max (p1).
    if indexMaxSetUp ~= 1 && indexMaxSetUp ~= 2
        permutation1 = sequence;   
        permutation1(indexMaxSetUp-2:indexMaxSetUp-1) = [permutation1(indexMaxSetUp-1), permutation1(indexMaxSetUp-2)];
        permutationMatrix = [permutationMatrix; permutation1];
    end
    %Permutazione dei due job adiacenti al setup max (p2).
    if indexMaxSetUp ~= 1
        permutation2 = sequence;   
        permutation2(indexMaxSetUp-1:indexMaxSetUp) = [permutation2(indexMaxSetUp), permutation2(indexMaxSetUp-1)];
        permutationMatrix = [permutationMatrix; permutation2];
    end
    %Permutazione dei due job seguenti il setup max (p3).
    if indexMaxSetUp ~= n
        permutation3 = sequence;   
        permutation3(indexMaxSetUp:indexMaxSetUp+1) = [permutation3(indexMaxSetUp+1), permutation3(indexMaxSetUp)];
        permutationMatrix = [permutationMatrix; permutation3];
    end 
    [currentLMax, ~, currentSequence] = latenessMax(permutationMatrix);
end

%% FUNZIONE: ALGORITMO SETUP MASSIMO
% @input startingSequence, sequenza di partenza da cui comincia l'esecuzione dell'algoritmo.
% @output lMax, lateness massima migliore (più piccola) raggiunta dall'algoritmo.
% @output sequence, sequenza migliore tra quelle testate, cioè con lateness massima lMax. 
function [lMax, sequence] = algorithmSetUpMax(startingSequence)
    lMax = latenessMax(startingSequence);
    sequence = startingSequence;
    [currentSequence, currentLMax] = singleIterationSetUpMax(startingSequence);  %eseguo la prima iterazione cominciando dalla sequenza di partenza. 
    while currentLMax < lMax   %continuo ad iterare solo se la lateness massima continua ad essere migliore.
        sequence = currentSequence;
        lMax = currentLMax;
        [currentSequence, currentLMax] = singleIterationSetUpMax(sequence);
    end
end

%% FUNZIONE: SINGOLA ITERAZIONE RANDOM
% @input sequence, vettore rappresentante una certa sequenza di job.
% @input numIterations, numero di nuove sequenze generate con scambi casuali. 
% @output currentSequence, sequenza migliore (con minore lateness max) tra le numIterations ottenute con scambi random di 2 elementi della sequence.
% @output currentLMax, lateness massima relativa a currentSequence.
function [currentSequence, currentLMax] = singleIterationRand(sequence, numIterations)
    global n
    switchesMatrix =[];
    for i = 1:numIterations
        switchSequence = sequence;
        a = randi([1,n]);
        b = randi([1,n-1]); %evito di generare a=b
        if b >= a
            b = b+1;
        end
        t = switchSequence(a);
        switchSequence(a) = switchSequence(b);
        switchSequence(b) = t;
        switchesMatrix = [switchesMatrix; switchSequence];
    end
    [currentLMax, ~, currentSequence] = latenessMax(switchesMatrix);
end

%% FUNZIONE: ALGORITMO RANDOM
% @input startingSequence, sequenza di partenza da cui comincia l'esecuzione dell'algoritmo.
% @output lMax, lateness massima migliore (più piccola) raggiunta dall'algoritmo.
% @output sequence, sequenza migliore tra quelle testate, cioè con lateness massima lMax. 
function [lMax, sequence] = algorithmRand(startingSequence, numIterations)
    lMax = latenessMax(startingSequence);
    sequence = startingSequence;
    [currentSequence, currentLMax] = singleIterationRand(startingSequence, numIterations);  %eseguo la prima iterazione cominciando dalla sequenza di partenza. 
    while currentLMax < lMax   %continuo ad iterare solo se la lateness massima continua ad essere migliore.
        sequence = currentSequence;
        lMax = currentLMax;
        [currentSequence, currentLMax] = singleIterationRand(sequence, numIterations);
    end
end

%% FUNZIONE: LATENESS MASSIMA
% @input K, matrice con sequenze di job (sulle righe) di cui voglio calcolare la migliore lateness max.
% @output lMax, lateness massima migliore tra le lMax delle singole sequenze.
% @output index, indice della posizione nella sequenza del job con lateness massima. 
% @output sequence, sequenza di K con lateness massima minima.
function [lMax, index, sequence] = latenessMax(K)
    global DueDate n s sinit p
    vectorLMax = zeros(size(K,1), 1);
    I = zeros(size(K,1), 1);
    for j = 1:size(K, 1)   %iterazione sulle sequenze.
        C = zeros(n, 1); 
        for q = 1:n   %iterazione sulle posizioni dei job.
            if q == 1
                C(K(j, q)) = sinit(K(j, q)) + p(K(j, q));   %calcolo del tempo di completamento del primo job.
            else
                C(K(j, q)) = C(K(j, q-1)) + s(K(j, q-1), K(j, q)) + p(K(j, q));   %calcolo del tempo di completamento dei job successivi.
            end
        end
        L = C - DueDate;  %vettore delle lateness di ogni job.
        [vectorLMax(j), I(j)] = max(L);
    end
    [lMax, indexSequence] = min(vectorLMax);
    sequence = K(indexSequence,:);
    jobLMax = I(indexSequence);
    index = find(sequence == jobLMax);
end

%% FUNZIONE: SEQUENZA RANDOM
% @output rs, vettore riga contenente sequenza casuale di n jobs.
function rs = randomSequence()
    global n
    t = randi([0,100], 1, n);
    [~,rs] = sort(t);
end

%% FUNZIONE: STAMPA RISULTATI FORMATTATA (euristica iterativa)
% @input a, sequenza di arrivo.
% @input la, lateness max sequenza di arrivo.
% @input lp, lateness max sequenza di partenza.
% @input t, tempo di esecuzione dell'algoritmo.
% @output percentageImprovement, miglioramento percentuale della lateness max rispetto alla lateness max di partenza.
function percentageImprovement = formatPrintI(a, la, lp, t)
    global n
    fprintf("La sequenza ottenuta è: ");
    for i = 1:n
        fprintf("\t %d", a(i));
    end
    fprintf("\n")
    fprintf("La lateness max ottenuta è: %d.\n", la)
    %Valuto il miglioramento sulla lateness massima.
    fprintf("Il miglioramento sulla lateness max è: %d.\n", la - lp)
    percentageImprovement = round((lp - la)*100/lp, 2);
    fprintf("Il miglioramento percentuale sulla lateness max è: %s per cento.\n", num2str(percentageImprovement))
    %Stampo il tempo di esecuzione dell'algoritmo.
    fprintf("Il tempo di esecuzione dell'algoritmo è: %f secondi.\n", t)
    fprintf("\n")
    fprintf("\n")
end

%% FUNZIONE: STAMPA RISULTATI FORMATTATA
% @input s, sequenza.
% @input l, lateness max.
% @input t, tempo di esecuzione dell'algoritmo.
function formatPrint(s, l, t)
    global n
    fprintf("La sequenza ottenuta è: ");
    for i = 1:n
        fprintf("\t %.0f", s(i));
    end
    fprintf("\n")
    fprintf("La lateness max ottenuta è: %.0f.\n", l)
    %Stampo il tempo di esecuzione dell'algoritmo.
    fprintf("Il tempo di esecuzione dell'algoritmo è: %f secondi.\n", t)
    fprintf("\n")
    fprintf("\n")
end
