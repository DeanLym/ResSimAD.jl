run C:\Users\yiminliu\mrst-2020a\startup
mrstModule add ad-core ad-blackoil ad-props mrst-gui deckformat
%% Read Eclipse input deck
inputfile = '../Eclipse/EXAMPLE1.DATA';
deck = readEclipseDeck(inputfile);
deck = convertDeckUnits(deck);
%% Init grid and rock
G = initEclipseGrid(deck);
G = computeGeometry(G);
rock = initEclipseRock(deck);
%% Init fluid
fluid = initDeckADIFluid(deck);
%% Init wells
W = processWells(G, rock, deck.SCHEDULE.control(1));
%% Set schedule
simTime = 1825*day;
nstep   = 64;
refine  = 5;
startSteps = repmat((simTime/(nstep + 1))/refine, refine, 1);
restSteps =  repmat(simTime/(nstep + 1), nstep, 1);
timesteps = [0.01*day; startSteps; restSteps];
schedule = simpleSchedule(timesteps, 'W', W);
%% Create mode
model = TwoPhaseOilWaterModel(G, rock, fluid);
%% Set initial state
p0 = 413.6854 * ones(450,1) * barsa;
sw0 = 0.1 * ones(450,1);
sat = [sw0, 1 - sw0];
state0 = initResSol(G, p0, sat);
%% Create solver
solver = BackslashSolverAD();
%% Begin simulation
[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule, 'LinearSolver', solver);
%% Get production rates
T = convertTo(cumsum(schedule.step.val), day);
[qw, qo, qg, bhp] = wellSolToVector(wellSols);
%% Save simulation results to file
dlmwrite('qo.txt', qo/(meter^3)*day*6.28981);
dlmwrite('qw.txt', qw/(meter^3)*day*6.28981);
dlmwrite('t.txt', T);
