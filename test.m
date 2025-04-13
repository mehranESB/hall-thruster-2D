
env = Enviroment();
env.Thruster.Cathode.IsInMiddle = 1;
env.Thruster.Cathode.Height = 0.10;
env.Thruster.Chambers{1}.AnodeVoltage = 100;
env.Thruster.Chambers{2}.AnodeVoltage = 100;

solverB = SolverB(env);
solverE = SolverE(env);

solverB = solve(solverB);
solverE = solve(solverE);

EC = EfficiencyCalculator(solverE, solverB);
D = EC.getDiffusionBohmCoef(10);
n = EC.getGasDensityUpperBound(0.2, 10000);



