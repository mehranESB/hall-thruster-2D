% create enviroment with one chamber 
env = Enviroment();
env.Thruster.Cathode.IsInMiddle = 1;
env.Thruster.Cathode.Height = 0.05;

chamber = Chamber;
env.Thruster.Chambers = {chamber};

I = [1, 2, 5, 10, 20, 30, 50, 70, 100, 150, 200];
D = zeros(1, length(I));
n = zeros(1, length(I));
for i = 1 : length(I)
    % modify current of chamber 
    chamber.MagnetPower = I(i);
    env.Thruster.Chambers = {chamber};

    % solve magneto static and electro static for this enviroment
    solverB = SolverB(env);
    solverE = SolverE(env);
    solverB = solve(solverB);
    solverE = solve(solverE);

    % compute efficiency regard to SPT-100 conditions 
    EC = EfficiencyCalculator(solverE, solverB);
    D(i) = EC.getDiffusionBohmCoef(10);
    n(i) = EC.getGasDensityUpperBound(0.2, 10000);
end
L = 170 * I;




