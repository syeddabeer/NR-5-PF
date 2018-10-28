% Matlab Code Program for Newton-Raphson Load Flow Analysis..
clear all
nbus = 5;                  % 5 Buses Given

disp('########## Calculation of Bus Admittance Matrix ###########')
%%%%%%%% Calculation of Bus Admittance Matrix from given R and X %%%%%%%
linedata = inputLineData(5);      % Calling inputLineData...
fb = linedata(:,1);             % From bus number...
tb = linedata(:,2);             % To bus number...
r = linedata(:,3);              % Resistance, R...
x = linedata(:,4);              % Reactance, X...
b = linedata(:,5);              % Ground Admittance, B/2...
a = linedata(:,6);              % Tap setting value..
z = r + i*x;                    % z matrix...
y = 1./z;                       % To get inverse of each element...
b = i*b;                        % Make B imaginary...
SC=inputBusData(5);
sc=SC(:,11);
nb = max(max(fb),max(tb));      % No. of buses...
nl = length(fb);                % No. of branches...
Y = zeros(nb,nb);               % Initialise YBus...
final=0
 
 % Formation of the Off Diagonal Elements...
 for k = 1:nl
     Y(fb(k),tb(k)) = Y(fb(k),tb(k)) - y(k)/a(k);
     Y(tb(k),fb(k)) = Y(fb(k),tb(k));
 end
 
 % Formation of Diagonal Elements....
 for m = 1:nb
     for n = 1:nl
         if fb(n) == m
             Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
         elseif tb(n) == m
             Y(m,m) = Y(m,m) + y(n) + b(n);
         end
     end
 end

 for m = 1:nb
    Y(m,m)=Y(m,m)+sc(m)*(1i); 
     
 end 

display('the Admittance Bus Ybus Matrix is calculated as')
display(Y)
display('##########################################################')


busd = inputBusData(nbus);      % Calling inputBusData..
BMva = 1;                   % We have inputted all the p.u. values
bus = busd(:,1);            % Bus Number..
type = busd(:,2);           % Type of Bus 1-Slack/Swing, 2-PV, 3-PQ..
V = busd(:,3);              % Specified Voltage..
del = busd(:,4);            % Voltage Angle..
Pg = busd(:,5)/BMva;        % PGi..
Qg = busd(:,6)/BMva;        % QGi..
Pl = busd(:,7)/BMva;        % PLi..
Ql = busd(:,8)/BMva;        % QLi..
Qmin = busd(:,9)/BMva;      % Minimum Reactive Power Limit..
Qmax = busd(:,10)/BMva;     % Maximum Reactive Power Limit..
P = Pg - Pl;                % Pi = PGi - PLi..
Q = Qg - Ql;                % Qi = QGi - QLi..
Psp = P;                    % P Specified..
Qsp = Q;                    % Q Specified..
G = real(Y);                % Conductance matrix..
B = imag(Y);                % Susceptance matrix..
pv = find(type == 2 | type == 1);   % PV Buses..
pq = find(type == 3);               % PQ Buses..
npv = length(pv);                   % No. of PV buses..
npq = length(pq);                   % No. of PQ buses..



Tol = 1;
Iter = 1;
while (Tol > 0.1)   % Convergence value given is 0.1 MVA   
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    % Calculation of P and Q
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end
   
    % Calculation of Mismatch vector for P and Q
    dPa = Psp-P;
    dQa = Qsp-Q;
    k = 1;
    dQ = zeros(npq,1);
    for i = 1:nbus
        if type(i) == 3
            dQ(k,1) = dQa(i);
            k = k+1;
        end
    end
    dP = dPa(2:nbus);
    M = [dP; dQ];       % Mismatch Vector
    
    % Jacobian
    % Jacobn1 - Derivative of Real Power Injected w.r.t. Angles..
    Jacobn1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    Jacobn1(i,k) = Jacobn1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                Jacobn1(i,k) = Jacobn1(i,k) - V(m)^2*B(m,m);
            else
                Jacobn1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % Jacobn2 - Derivative of Real Power Injected w.r.t. Voltage
    Jacobn2 = zeros(nbus-1,npq);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    Jacobn2(i,k) = Jacobn2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                Jacobn2(i,k) = Jacobn2(i,k) + V(m)*G(m,m);
            else
                Jacobn2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % Jacobn3 - Derivative of Reactive Power Injected w.r.t. Angles..
    Jacobn3 = zeros(npq,nbus-1);
    for i = 1:npq
        m = pq(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    Jacobn3(i,k) = Jacobn3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                Jacobn3(i,k) = Jacobn3(i,k) - V(m)^2*G(m,m);
            else
                Jacobn3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % Jacobn4 - Derivative of Reactive Power Injected w.r.t. V..
    Jacobn4 = zeros(npq,npq);
    for i = 1:npq
        m = pq(i);
        for k = 1:npq
            n = pq(k);
            if n == m
                for n = 1:nbus
                    Jacobn4(i,k) = Jacobn4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                end
                Jacobn4(i,k) = Jacobn4(i,k) - V(m)*B(m,m);
            else
                Jacobn4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    Jacobian = [Jacobn1 Jacobn2; Jacobn3 Jacobn4];     % Jacobian Matrix..
    X = inv(Jacobian)*M;           % Correction Vector
    dTh = X(1:nbus-1);      % Change in Angle..
    dV = X(nbus:end);       % Change in Voltage..
    
    % Updating State Vectors..
    del(2:nbus) = dTh + del(2:nbus);    % Voltage Angle..
    k = 1;
    for i = 2:nbus
        if type(i) == 3
            V(i) = dV(k) + V(i);        % Voltage Magnitude..
            k = k+1;
        end
    end 
    
    Tol = max(abs(M)); % Tolerance..
        scanner1 = ['The values (in p.u.) for iteration :: ', num2str(Iter),' :: are given as belows:'];
        display(scanner1);
        disp('#########################################################################################');
        disp('-----------------------------------------------------------------------------------------');
        disp('                              Newton Raphson Loadflow Analysis ');
        disp('-----------------------------------------------------------------------------------------');
        loadflow(nbus,V,del,BMva,Y,final);  %print after each iteration
    Iter = Iter + 1;
end

disp('*****************************************************************************************');
disp('-----------------------------------------------------------------------------------------');
disp('                     Converged Values (in p.u.) after N-R method on 5-bus system ');
disp('-----------------------------------------------------------------------------------------');
loadflow(nbus,V,del,BMva,Y,final);
disp('*****************************************************************************************');

fprintf('\n');
disp('Final Voltage values (in p.u.) are as below:');
disp(' Bus     V  ');
for m = 1:5
    fprintf('%3g', m); fprintf('  %8.4f', V(m)); 
    fprintf('\n');
end
fprintf('\n');

fprintf('\n');
disp('Final Angle values (in degrees) are as below:');
disp(' Bus     Angle ');
for m = 1:5
    fprintf('%3g', m); fprintf('  %8.4f', (180/pi)*del(m)); 
    fprintf('\n');
end
fprintf('\n');

final = 1;
disp('*****************************************************************************************');
disp('-----------------------------------------------------------------------------------------');
disp('          Final Power Values (in M.W. and MVAR) after N-R method on 5-bus system         ');
disp('-----------------------------------------------------------------------------------------');
loadflow(nbus,V,del,BMva,Y,final);
disp('*****************************************************************************************');



