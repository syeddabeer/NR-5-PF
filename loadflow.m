% Program for Bus Power Injections, Line & Power flows (p.u)...
function [Pinj,Qinj,Pgen,Qgen,Pload,Qload,Y] = loadflow(nb,V,del,BMva,Y,final)
lined = inputLineData(nb);          % Get linedats..
busd = inputBusData(nb);            % Get busdatas..
Vm = V.*cos(del) + j*V.*sin(del);   % Converting polar to rectangular..
Del = 180/pi*del;               % Bus Voltage Angles in Degree...
fb = lined(:,1);                % From bus number...
tb = lined(:,2);                % To bus number...
nl = length(fb);                % No. of Branches..
Pload = busd(:,7);                 % PLi..
Qload = busd(:,8);                 % QLi..
Iij = zeros(nb,nb);
Sij = zeros(nb,nb);
Si = zeros(nb,1);

% Bus Current Injections..
 I = Y*Vm;
 Im = abs(I);
 Ia = angle(I);
 
%Line Current Flows..
for m = 1:nl
    p = fb(m); q = tb(m);
    Iij(p,q) = -(Vm(p) - Vm(q))*Y(p,q); % Y(m,n) = -y(m,n)..
    Iij(q,p) = -Iij(p,q);
end
Iij = sparse(Iij);
Iijm = abs(Iij);
Iija = angle(Iij);

% Line Power Flows..
for m = 1:nb
    for n = 1:nb
        if m ~= n
            Sij(m,n) = Vm(m)*conj(Iij(m,n))*BMva;
        end
    end
end
Sij = sparse(Sij);
Pij = real(Sij);
Qij = imag(Sij);
 
% Line Losses..
Lij = zeros(nl,1);
for m = 1:nl
    p = fb(m); q = tb(m);
    Lij(m) = Sij(p,q) + Sij(q,p);
end
Lpij = real(Lij);
Lqij = imag(Lij);

% Bus Power Injections..
for i = 1:nb
    for k = 1:nb
        Si(i) = Si(i) + conj(Vm(i))* Vm(k)*Y(i,k)*BMva;
    end
end
Pinj = real(Si);
Qinj = -imag(Si);
Pgen = Pinj+Pload;
Qgen = Qinj+Qload;
 
if final == 0
disp('| Bus |    V   |  Angle  |  Generation (p.u.) |    Load (p.u.)     | Injected Power(G-L)|');
disp('| No  |   pu   |  Degree |    MW   |   MVar   |    MW   |  Mvar    |     MW     |  MVar | ');
for m = 1:nb
    disp('-----------------------------------------------------------------------------------------');
    fprintf('%3g', m); fprintf('  %8.4f', V(m)); fprintf('   %8.4f', Del(m)); 
    fprintf('  %8.3f', Pgen(m)); fprintf('   %8.3f', Qgen(m)); 
    fprintf('  %8.3f', Pload(m)); fprintf('   %8.3f', Qload(m)); 
    fprintf('  %8.3f', Pinj(m)); fprintf('   %8.3f', Qinj(m)); 
    fprintf('\n');
end
disp('-----------------------------------------------------------------------------------------');
fprintf(' Total                  '); 
fprintf('  %8.3f', sum(Pinj+Pload)); fprintf('   %8.3f', sum(Qinj+Qload));
fprintf('  %8.3f', sum(Pload)); fprintf('   %8.3f', sum(Qload)); 
fprintf('  %8.3f', sum(Pinj)); fprintf('   %8.3f', sum(Qinj));
fprintf('\n');
disp('-----------------------------------------------------------------------------------------');
disp('#########################################################################################');
end

if final == 1
disp('| Bus |   Generation     |       Load         | Injected Power(G-L)|');
disp('| No  |   MW   |   MVar  |    MW   |  Mvar    |     MW     |  MVar | ');
for m = 1:nb
    disp('-----------------------------------------------------------------------------------------');
    fprintf('%3g', m);  
    fprintf('  %8.3f', round(Pgen(m)*100)); fprintf('   %8.3f', round(Qgen(m)*100)); 
    fprintf('  %8.3f', round(Pload(m)*100)); fprintf('   %8.3f', round(Qload(m)*100)); 
    fprintf('  %8.3f', round(Pinj(m)*100)); fprintf('   %8.3f', round(Qinj(m)*100)); 
    fprintf('\n');
end
disp('-----------------------------------------------------------------------------------------');
fprintf(' Total'); 
fprintf('  %8.3f', sum(Pgen*100)); fprintf('   %8.3f', sum(Pgen*100));
fprintf('  %8.3f', sum(Pload*100)); fprintf('   %8.3f', sum(Qload*100)); 
fprintf('  %8.3f', sum(Pinj*100)); fprintf('   %8.3f', sum(Qinj*100));
fprintf('\n');
disp('-----------------------------------------------------------------------------------------');
disp('#########################################################################################');
end
    