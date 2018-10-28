% Input given for Buses Data. 
% Unknown Power is assigned 0 value
% Unknown Voltage is assigned value of 1
% Unknown Angle is assigned value of 0
% Type of Bus....
% 1 - Slack/Swing Bus. 2 - PV Bus. 3 - PQ Bus.
function busdt = inputBusData(num)
%         |Bus | Type | Vsp | theta | 	PGi |   QGi |  PLi |     QLi |     Qmin |            Qmax |    Shunt Capacitance
busdt =    [1     1    1.000    0      0       0        0        0          0                0          0;
            2     3    1.000    0      0       0        8.00     2.80       0                0          0; 
            3     2    1.050    0      5.20    0        0.80     0.40      -2.8             4.0         0;
            4     3    1.000    0      0       0        0        0          0               0           0;
            5     3    1.000    0      0       0        0        0          0               0           0;
           ];
end