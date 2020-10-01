/*
 * Scilab environment settings
 */

clear;
clc;

/*
 * Scilab custom matrix functions
 */

function Matrix = addZerosColumn(Matrix)
  Matrix = [Matrix zeros(size(Matrix, 1),1)];
endfunction

function Matrix = addConstantColumn(Matrix, constant)
  Matrix = addZerosColumn(Matrix);
  Matrix(:,$) = constant;
endfunction

function vector = fixInf(vector)
  for index=1:size(vector,1)
    element = vector(index,1);
    if element == %inf then
      vector(index,1) = 0; 
    end
  end
endfunction

/*
 * Electromagnetic theory functions
 */

function Zb = getBaseImpedance(Vb, Sb)
  Zb = Vb.^2./Sb;
endfunction

function Ib = getBaseCurrent(Vb, Sb)
  Ib = Sb./(Vb.*sqrt(3));
endfunction

function valueInPU = value2PU(value, base)
  valueInPU = value./base;
endfunction

function value = PU2value(PU, base)
  value = PU.*base;
endfunction

function w = getAngularVelocity(f)
  w = 2.*%pi.*f;
endfunction

function X = getReactance(w, L, C)
  Xl = w.*L
  Xc = zeros(size(w, 1),1);
  Xc = 1./(w.*C);
  Xc = fixInf(Xc);
  X = Xl-Xc;
endfunction

function Y = getAdmittance(Z)
  Y = 1./Z;
  Y = fixInf(Y);
endfunction

function G = getConductance(R, X)
  G = R./(R.^2+X.^2);
endfunction

function B = getSusceptance(R, X)
  B = (-1).*X./(R.^2+X.^2);
endfunction

/*
 * Scilab math custom functions
 */
 
function Angle = angle(complex)
  Angle = atan(imag(complex),real(complex));
endfunction

/*
 * Power flow functions
 */
 
function parallel = isParallel(nTLs, startBus, endBus)
  parallel = zeros(nTLs,1);
  for TL=2:nTLs
    if startBus(TL-1)==startBus(TL) & endBus(TL-1)==endBus(TL) then
      parallel(TL-1,1) = 1;
      parallel(TL,1) = 1;
    end
  end  
endfunction

function Ybus = ybus(nBus, nTLs, startBus, endBus, Yseries, Yshunt)
  parallel = isParallel(nTLs, startBus, endBus);
  Ybus = zeros(nBus,nBus);
  for i=1:nBus
    for TL=1:nTLs
      diagonalPrincipal = ((startBus(TL)==i)||(endBus(TL)== i));
      if diagonalPrincipal then
        Ybus(i,i) = Ybus(i,i) + Yseries(TL) + Yshunt(TL);
      else
        if parallel(TL,1) then
          Ybus(startBus(TL),endBus(TL)) = -2*Yseries(TL);
          Ybus(endBus(TL),startBus(TL)) = -2*Yseries(TL);
        else
          Ybus(startBus(TL),endBus(TL)) = -Yseries(TL);
          Ybus(endBus(TL),startBus(TL)) = -Yseries(TL);
        end
      end
    end
  end
endfunction

function [iteration, tolerance, V] = voltageByGaussSeidel(nbus, typeBus, V, P, Q, Y)
  tolerance = 1;
  toleranceMax = 1e-4;
  iteration = 0;
  iterationMax = 1e4;
  Vprev = V;
  accelerationFactor = 1.5;
  while ( (tolerance>toleranceMax) & (iteration<iterationMax) )
    for i = 2:nbus
      // Sum Yik * Vk
      sumYV = 0;
      for k = 1:nbus
        if i ~= k then
          sumYV = sumYV + Y(i,k)* V(k);
        end
      end
      // Computing Qi for PV bus
      if typeBus(i) == 2 then
        Q(i) = -imag(conj(V(i))*(sumyv + Y(i,i)*V(i)));
        // Checking for Qi Violation.
        if (Q(i) > Qmax(i)) | (Q(i) < Qmin(i)) then
          if Q(i) < Qmin(i) then
            Q(i) = Qmin(i);
          else
            Q(i) = Qmax(i);
          end
          typeBus(i) = 1;
        end
      end
      V(i) = (1/Y(i,i))*((P(i)-%i*Q(i))/conj(V(i))-sumYV);
      // For PV Buses, Voltage Magnitude remains same, but Angle changes
      if typeBus(i) == 2 then
        V(i) = abs(Vprev(k))*cos(angle(V(k))) + %i*abs(Vprev(k))*sin(angle(V(k)));
      end
    end
    V = Vprev+accelerationFactor.*(V-Vprev);
    tolerance = max(abs(abs(V)-abs(Vprev)));
    iteration = iteration + 1;
    Vprev = V;
  end
endfunction

function [iteration, tolerance, V] = voltageByNewtonRaphson(nbus, typeBus, V, P, Q, Y)
  tolerance = 1;
  toleranceMax = 1e-6;
  iteration = 0;
  iterationMax = 1e6;
  theta = zeros(nbus,1);
  Pesp = P; 
  Qesp = Q;
  G = real(Y);
  B = imag(Y);
  pq = find(typeBus==1);
  nPQ = size(pq,2);  
  while ((tolerance > toleranceMax) && (iteration<iterationMax))
    P = zeros(nbus,1);
    Q = zeros(nbus,1);
    //Calc P e Q
    for i = 1:nbus    
      for k = 1:nbus
        P(i)= P(i) + V(i)*V(k)*(G(i,k)*cos(theta(i)-theta(k)) + B(i,k)*sin(theta(i)-theta(k)));
        Q(i)= Q(i) + V(i)*V(k)*(G(i,k)*sin(theta(i)-theta(k)) - B(i,k)*cos(theta(i)-theta(k)));
      end
    end
    // Calc delta P e Q
    dPaux = Pesp-P;
    dQaux = Qesp-Q;
    k = 1;
    dQ = zeros(nPQ, 1);
    for i = 1:nbus
      if typeBus(i) == 1    
        dQ(k,1) = dQaux(i);
        k = k+1;
      end
    end
    //Delta P
    dP = dPaux(2:nbus);
    // Delta Power Matrix
    F = [dP;
         dQ
        ];
    /*
     * Jacobiana Matrix
     */
    // H Matrix - Active power partial derivative by angles
    H = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
      m = i+1;
      for k = 1:(nbus-1)
        n = k + 1;
        if n == m
          for n = 1:nbus
            H(i,k) = H(i,k) + V(m)*V(n)*(-G(m,n)*sin(theta(m)-theta(n)) + B(m,n)*cos(theta(m)-theta(n)));
          end
          H(i,k) = H(i,k) - V(m)^2*B(m,m);
        else
          H(i,k) = V(m)*V(n)*(G(m,n)*sin(theta(m)-theta(n)) - B(m,n)*cos(theta(m)-theta(n)));
        end
      end
    end
    // N Matrix - Active power partial derivative by voltage
    N = zeros(nbus-1,nPQ);
    for i = 1:(nbus-1)
      m = i+1;
      for k = 1:nPQ
        n = pq(k);
        if n == m
          for n = 1:nbus
            N(i,k) = N(i,k) + V(n)*(G(m,n)*cos(theta(m)-theta(n)) + B(m,n)*sin(theta(m)-theta(n)));
          end
          N(i,k) = N(i,k) + V(m)*G(m,m);
        else
          N(i,k) = V(m)*(G(m,n)*cos(theta(m)-theta(n)) + B(m,n)*sin(theta(m)-theta(n)));
        end
      end
    end
    // M Matrix - Reactive power partial derivative by angles
    M = zeros(nPQ,nbus-1);
    for i = 1:nPQ
      n = pq(i);
      for k = 1:(nbus-1)
        n = k+1;
        if n == m
          for n = 1:nbus
            M(i,k) = M(i,k) + V(m)*V(n)*(G(m,n)*cos(theta(m)-theta(n)) + B(m,n)*sin(theta(m)-theta(n)));
          end
          M(i,k) = M(i,k) - V(m)^2*G(m,m);
        else
          M(i,k) = V(m)*V(n)*(-G(m,n)*cos(theta(m)-theta(n)) - B(m,n)*sin(theta(m)-theta(n)));
        end
      end
    end
    // L Matrix - Reactive power partial derivative by voltage
    L = zeros(nPQ,nPQ);
    for i = 1:nPQ
      m = pq(i);
      for k = 1:nPQ
        n = pq(k);
        if n == m
          for n = 1:nbus
            L(i,k) = L(i,k) + V(n)*(G(m,n)*sin(theta(m)-theta(n)) - B(m,n)*cos(theta(m)-theta(n)));
          end
          L(i,k) = L(i,k) - V(m)*B(m,m);
        else
          L(i,k) = V(m)*(G(m,n)*sin(theta(m)-theta(n)) - B(m,n)*cos(theta(m)-theta(n)));
        end
      end
    end  
    // Jacobiana Matrix
    J = [ H N;
          M L ];
    //X = inv(J)*F
    X = J\F;
    dTheta = X(1:nbus-1);
    dV = X(nbus:$);
    // Update status vector
    theta(2:nbus) = dTheta + theta(2:nbus);
    k = 1;
    for i = 2:nbus
      if typeBus(i) == 1
        V(i) = dV(k) + V(i);
        k = k+1;
      end
    end
    tolerance = max(abs(F));
    iteration = iteration + 1;
  end
endfunction

function [iteration, tolerance, V] = fastDecoupledNewtonRaphson(bus, TLs, V)
endfunction

/*
 * TLs: get and add power flow input parameters
 */

function TLs = addBaseImpedanceColumn(TLs)
  TLs = addZerosColumn(TLs);
  Vb = TLs(:,7);
  Sb = TLs(:,12);
  Zb = getBaseImpedance(Vb, Sb);
  TLs(:,$) = Zb;
endfunction

function TLs = addBaseCurrentColumn(TLs)
  TLs = addZerosColumn(TLs);
  Vb = TLs(:,7);
  Sb = TLs(:,12);
  Ib = getBaseCurrent(Vb, Sb);
  TLs(:,$) = Ib;
endfunction

function TLs = addAngularVelocityColumn(TLs)
  TLs = addZerosColumn(TLs);
  f = TLs(:,10);
  w = getAngularVelocity(f);
  TLs(:,$) = w;
endfunction

function TLs = addSerieAdmittanceInPUColumn(TLs)
  TLs = addZerosColumn(TLs);
  R = TLs(:,4);
  X = TLs(:,5);
  Zb = TLs(:,13);
  Rpu = value2PU(R,Zb);
  Xpu = value2PU(X,Zb);
  Gpu = getConductance(Rpu, Xpu);
  Bpu = getSusceptance(Rpu, Xpu);
  Ypu = Gpu+%i.*Bpu;
  TLs(:,$) = Ypu;
endfunction

function TLs = addShuntAdmittanceInPUColumn(TLs)
  TLs = addZerosColumn(TLs);
  w = TLs(:,11);
  C = TLs(:,6);
  L = zeros(size(TLs,1),1);
  X = getReactance(w,L,C);
  Zs = %i*X 
  Zb = TLs(:,13);
  Zspu = value2PU(Zs,Zb);
  Yspu = getAdmittance(Zspu);
  TLs(:,$) = Yspu./2;
endfunction   

/*
 * TLs: get and add power flow output parameters
 */

function TLs = addSerieCurrentInPUColumn(TLs, V, Ybus)
  TLs = addZerosColumn(TLs);
  for TL = 1:size(TLs, 1)
    startBus = TLs(TL,2);
    endBus = TLs(TL,3);
    i = startBus;
    j = endBus;
    I = (V(i)-V(j))*Ybus(i,j);
    TLs(TL,$) = I;
  end
endfunction

function TLs = addShuntCurrentInPUColumn(TLs, V)
  TLs = addZerosColumn(TLs);
  for TL = 1:size(TLs, 1)
    bus = TLs(TL,2);
    i = bus;
    Ys = TLs(TL,16)/2;
    Is = V(i)*Ys;
    TLs(TL,$) = Is;
  end
endfunction

function TLs = addNetCurrentInPUColumn(TLs)
  nTLs = size(TLs,1);
  startBus = TLs(:,2);
  endBus = TLs(:,3);
  parallel = isParallel(nTLs, startBus, endBus);
  TLs = addZerosColumn(TLs);
  I = TLs(:,17);
  Is = TLs(:,18);
  for TL = 1:nTLs
    Itl(TL,1) = I(TL,1)+Is(TL,1);
    if parallel(TL,1) then
      Itl(TL,1) = Itl(TL,1)/2;
    end
  end
  TLs(:,$) = Itl;
endfunction

function TLs = addActivePowerOnTLInPU(TLs, V)
  TLs = addZerosColumn(TLs);
  for TL = 1:size(TLs, 1)
    bus = TLs(TL,2);
    i = bus;
    Itl = TLs(TL,19);
    Stl = V(i)*conj(Itl);
    Ptl = real(Stl);
    TLs(TL,$) = Ptl;
  end
endfunction

function TLs = addReactivePowerOnTLInPU(TLs, V)
  TLs = addZerosColumn(TLs);
  for TL = 1:size(TLs, 1)
    bus = TLs(TL,2);
    i = bus;
    Itl = TLs(TL,19);
    Stl = V(i)*conj(Itl);
    Qtl = imag(Stl);
    TLs(TL,$) = Qtl;
  end
endfunction

/*
 * Bus: add flow parameters
 */

function bus = updateVoltageInPUColumn(bus, V)
  bus(:,6) = abs(V);
  bus(:,7) = (180/%pi)*angle(V);
endfunction

function bus = addNetActivePowerInPUColumn(bus)
  bus = addZerosColumn(bus);
  Pg = bus(:,2);
  Pc = bus(:,4);
  Pl = Pg-Pc;
  Sb = bus(:,9);
  Plpu = value2PU(Pl,Sb);
  bus(:,$) = Plpu;
endfunction

function bus = addNetReactivePowerInPUColumn(bus)
  bus = addZerosColumn(bus);
  Qg = bus(:,3);
  Qc = bus(:,5);
  Ql = Qg-Qc;
  Sb = bus(:,9);
  Qlpu = value2PU(Ql,Sb);
  bus(:,$) = Qlpu;
endfunction

function bus = addCurrentInPUColumn(bus, V, Y)
  bus = addZerosColumn(bus);
  for k = 1:size(Y, 1)
    I = 0;
    for N = 1:size(Y, 1)
      I = I + Y(k, N)*V(N); //I = (Y11V1+Y12V2+Y13V3+...+Y1nVn)
    end
    bus(k,$) = I;
  end
endfunction

function bus = addActivePowerInPUColumn(bus, V)
  bus = addZerosColumn(bus);
  I = bus(:,12);
  S = V.*conj(I);
  P = real(S)+bus(:,4);
  bus(:,$) = P;
endfunction

function bus = addReactivePowerInPUColumn(bus, V)
  bus = addZerosColumn(bus);
  I = bus(:,12);
  S = V.*conj(I);
  Q = imag(S)+bus(:,5);
  bus(:,$) = Q;
endfunction

/*
 * Transmission Lines Parameters
 */

//                                               Series Impedance          Shunt Impedance     base                   Tap
//      Nº      startBus      endBus        resistance       reactance       capacitance      voltage      Imax      Trans.
TLs = [ 01         01           02            16.9587         24.9762         4.1690e-7         69e3       326.0     0.000;
        02         01           03             5.1138         12.4038         2.5360e-7         69e3       435.0     0.000;
        03         01           03             5.1138         12.4038         2.5360e-7         69e3       435.0     0.000;
        04         01           04             6.5702         28.8326         6.3310e-7         69e3       636.0     0.000;
        05         01           05             0.0000         36.0408         0.0000            69e3        84.0     1.000;
        06         01           05             0.0000         36.0408         0.0000            69e3        84.0     1.000;
        07         01           11            10.4561         25.3609         5.1830e-7         69e3       435.0     0.000;
        08         02           04             8.0937         11.9215         1.9900e-7         69e3       326.0     0.000;
        09         02           06             0.0000         61.3217         0.0000            69e3        42.0     1.000;
        10         02           07             0.0000         35.7075         0.0000            69e3        84.0     1.000;
        11         02           12             2.6233          6.3654         1.2990e-7         69e3       435.0     0.000;
        12         03           08             0.0000         35.0886         0.0000            69e3        84.0     1.000;
        13         03           08             0.0000         35.0886         0.0000            69e3        84.0     1.000;
        14         03           10             2.2043          5.3514         1.0920e-7         69e3       175.0     0.000;
        15         04           09             0.0000         36.1360         0.0000            69e3        90.0     1/0.975;
        16         04           09             0.0000         36.1360         0.0000            69e3        90.0     1/0.975;
        17         04           13             4.8943         12.8261         2.3240e-7         69e3       402.0     0.000;
        18         10           02             2.6804          6.5035         1.3260e-7         69e3       435.0     0.000;
        19         13           14             0.0000         41.2303         0.0000            69e3       84.0     1/0.950;
      ];
// Column 10: Frequency [Hz]
TLs = addConstantColumn(TLs, 60);
// Column 11: Angular Velocity [rad/s]
TLs = addAngularVelocityColumn(TLs);
// Column 12: Base Power [VA]
TLs = addConstantColumn(TLs, 100e6);
// Column 13: Base Impedance [Ohm]
TLs = addBaseImpedanceColumn(TLs);
// Column 14: Base Current [A]
TLs = addBaseCurrentColumn(TLs);
// Column 15: Serie Admittance (p.u)
TLs = addSerieAdmittanceInPUColumn(TLs);
// Column 16: Shunt Admittance (p.u)
TLs = addShuntAdmittanceInPUColumn(TLs);

/*
 * Bus parameters
 */

//      Nº     Pg       Qg      Pc       Qc         V      theta    type
bus = [ 01    0.00     0.00    0.00     0.00      1.043     0.00     3;
        02    0.00     0.00    0.00     0.00      1.000     0.00     1;
        03    0.00     0.00    0.00     0.00      1.000     0.00     1;
        04    0.00     0.00    0.00    -5.00e6    1.000     0.00     1;
        05    0.00     0.00   18.21e6   1.90e6    1.000     0.00     1;
        06    0.00     0.00    2.85e6  -0.50e6    1.000     0.00     1;
        07    0.00     0.00    6.62e6  -0.87e6    1.000     0.00     1;
        08    0.00     0.00   15.21e6   2.09e6    1.000     0.00     1;
        09    0.00     0.00   18.12e6   4.02e6    1.000     0.00     1;
        10    0.00     0.00   17.24e6   5.67e6    1.000     0.00     1;
        11    0.00     0.00    7.00e6   2.30e6    1.000     0.00     1;
        12    0.00     0.00    5.00e6   1.64e6    1.000     0.00     1;
        13    0.00     0.00    0.00     0.00      1.000     0.00     1;
        14    0.00     0.00    5.77e6   1.53e6    1.000     0.00     1;
      ];
// Column 09: Base Power [VA]
bus = addConstantColumn(bus, 100e6);
// Column 10: Net Active Power (p.u)
bus = addNetActivePowerInPUColumn(bus);
// Column 11: Net Reactive Power (p.u)
bus = addNetReactivePowerInPUColumn(bus);

/*
 * Admittance matrix
 */
//             nBus         nTLs      startBus   endBus    Yseries    Yshunt
Ybus = ybus(size(bus,1), size(TLs,1), TLs(:,2), TLs(:,3), TLs(:,15), TLs(:,16));

/*
 * Bus voltage
 */
Vref = bus(:,6);
//                                                 nBus        type             NetP       NetQ
[iteration, tolerance, V] = voltageByGaussSeidel(size(bus,1), bus(:,8), Vref, bus(:,10), bus(:,11), Ybus);
bus = updateVoltageInPUColumn(bus, V);
[iteration2, tolerance2, V2] = voltageByNewtonRaphson(size(bus,1), bus(:,8), Vref, bus(:,10), bus(:,11), Ybus);

/* 
 * Power flow in the system transmission line
 */

// Column 17: Serie Current (p.u)
TLs = addSerieCurrentInPUColumn(TLs, V, Ybus);
// Column 18: Shunt Current (p.u)
TLs = addShuntCurrentInPUColumn(TLs, V);
// Column 19: Net Current (p.u)
TLs = addNetCurrentInPUColumn(TLs);
// Column 20: Active Power (p.u)
TLs = addActivePowerOnTLInPU(TLs, V);
// Column 21: Reactive Power (p.u)
TLs = addReactivePowerOnTLInPU(TLs, V);

/* 
 * Power flow in the system bars
 */

// Column 12: Current (p.u)
bus = addCurrentInPUColumn(bus, V, Ybus)
// Column 13: Active Power (p.u)
bus = addActivePowerInPUColumn(bus, V);
// Column 14: Reactive Power (p.u)
bus = addReactivePowerInPUColumn(bus, V);

/* 
 * Show Power flow results
 */

printf('\n Bus Voltage: \n');
printf('\n        Methods  Iterations  Max Tolerance');
printf('\n   Gauss-Seidel          %d       %lf',iteration,tolerance);
printf('\n Newton-Raphson          %d       %lf',iteration2,tolerance2);
printf('\n');
for k=1:size(bus, 1)
  printf('\n\tV%02d= %.4f ∟ %.4fº p.u\n', k, bus(k,6), bus(k,7));
end
printf('\n Load Flow\n');
printf('\n  Slack bus:\n');
printf('\n\t  P [MW] \t Q [MVAr] \n');
printf('\n\t %.4f \t %.4f \n\n',bus(1,13)*100,bus(1,14)*100);
printf('  Transmission Lines:\n');
printf('\n\t LINE   \t  P [MW] \t Q [MVAr] \t  I [A] \t LOAD\n');
for TL = 1:size(TLs, 1)
  printf('\n\t %02d-%02d  \t %.4f \t %.4f \t %.4f \t %.1f%%\n',TLs(TL,2), TLs(TL,3), TLs(TL,20)*100, TLs(TL,21)*100,PU2value(abs(TLs(TL,19)),TLs(TL,14)), PU2value(abs(TLs(TL,19)),TLs(TL,14))/TLs(TL,8)*100);  
end
