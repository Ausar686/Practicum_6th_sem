%% TEST 1
clear, clc
A = [1 0; 0 2];
B = [1 0; 0 4];
t0 = 0;
Q = [1 0; 0 1];
p = [0;0];
a = 3;
alpha = -2;
gamma = 1;
%% TEST 2
clear, clc
A = [1 2; 1 2];
B = [1 0; 0 1];
t0 = 0;
Q = [2 0; 0 2];
p = [0;0];
a = 4;
alpha = -5;
gamma = 3;
%% TEST 3
clear, clc
A = [0 0; 0 0];
B = [1 0; 0 1];
t0 = 0;
Q = [9 0; 0 4];
p = [0;0];
a = 4;
alpha = 0;
gamma = 0.001;
%% TEST 4
clear, clc
A = [1 0; 0 1];
B = [1 0; 0 1];
t0 = 0;
Q = [1 0; 0 1];
p = [0;0];
a = 4;
alpha = 0;
gamma = 10;
%% TEST 5
clear, clc
A = [1 -3; 1 5];
B = [2 0; 3 4];
t0 = 0;
Q = [1 0; 0 1];
p = [0;0];
a = 4;
alpha = 0;
gamma = 1;
%% TEST 6
clear, clc
A = [0.01 -1; 2 0.01];
B = [1 0; 0 1];
t0 = 0;
Q = [1 0; 0 1];
p = [0;0];
a = 3;
alpha = -5.22;
gamma = 0.0001;
%% MAIN PART
% ����� ���������:
% 1) ��������� ����������, 
% 2) ���������� ���������� 
% 3) ��������� ������ (����������)
% 4) ��������� ������ (���������)
% 5) ���������� ��������
% 6) �����

% ���������� ��������
constTimeLimit = 5;
constTimePenalty = 0.1;
constEps = 0.0001;
constArrayLength = 500;
constNSteps = 30;

%������ ��������� ������
disp("Welcome to linear-OC-problem-solver!");
disp("Options available are:");
disp("1) Settings");
disp("2) Computation");
disp("3) Quality increase (global)");
disp("4) Quality increase (local)");
disp("5) Plot generation");
disp("6) Exit");
disp("NOTE: You have to complete 1) and 2) in order to do 3-5");
optionStart = input("Input your option: ");

% ��������� ���������
if optionStart == 1
    A = input("A = "); % ������� 2�2
    B = input("B = "); % ������� 2�2
    t0 = input("t0 = "); % �����
    Q = input("Q = "); % ������� 2�2
    p = input("p = "); % ������-������� 1�2
    a = input("a = "); % ����� 
    alpha = input("alpha = "); % ����� > 0
    gamma = input("gamma = "); % �����

    
% ���������� �����������
elseif optionStart == 2
    resultTime = t0 + constTimeLimit;
    nSteps = constNSteps;
    angleStep = 2 .* pi / nSteps;
    d = a^2 - 4;
    arrayX1 = [];
    arrayX2 = [];
    arrayU1 = [];
    arrayU2 = [];
    arrayPsi1 = [];
    arrayPsi2 = [];
    arrayT = [];
    
    % �������� ������������ �������
    if gamma < 0
        disp("X0 is empty. No solutions");
        return
    end
    if a^2 - 4 < 0
        disp("X1 is empty. No solutions.")
        return
    end
    if abs(det(Q)) < eps
        disp("ERROR. det(Q) = 0.");
        return
    end
    
    % ������� �������� (����� ��� ����� ������� �������� � ������ 5)
    xBorders = [-2*sqrt(gamma) - 2, max(2*sqrt(gamma), a + 2 * sqrt(d) + 1)];
    yBorders = [ min(alpha - 2*gamma, 0), max(alpha + 2*gamma, 1 /(a - sqrt(d)) + 1)];
    xLow = xBorders(1) - 10;
    xUp = xBorders(2) + 10;
    yLow = yBorders(1) - 10;
    yUp = yBorders(2) + 10;
    
    % ����������� ������������� �������� X0 � X1
    figure;
    hold on;
    x0Line = @(x, y) (x.^2 + (y - alpha).^4 <= gamma);
    x_ = linspace(-sqrt(gamma) - 1, sqrt(gamma) + 1, 1001);
    y_ = linspace(alpha - sqrt(gamma) - 1, alpha + sqrt(gamma) + 1, 1001);
    [X_, Y_] = meshgrid(x_, y_);
    Z_ = x0Line(X_, Y_);
    contourf(X_, Y_, Z_, [1 1]);
    x1Line = @(x, y) (x .* y >= 1) .* (x + y <= a) .* (x >= 0);
    x_ = linspace(0, (a + sqrt(d)) / 2 + 2, 1001);
    y_ = linspace(0, (a + sqrt(d)) / 2 + 2, 1001);
    [X_, Y_] = meshgrid(x_, y_);
    Z_ = x1Line(X_, Y_);
    contourf(X_, Y_, Z_, [1 1]);    
    hold off;
    
    % ���� ��������� ������������, ������� ����� ����� 0
    if max(max(x0Line(X_, Y_) .* x1Line(X_, Y_))) > 0 % �������� �� ����������� ������
        disp("X0 intersects X1. The time required is equal 0.");
        isIntersected = 1;
        hold off;
        return
    else
        isIntersected = 0;
    end
    
    % ������������� (��� ���������� ������ �������������)
    if rank(horzcat(A, A * B)) == 2
        isControllable = 1;
    else
        isControllable = 0;
        [V, J] = jordan(B);
        for i = 1:size(J, 1)
            if abs(J(i, i)) < eps
                J(i, i) = J(i, i) + eps;
            end
        end
        B = V * J * inv(V);
    end
    
    % �������� ����: ������� ������� ���������� ������ �� ��������� �����
    for angle = 0:angleStep:2*pi
        
        % ���������� �������������� ��������� ����� x0 (������ �������)
        psi0 = [cos(angle); sin(angle)]; % ������ psi(t0), ������-�������
        f = @(y) -(abs(psi0(1)) * sqrt(gamma - y .^ 4) + abs(psi0(2)) .* y); % ��������� ������������ �� ������� X0
        supportFunctionX0 = -(f(fminbnd(f, 0, gamma.^(1/4)))) + psi0(2) .* alpha; % �������� ������� ������� � X0 � ����������� psi(t0)
        h = @(x) abs(sum(psi0 .* x) - supportFunctionX0) + (x(1).^2 + (x(2) - alpha).^4 > gamma); % ������� ��� ������ ��������������� x0
        x0 = fminsearch(h, [0; alpha], optimset('MaxFunEvals', 1000)); % �������������� ��������� ����� � ��������� ������������ � ������ "����"
        
        % ���������� ��������������� ���������� u (������ �������)
        psiFunc = @(t) cell2mat(arrayfun(@(x) expm(-A' * (x - t0)) * psi0, t, 'UniformOutput', false)); % ������� 2xn �������� psi(t), ��� t - ������ 1xn 
        btPsi = @(t) B' * psiFunc(t); % ������� 2xn �������� B' * psi(t), ��� t - ������ 1xn
        supportFunctionP = @(t) sum(btPsi(t) .* p) + sqrt(sum(btPsi(t) .* (Q * btPsi(t)))); % ������ 1�n �������� ������� ������� � ����������� B'*psi(t) � P 
        u = @(t) p + (Q * btPsi(t)) ./ sqrt(sum(btPsi(t) .* (Q * btPsi(t))));
        %u = @(t) cell2mat(arrayfun(@(x) fminsearch(@(w) abs(sum(btPsi(x) .* w) - supportFunctionP(x)) + (sum((w-p) .* (Q * (w-p))) > 1), p), t, 'UniformOutput', false));
        
        % ������� ����������������� ���������
        tspan = [t0 resultTime]; % ��������� ��������� ����������
        options = odeset('Events', @(t, x) borderFunc(t, x, xLow, xUp, yLow, yUp)); % ��������� �������� ������ �� �������
        [t, x] = ode45(@(t, x) A * x + B * u(t), tspan, x0, options);
        
        % �������� ������� ����������� �������
        idxSolutionInX1 = find(x1Line(x(:, 1), x(:, 2)) > 0);
        inX1Size = size(idxSolutionInX1);
        if inX1Size(1) > 0
            idxFirst = idxSolutionInX1(1);
            timeCurrent = t(idxFirst);
            %flagGoodDirection = 1;
        else
            tmp = size(t);
            idxFirst = min(tmp(1), constArrayLength);
            timeCurrent = resultTime + constTimePenalty;
            %flagGoodDirection = 0;
        end
        
        % ������� "��������" ���������� � ��������� ���������� ����� ��������� � X1 � �������������� �����������
        t = t(1:idxFirst);
        x = x(1:idxFirst, :);
        t(idxFirst:constArrayLength) = t(idxFirst);
        x(idxFirst:constArrayLength, 1) = x(idxFirst, 1);
        x(idxFirst:constArrayLength, 2) = x(idxFirst, 2);
        
        % ���������� �������� �����������
        u_ = u(t');
        psi = psiFunc(t');
        arrayT = [arrayT t];
        arrayX1 = [arrayX1 x(:, 1)];
        arrayX2 = [arrayX2 x(:, 2)];
        arrayU1 = [arrayU1 u_(1, :)'];
        arrayU2 = [arrayU2 u_(2, :)'];
        arrayPsi1 = [arrayPsi1 psi(1, :)'];
        arrayPsi2 = [arrayPsi2 psi(2, :)'];
        
        % ���������� ������ �����������
        if timeCurrent < resultTime
            resultTime = timeCurrent;
            resultAngle = angle;
            resultX = x';
            resultPsi = psi;
            resultT = t';
            resultU = u_;
        end
        
        % ����������� ��� �������� ��������� ���������� ���������
        disp(num2str(angle) + " completed.");
    end
    
    resultAngleStep = angleStep;
    disp("Computation completed.");



% ���������� ��������� �����������
elseif optionStart == 3
    angleDivisor = input('Global angle divisor = ');
    nSteps = nSteps * angleDivisor;
    angleStep = angleStep / angleDivisor;
    resultTimePrevious = resultTime;
    resultAnglePrevious = resultAngle;
    resultXPrevious = resultX;

    iterCounter = 0;
    
    % �������� ����: ������� ������� ���������� ������ �� ��������� �����
    % ��������� ����� �� ������� - ������ ��������� "������" ������
    for angle = 0:angleStep:2*pi
        
        if mod(iterCounter, angleDivisor) == 0 % �� ������������� ��� ���������� ��������
            iterCounter = iterCounter + 1;
            continue
        end
        
        % ���������� �������������� ��������� ����� x0 (������ �������)
        psi0 = [cos(angle); sin(angle)]; % ������ psi(t0), ������-�������
        f = @(y) -(abs(psi0(1)) * sqrt(gamma - y .^ 4) + abs(psi0(2)) .* y); % ��������� ������������ �� ������� X0
        supportFunctionX0 = -(f(fminbnd(f, 0, gamma.^(1/4)))) + psi0(2) .* alpha; % �������� ������� ������� � X0 � ����������� psi(t0)
        h = @(x) abs(sum(psi0 .* x) - supportFunctionX0) + (x(1).^2 + (x(2) - alpha).^4 > gamma); % ������� ��� ������ ��������������� x0
        x0 = fminsearch(h, [0; alpha], optimset('MaxFunEvals', 1000)); % �������������� ��������� ����� � ��������� ������������ � ������ "����"
        
        % ���������� ��������������� ���������� u (������ �������)
        psiFunc = @(t) cell2mat(arrayfun(@(x) expm(-A' * (x - t0)) * psi0, t, 'UniformOutput', false)); % ������� 2xn �������� psi(t), ��� t - ������ 1xn 
        btPsi = @(t) B' * psiFunc(t); % ������� 2xn �������� B' * psi(t), ��� t - ������ 1xn
        supportFunctionP = @(t) sum(btPsi(t) .* p) + sqrt(sum(btPsi(t) .* (Q * btPsi(t)))); % ������ 1�n �������� ������� ������� � ����������� B'*psi(t) � P 
        u = @(t) p + (Q * btPsi(t)) ./ sqrt(sum(btPsi(t) .* (Q * btPsi(t))));
        %u = @(t) cell2mat(arrayfun(@(x) fminsearch(@(w) abs(sum(btPsi(x) .* w) - supportFunctionP(x)) + (sum((w-p) .* (Q * (w-p))) > 1), p), t, 'UniformOutput', false));
        
        % ������� ����������������� ���������
        tspan = [t0 resultTime]; % ��������� ��������� ����������
        options = odeset('Events', @(t, x) borderFunc(t, x, xLow, xUp, yLow, yUp)); % ��������� �������� ������ �� �������
        [t, x] = ode45(@(t, x) A * x + B * u(t), tspan, x0, options);
        
        % �������� ������� ����������� �������
        idxSolutionInX1 = find(x1Line(x(:, 1), x(:, 2)) > 0);
        inX1Size = size(idxSolutionInX1);
        if inX1Size(1) > 0
            idxFirst = idxSolutionInX1(1);
            timeCurrent = t(idxFirst);
        else
            tmp = size(t);
            idxFirst = min(tmp(1), constArrayLength);
            timeCurrent = resultTime + constTimePenalty;
        end
        
        % ������� "��������" ���������� � ��������� ���������� ����� ��������� � X1
        t = t(1:idxFirst);
        x = x(1:idxFirst, :);
        t(idxFirst:constArrayLength) = t(idxFirst);
        x(idxFirst:constArrayLength, 1) = x(idxFirst, 1);
        x(idxFirst:constArrayLength, 2) = x(idxFirst, 2);
        
        % ���������� �������� �����������
        u_ = u(t');
        psi = psiFunc(t');
        arrayT = [arrayT t];
        arrayX1 = [arrayX1 x(:, 1)];
        arrayX2 = [arrayX2 x(:, 2)];
        arrayU1 = [arrayU1 u_(1, :)'];
        arrayU2 = [arrayU2 u_(2, :)'];
        arrayPsi1 = [arrayPsi1 psi(1, :)'];
        arrayPsi2 = [arrayPsi2 psi(2, :)'];
        
        % ���������� ������ �����������
        if timeCurrent < resultTime
            resultTime = timeCurrent;
            resultAngle = angle;
            resultX = x';
            resultPsi = psi;
            resultT = t';
        end
        
        % ����������� ��� �������� ��������� ���������� ���������
        disp(num2str(angle) + " completed.");
        iterCounter = iterCounter + 1;
    end
    
    if resultTime < resultTimePrevious
        disp("Quality increased.");
    else
        disp("Quality remains the same.");
    end
    
    resultAngleStep = angleStep;
    
    disp("Computation completed.");

 
% ��������� �������� � ����������� ����������� ���������� (��������)   
elseif optionStart == 4
    resultTimePrevious = resultTime;
    resultAnglePrevious = resultAngle;
    resultXPrevious = resultX;
    localAngleStep = resultAngleStep / 10;
    margin = 1;
    arrayX1 = [];
    arrayX2 = [];
    arrayU1 = [];
    arrayU2 = [];
    arrayPsi1 = [];
    arrayPsi2 = [];
    arrayT = [];
    for angle = resultAngle - margin * resultAngleStep:localAngleStep:resultAngle + margin * resultAngleStep
        % ���������� �������������� ��������� ����� x0 (������ �������)
        psi0 = [cos(angle); sin(angle)]; % ������ psi(t0), ������-�������
        f = @(y) -(abs(psi0(1)) * sqrt(gamma - y .^ 4) + abs(psi0(2)) .* y); % ��������� ������������ �� ������� X0
        supportFunctionX0 = -(f(fminbnd(f, 0, gamma.^(1/4)))) + psi0(2) .* alpha; % �������� ������� ������� � X0 � ����������� psi(t0)
        h = @(x) abs(sum(psi0 .* x) - supportFunctionX0) + (x(1).^2 + (x(2) - alpha).^4 > gamma); % ������� ��� ������ ��������������� x0
        x0 = fminsearch(h, [0; alpha], optimset('MaxFunEvals', 1000)); % �������������� ��������� ����� � ��������� ������������ � ������ "����"
        
        % ���������� ��������������� ���������� u (������ �������)
        psiFunc = @(t) cell2mat(arrayfun(@(x) expm(-A' * (x - t0)) * psi0, t, 'UniformOutput', false)); % ������� 2xn �������� psi(t), ��� t - ������ 1xn 
        btPsi = @(t) B' * psiFunc(t); % ������� 2xn �������� B' * psi(t), ��� t - ������ 1xn
        supportFunctionP = @(t) sum(btPsi(t) .* p) + sqrt(sum(btPsi(t) .* (Q * btPsi(t)))); % ������ 1�n �������� ������� ������� � ����������� B'*psi(t) � P 
        u = @(t) p + (Q * btPsi(t)) ./ sqrt(sum(btPsi(t) .* (Q * btPsi(t))));
        %u = @(t) cell2mat(arrayfun(@(x) fminsearch(@(w) abs(sum(btPsi(x) .* w) - supportFunctionP(x)) + (sum((w-p) .* (Q * (w-p))) > 1), p), t, 'UniformOutput', false));
        
        % ������� ����������������� ���������
        tspan = [t0 resultTime]; % ��������� ��������� ����������
        options = odeset('Events', @(t, x) borderFunc(t, x, xLow, xUp, yLow, yUp)); % ��������� �������� ������ �� �������
        [t, x] = ode45(@(t, x) A * x + B * u(t), tspan, x0, options);
        
        % �������� ������� ����������� �������
        idxSolutionInX1 = find(x1Line(x(:, 1), x(:, 2)) > 0);
        inX1Size = size(idxSolutionInX1);
        if inX1Size(1) > 0
            idxFirst = idxSolutionInX1(1);
            timeCurrent = t(idxFirst);
        else
            tmp = size(t);
            idxFirst = min(tmp(1), constArrayLength);
            timeCurrent = resultTime + constTimePenalty;
        end
        
        % ������� "��������" ���������� � ��������� ���������� ����� ��������� � X1
        t = t(1:idxFirst);
        x = x(1:idxFirst, :);
        t(idxFirst:constArrayLength) = t(idxFirst);
        x(idxFirst:constArrayLength, 1) = x(idxFirst, 1);
        x(idxFirst:constArrayLength, 2) = x(idxFirst, 2);
        
        % ���������� �������� �����������
        u_ = u(t');
        psi = psiFunc(t');
        arrayT = [arrayT t];
        arrayX1 = [arrayX1 x(:, 1)];
        arrayX2 = [arrayX2 x(:, 2)];
        arrayU1 = [arrayU1 u_(1, :)'];
        arrayU2 = [arrayU2 u_(2, :)'];
        arrayPsi1 = [arrayPsi1 psi(1, :)'];
        arrayPsi2 = [arrayPsi2 psi(2, :)'];
        
        % ���������� ������ �����������
        if timeCurrent < resultTime
            resultTime = timeCurrent;
            resultAngle = angle;
            resultX = x';
            resultPsi = psi;
            resultT = t';
            resultU = u_;
        end
        
        % ����������� ��� �������� ��������� ���������� ���������
        disp(num2str(angle) + " completed.");
    end
    
    if resultTime < resultTimePrevious
        disp("Quality increased.");
    else
        disp("Quality remains the same.");
    end
    
    disp("Computation completed.");


% ���������� ��������
elseif optionStart == 5
    disp("What kind of plot do you want?");
    disp("1) x1, x2");
    disp("2) t, x1");
    disp("3) t, x2");
    disp("4) t, u1");
    disp("5) t, u2");
    disp("6) t, psi1");
    disp("7) t, psi2");
    disp("8) u1, u2");
    disp("9) Transversality Condition in X1");
    optionPlot = input("Plot � ");
    disp("Do you want to see sub-optiomal solutions? 1 - yes, 2 - no");
    optionSubOptimal = input("Sub-optimal solutions? ");
    figure;
    hold on;
    if optionPlot == 1
        x0Line = @(x, y) (x.^2 + (y - alpha).^4 <= gamma);
        x_ = linspace(-sqrt(gamma) - 1, sqrt(gamma) + 1, 1001);
        y_ = linspace(alpha - sqrt(gamma) - 1, alpha + sqrt(gamma) + 1, 1001);
        [X_, Y_] = meshgrid(x_, y_);
        Z_ = x0Line(X_, Y_);
        contourf(X_, Y_, Z_, [1 1]);
        x1Line = @(x, y) (x .* y >= 1) .* (x + y <= a) .* (x >= 0);
        x_ = linspace(0, (a + sqrt(d)) / 2 + 2, 1001);
        y_ = linspace(0, (a + sqrt(d)) / 2 + 2, 1001);
        [X_, Y_] = meshgrid(x_, y_);
        Z_ = x1Line(X_, Y_);
        contourf(X_, Y_, Z_, [1 1]);
        if optionSubOptimal == 1
            pl = plot(arrayX1, arrayX2, 'b');
            if resultTime < t0 + constTimeLimit
                plot(resultX(1, :), resultX(2, :), 'r*');
            end
        else
            pl = plot(resultX(1, :), resultX(2, :));
        end
        xlabel('x1');
        ylabel('x2');
        %legend;
    elseif optionPlot == 2
        if optionSubOptimal == 1
            pl = plot(arrayT, arrayX1);
            if resultTime < t0 + constTimeLimit
                plot(resultT, resultX(1, :), 'r*');
            end
        else
            pl = plot(resultT, resultX(1, :)); 
        end
        xlabel('t');
        ylabel('x1');
        %legend;
    elseif optionPlot == 3
        if optionSubOptimal == 1
            pl = plot(arrayT, arrayX2);
            if resultTime < t0 + constTimeLimit
                plot(resultT, resultX(2, :), 'r*');
            end
        else
            pl = plot(resultT, resultX(2, :)); 
        end
        xlabel('t');
        ylabel('x2');
        %legend;
    elseif optionPlot == 4
        if optionSubOptimal == 1
            pl = plot(arrayT, arrayU1);
            if resultTime < t0 + constTimeLimit
                plot(resultT, resultU(1, :), 'r*');
            end
        else
            pl = plot(resultT, resultU(1, :)); 
        end
        xlabel('t');
        ylabel('u1');
        %legend;
    elseif optionPlot == 5
        if optionSubOptimal == 1
            plot(arrayT, arrayU2);
            if resultTime < t0 + constTimeLimit
                pl = plot(resultT, resultU(2, :), 'r*');
            end
        else
            pl = plot(resultT, resultU(2, :)); 
        end
        xlabel('t');
        ylabel('u2');
        %legend;
    elseif optionPlot == 6
        if optionSubOptimal == 1
            pl = plot(arrayT, arrayPsi1);
            if resultTime < t0 + constTimeLimit
                plot(resultT, resultPsi(1, :), 'r*');
            end
        else
            pl = plot(resultT, resultPsi(1, :)); 
        end
        xlabel('t');
        ylabel('psi1');
        %legend;
    elseif optionPlot == 7
        if optionSubOptimal == 1
            pl = plot(arrayT, arrayPsi2);
            if resultTime < t0 + constTimeLimit
                plot(resultT, resultPsi(2, :), 'r*');
            end
        else
            pl = plot(resultT, resultPsi(2, :)); 
        end
        xlabel('t');
        ylabel('psi2');
    elseif optionPlot == 8
        if optionSubOptimal == 1
            pl = plot(arrayU1, arrayU2);
            if resultTime < t0 + constTimeLimit
                plot(resultU(1, :), resultU(2, :), 'r*');
            end
        else
            pl = plot(resultU(1, :), resultU(2, :)); 
        end
        xlabel('u1');
        ylabel('u2');
        %legend;
    elseif optionPlot == 9
        if resultTime == t0 + constTimeLimit
            disp("No optimal control found");
            return
        end
        x0Line = @(x, y) (x.^2 + (y - alpha).^4 <= gamma);
        x_ = linspace(-sqrt(gamma) - 1, sqrt(gamma) + 1, 1001);
        y_ = linspace(alpha - sqrt(gamma) - 1, alpha + sqrt(gamma) + 1, 1001);
        [X_, Y_] = meshgrid(x_, y_);
        Z_ = x0Line(X_, Y_);
        contourf(X_, Y_, Z_, [1 1]);
        x1Line = @(x, y) (x .* y >= 1) .* (x + y <= a) .* (x >= 0);
        x_ = linspace(0, (a + sqrt(d)) / 2 + 2, 1001);
        y_ = linspace(0, (a + sqrt(d)) / 2 + 2, 1001);
        [X_, Y_] = meshgrid(x_, y_);
        Z_ = x1Line(X_, Y_);
        contourf(X_, Y_, Z_, [1 1]);
        
        scp = -(sum(resultPsi(:, end) .* resultX(:, end)));
        tmp = @(x) (-1.7) * sum(-resultPsi(:, end) .* x) - x1Line(x(1), x(2));
        x_sup = fminsearch(tmp, [a/2; ((a.^2)/(4*a))], optimset('MaxFunEvals', 1000));
        spf = sum(-resultPsi(:, end) .* x_sup);
        normal = @(x, x0) x0.^2 .*(x - x0) + (1 / x0);
        x_norm = [0, resultX(1, end)];
        y_norm = [normal(0, resultX(1, end)), resultX(2, end)];
        %x_plot0 = [resultX(1, 1), resultX(1, 1) + resultPsi(1, 1)];
        %y_plot0 = [resultX(2, 1), resultX(2, 1) + resultPsi(2, 1)];
        x_plot1 = [resultX(1, end), resultX(1, end) - resultPsi(1, end)];
        y_plot1 = [resultX(2, end), resultX(2, end) - resultPsi(2, end)];
        pl = plot(resultX(1, :), resultX(2, :), x_plot1, y_plot1, resultX(1, end), resultX(2, end), 'y*', x_norm, y_norm);
        xlabel('x1');
        ylabel('x2');
        legend('������� �0', '������� �1', '����������� ����������', '������ -psi(t1)', '�������� �����', '������ �������');
        disp(abs(spf - scp) / sqrt(sum(resultPsi(:, end) .* resultPsi(:, end))));
    else
        disp("ERROR. Wrong option");
    end
    for i = 1:size(pl)
        pl(i).LineWidth = 2;
    end
    hold off;
    disp("Plot is ready.");
    

% ����� �� ���������    
elseif optionStart == 6
    disp("Goodbye!");

% �������� ����
else
    disp("ERROR. Wrong option.");
end

function [value,isterminal,direction] = borderFunc(t,x, xLow, xUp, yLow, yUp)
isterminal = 1;
direction = 0;
if x(1) > xUp || x(1) < xLow
    value = 0;
elseif x(2) > yUp || x(2) < yLow
    value = 0;
else
    value = 1;
end
end