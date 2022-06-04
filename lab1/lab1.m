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
% Опции алгоритма:
% 1) Настройка параметров, 
% 2) Проведение вычислений 
% 3) Улучшение ответа (глобальное)
% 4) Улучшение ответа (локальное)
% 5) Построение графиков
% 6) Выход

% Объявление констант
constTimeLimit = 5;
constTimePenalty = 0.1;
constEps = 0.0001;
constArrayLength = 500;
constNSteps = 30;

%Начало обработки команд
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

% Настройка алгоритма
if optionStart == 1
    A = input("A = "); % Матрица 2х2
    B = input("B = "); % Матрица 2х2
    t0 = input("t0 = "); % Число
    Q = input("Q = "); % Матрица 2х2
    p = input("p = "); % Вектор-столбец 1х2
    a = input("a = "); % Число 
    alpha = input("alpha = "); % Число > 0
    gamma = input("gamma = "); % Число

    
% Вычисление результатов
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
    
    % Проверка корректности условия
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
    
    % Границы графиков (нужны для более удобных графиков в пункте 5)
    xBorders = [-2*sqrt(gamma) - 2, max(2*sqrt(gamma), a + 2 * sqrt(d) + 1)];
    yBorders = [ min(alpha - 2*gamma, 0), max(alpha + 2*gamma, 1 /(a - sqrt(d)) + 1)];
    xLow = xBorders(1) - 10;
    xUp = xBorders(2) + 10;
    yLow = yBorders(1) - 10;
    yUp = yBorders(2) + 10;
    
    % Графическое представление множеств X0 и X1
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
    
    % Если множества пересекаются, искомое время равно 0
    if max(max(x0Line(X_, Y_) .* x1Line(X_, Y_))) > 0 % Проверка на пересечение границ
        disp("X0 intersects X1. The time required is equal 0.");
        isIntersected = 1;
        hold off;
        return
    else
        isIntersected = 0;
    end
    
    % Регуляризация (при отсутствии вполне управляемости)
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
    
    % Основной цикл: перебор решений сопряжённой задачи на единичной сфере
    for angle = 0:angleStep:2*pi
        
        % Нахождение субоптимальной начальной точки x0 (вектор столбец)
        psi0 = [cos(angle); sin(angle)]; % Вектор psi(t0), вектор-столбец
        f = @(y) -(abs(psi0(1)) * sqrt(gamma - y .^ 4) + abs(psi0(2)) .* y); % Скалярное произведение на границе X0
        supportFunctionX0 = -(f(fminbnd(f, 0, gamma.^(1/4)))) + psi0(2) .* alpha; % Значение опорной функции к X0 в направлении psi(t0)
        h = @(x) abs(sum(psi0 .* x) - supportFunctionX0) + (x(1).^2 + (x(2) - alpha).^4 > gamma); % Функция для поиска субоптимального x0
        x0 = fminsearch(h, [0; alpha], optimset('MaxFunEvals', 1000)); % Субоптимальная начальная точка с начальным приближением в центре "шара"
        
        % Нахождение субоптимального управления u (вектор столбец)
        psiFunc = @(t) cell2mat(arrayfun(@(x) expm(-A' * (x - t0)) * psi0, t, 'UniformOutput', false)); % Матрица 2xn значений psi(t), где t - строка 1xn 
        btPsi = @(t) B' * psiFunc(t); % Матрица 2xn значение B' * psi(t), где t - строка 1xn
        supportFunctionP = @(t) sum(btPsi(t) .* p) + sqrt(sum(btPsi(t) .* (Q * btPsi(t)))); % Строка 1хn значений опорной функции в направлении B'*psi(t) к P 
        u = @(t) p + (Q * btPsi(t)) ./ sqrt(sum(btPsi(t) .* (Q * btPsi(t))));
        %u = @(t) cell2mat(arrayfun(@(x) fminsearch(@(w) abs(sum(btPsi(x) .* w) - supportFunctionP(x)) + (sum((w-p) .* (Q * (w-p))) > 1), p), t, 'UniformOutput', false));
        
        % Решение дифференциального уравнения
        tspan = [t0 resultTime]; % Доступный временной промежуток
        options = odeset('Events', @(t, x) borderFunc(t, x, xLow, xUp, yLow, yUp)); % Настройка проверки выхода за границы
        [t, x] = ode45(@(t, x) A * x + B * u(t), tspan, x0, options);
        
        % Проверка времени полученного решения
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
        
        % Обрезка "ненужной" информации о поведении траектории после попадания в X1 и стандартизация результатов
        t = t(1:idxFirst);
        x = x(1:idxFirst, :);
        t(idxFirst:constArrayLength) = t(idxFirst);
        x(idxFirst:constArrayLength, 1) = x(idxFirst, 1);
        x(idxFirst:constArrayLength, 2) = x(idxFirst, 2);
        
        % Пополнение массивов результатов
        u_ = u(t');
        psi = psiFunc(t');
        arrayT = [arrayT t];
        arrayX1 = [arrayX1 x(:, 1)];
        arrayX2 = [arrayX2 x(:, 2)];
        arrayU1 = [arrayU1 u_(1, :)'];
        arrayU2 = [arrayU2 u_(2, :)'];
        arrayPsi1 = [arrayPsi1 psi(1, :)'];
        arrayPsi2 = [arrayPsi2 psi(2, :)'];
        
        % Обновление лучших результатов
        if timeCurrent < resultTime
            resultTime = timeCurrent;
            resultAngle = angle;
            resultX = x';
            resultPsi = psi;
            resultT = t';
            resultU = u_;
        end
        
        % Уведомление для контроля прогресса выполнения программы
        disp(num2str(angle) + " completed.");
    end
    
    resultAngleStep = angleStep;
    disp("Computation completed.");



% Глобальное улучшение результатов
elseif optionStart == 3
    angleDivisor = input('Global angle divisor = ');
    nSteps = nSteps * angleDivisor;
    angleStep = angleStep / angleDivisor;
    resultTimePrevious = resultTime;
    resultAnglePrevious = resultAngle;
    resultXPrevious = resultX;

    iterCounter = 0;
    
    % Основной цикл: перебор решений сопряжённой задачи на единичной сфере
    % Стартовый лимит по времени - лучший результат "старой" версии
    for angle = 0:angleStep:2*pi
        
        if mod(iterCounter, angleDivisor) == 0 % Не пересчитываем уже полученные значения
            iterCounter = iterCounter + 1;
            continue
        end
        
        % Нахождение субоптимальной начальной точки x0 (вектор столбец)
        psi0 = [cos(angle); sin(angle)]; % Вектор psi(t0), вектор-столбец
        f = @(y) -(abs(psi0(1)) * sqrt(gamma - y .^ 4) + abs(psi0(2)) .* y); % Скалярное произведение на границе X0
        supportFunctionX0 = -(f(fminbnd(f, 0, gamma.^(1/4)))) + psi0(2) .* alpha; % Значение опорной функции к X0 в направлении psi(t0)
        h = @(x) abs(sum(psi0 .* x) - supportFunctionX0) + (x(1).^2 + (x(2) - alpha).^4 > gamma); % Функция для поиска субоптимального x0
        x0 = fminsearch(h, [0; alpha], optimset('MaxFunEvals', 1000)); % Субоптимальная начальная точка с начальным приближением в центре "шара"
        
        % Нахождение субоптимального управления u (вектор столбец)
        psiFunc = @(t) cell2mat(arrayfun(@(x) expm(-A' * (x - t0)) * psi0, t, 'UniformOutput', false)); % Матрица 2xn значений psi(t), где t - строка 1xn 
        btPsi = @(t) B' * psiFunc(t); % Матрица 2xn значение B' * psi(t), где t - строка 1xn
        supportFunctionP = @(t) sum(btPsi(t) .* p) + sqrt(sum(btPsi(t) .* (Q * btPsi(t)))); % Строка 1хn значений опорной функции в направлении B'*psi(t) к P 
        u = @(t) p + (Q * btPsi(t)) ./ sqrt(sum(btPsi(t) .* (Q * btPsi(t))));
        %u = @(t) cell2mat(arrayfun(@(x) fminsearch(@(w) abs(sum(btPsi(x) .* w) - supportFunctionP(x)) + (sum((w-p) .* (Q * (w-p))) > 1), p), t, 'UniformOutput', false));
        
        % Решение дифференциального уравнения
        tspan = [t0 resultTime]; % Доступный временной промежуток
        options = odeset('Events', @(t, x) borderFunc(t, x, xLow, xUp, yLow, yUp)); % Настройка проверки выхода за границы
        [t, x] = ode45(@(t, x) A * x + B * u(t), tspan, x0, options);
        
        % Проверка времени полученного решения
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
        
        % Обрезка "ненужной" информации о поведении траектории после попадания в X1
        t = t(1:idxFirst);
        x = x(1:idxFirst, :);
        t(idxFirst:constArrayLength) = t(idxFirst);
        x(idxFirst:constArrayLength, 1) = x(idxFirst, 1);
        x(idxFirst:constArrayLength, 2) = x(idxFirst, 2);
        
        % Пополнение массивов результатов
        u_ = u(t');
        psi = psiFunc(t');
        arrayT = [arrayT t];
        arrayX1 = [arrayX1 x(:, 1)];
        arrayX2 = [arrayX2 x(:, 2)];
        arrayU1 = [arrayU1 u_(1, :)'];
        arrayU2 = [arrayU2 u_(2, :)'];
        arrayPsi1 = [arrayPsi1 psi(1, :)'];
        arrayPsi2 = [arrayPsi2 psi(2, :)'];
        
        % Обновление лучших результатов
        if timeCurrent < resultTime
            resultTime = timeCurrent;
            resultAngle = angle;
            resultX = x';
            resultPsi = psi;
            resultT = t';
        end
        
        % Уведомление для контроля прогресса выполнения программы
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

 
% Улучшение качества в окрестности полученного результата (локально)   
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
        % Нахождение субоптимальной начальной точки x0 (вектор столбец)
        psi0 = [cos(angle); sin(angle)]; % Вектор psi(t0), вектор-столбец
        f = @(y) -(abs(psi0(1)) * sqrt(gamma - y .^ 4) + abs(psi0(2)) .* y); % Скалярное произведение на границе X0
        supportFunctionX0 = -(f(fminbnd(f, 0, gamma.^(1/4)))) + psi0(2) .* alpha; % Значение опорной функции к X0 в направлении psi(t0)
        h = @(x) abs(sum(psi0 .* x) - supportFunctionX0) + (x(1).^2 + (x(2) - alpha).^4 > gamma); % Функция для поиска субоптимального x0
        x0 = fminsearch(h, [0; alpha], optimset('MaxFunEvals', 1000)); % Субоптимальная начальная точка с начальным приближением в центре "шара"
        
        % Нахождение субоптимального управления u (вектор столбец)
        psiFunc = @(t) cell2mat(arrayfun(@(x) expm(-A' * (x - t0)) * psi0, t, 'UniformOutput', false)); % Матрица 2xn значений psi(t), где t - строка 1xn 
        btPsi = @(t) B' * psiFunc(t); % Матрица 2xn значение B' * psi(t), где t - строка 1xn
        supportFunctionP = @(t) sum(btPsi(t) .* p) + sqrt(sum(btPsi(t) .* (Q * btPsi(t)))); % Строка 1хn значений опорной функции в направлении B'*psi(t) к P 
        u = @(t) p + (Q * btPsi(t)) ./ sqrt(sum(btPsi(t) .* (Q * btPsi(t))));
        %u = @(t) cell2mat(arrayfun(@(x) fminsearch(@(w) abs(sum(btPsi(x) .* w) - supportFunctionP(x)) + (sum((w-p) .* (Q * (w-p))) > 1), p), t, 'UniformOutput', false));
        
        % Решение дифференциального уравнения
        tspan = [t0 resultTime]; % Доступный временной промежуток
        options = odeset('Events', @(t, x) borderFunc(t, x, xLow, xUp, yLow, yUp)); % Настройка проверки выхода за границы
        [t, x] = ode45(@(t, x) A * x + B * u(t), tspan, x0, options);
        
        % Проверка времени полученного решения
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
        
        % Обрезка "ненужной" информации о поведении траектории после попадания в X1
        t = t(1:idxFirst);
        x = x(1:idxFirst, :);
        t(idxFirst:constArrayLength) = t(idxFirst);
        x(idxFirst:constArrayLength, 1) = x(idxFirst, 1);
        x(idxFirst:constArrayLength, 2) = x(idxFirst, 2);
        
        % Пополнение массивов результатов
        u_ = u(t');
        psi = psiFunc(t');
        arrayT = [arrayT t];
        arrayX1 = [arrayX1 x(:, 1)];
        arrayX2 = [arrayX2 x(:, 2)];
        arrayU1 = [arrayU1 u_(1, :)'];
        arrayU2 = [arrayU2 u_(2, :)'];
        arrayPsi1 = [arrayPsi1 psi(1, :)'];
        arrayPsi2 = [arrayPsi2 psi(2, :)'];
        
        % Обновление лучших результатов
        if timeCurrent < resultTime
            resultTime = timeCurrent;
            resultAngle = angle;
            resultX = x';
            resultPsi = psi;
            resultT = t';
            resultU = u_;
        end
        
        % Уведомление для контроля прогресса выполнения программы
        disp(num2str(angle) + " completed.");
    end
    
    if resultTime < resultTimePrevious
        disp("Quality increased.");
    else
        disp("Quality remains the same.");
    end
    
    disp("Computation completed.");


% Построение графиков
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
    optionPlot = input("Plot № ");
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
        legend('Граница Х0', 'Граница Х1', 'Оптимальная траектория', 'Вектор -psi(t1)', 'Конечная точка', 'Вектор нормали');
        disp(abs(spf - scp) / sqrt(sum(resultPsi(:, end) .* resultPsi(:, end))));
    else
        disp("ERROR. Wrong option");
    end
    for i = 1:size(pl)
        pl(i).LineWidth = 2;
    end
    hold off;
    disp("Plot is ready.");
    

% Выход из программы    
elseif optionStart == 6
    disp("Goodbye!");

% Неверный ввод
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