%% TASK 1
a = input('');
b = input('');
n = input('');
if ~isnumeric(a) || ~isnumeric(b) || ~isnumeric(n) || n <= 1 || mod(n, 1) >= 0.0000001 % Проверка корректности условия
    clear;
    error('Wrong arguments')
end
if a > b % Упорядочивание концов
    t = a;
    a = b;
    b = t;
end
dots = linspace(a, b, n); % Сетка точек
values = func_1_1(dots); % Значения в сетке
[y_max, max_index] = max(values); % y_max - максимум, max_index - индекс максимума в массиве
[y_min, min_index] = min(values); % Аналогично для минимума
x_min = dots(min_index); % Находим абсциссы точек мин.
x_max = dots(max_index); % Находим абсциссы точек макс.
disp("Глобальный минимум - " + num2str(y_min)); % Глобальный мин.(значение)
disp("Глобальный максимум - " + num2str(y_max)); % Глобальный макс.(значение)
plot(dots, values);
xlabel("x");
ylabel("x^2 cos(|x|)");
legend("y = x^2 cos(|x|)");
repmat(y_max, size(x_max)); % Выравниваем размеры y_max и x_max
repmat(y_min, size(x_min)); % Аналогично для мин.
text(x_min, y_min, '\leftarrow y_{min}') % Помечаем на графике точки минимума стрелочкой
text(x_max, y_max, '\leftarrow y_{max}') % Аналогично для максимума
clear; % Очистка переменных блока
%% TASK 2
n = input('');
disp("Проверка на простоту числа n выдала результат: " + num2str(isprime(n)));
disp("Поиск чисел, удовлетворяющих условию задачи...")
if n < 7
    disp("Таких чисел нет!");
else
    %x = 1:n / 7:2; % Массив нечётных целых чисел от 1 до n/7
    %disp(7 * x); % Искомые числа
    disp(7:14:n);
end
res1 = (ones(n) .* (2:(n + 1)))'; % Матрица n*n из единиц, потом умножаем по строкам на 2,3,4...
disp(res1);
res2 = reshape(1:((n + 1) ^ 2), [n + 1, n + 1])'; % Нужны элементы от 1 до (n+1)^2, далее меняем размер и транспонируем
disp(res2);
res3 = reshape(res2(n^2 : (n + 1)^2), [n + 1, 2]); % Явно вычисляем, какие значения будут в этих двух столбцах и выводим их в нужном формате
disp(res3);
clear;
%% TASK 3
matrix = round(rand(5, 8) * 21 - 10.5); % Задаём матрицу с равн. распред. значениями на [0, 1],
% Домножаем на 21 (т.е. уже U[0, 21]) и вычитаем 10.5 (т.е. U[-10.5, 10.5])
% Далее в каждое целое число переходит промежуток (n-1/2, n+1/2), границы
% неважны, тк вероятность попадания в точку равна 0.
res1 = max(diag(matrix)); % Максимальный элемент на диагонали
sums = sum(matrix, 1); % Построчные суммы
prods = prod(matrix, 1); % Построчные произведения
res2 = max(sums ./ prods); % Максимальное отношение построчных сумм к построчным произведениям
res3 = min(sums ./ prods); % Минимальное отношение построчных сумм к построчным произведениям
new_matrix = sortrows(matrix, 'descend'); % Сортирует в обратном лексикографическом порядке
disp("Наибольший диагональный элемент - " + num2str(res1));
disp("Наибольшее отношение суммы к произведениям: " + num2str(res2));
disp("Наименьшее отношение суммы к произведениям: " + num2str(res3));
disp(new_matrix);
clear;
%% TASK 4
n = input('');
m = input('');
matrix = rand(n, m);
[r, g, b] = func_1_4(matrix, n, m);
disp("MATRIX = ");
disp(matrix);
disp("R = ");
disp(r);
disp("G = ");
disp(g);
disp("B = ");
disp(b);
clear;
%% TASK 5
n = input('');
m = input('');
x = rand(1, n); %Случайный тестовый вектор x
y = rand(1, m); %Случайный тестовый вектор y
res = func_1_5(x, y, n, m);
disp(x);
disp(y);
disp(res);
clear;
%% TASK 6
n = input('');
x = rand(3, n); % Матрица случайных векторов
x_new = reshape(repmat(x, n, 1), [3, n * n]); % Делаем аналог задачи 5: декартово произведение  мн-ва векторов на себя
y_new = repmat(x, 1, n); 
crosses = cross(x_new, y_new); % Из декартова произведения найдём векторные произведения
abses = reshape(sqrt(sum(crosses .* crosses, 1)), [n, n]); % А из векторных произведений - их модули. Далее нужно лишь сделать матрицу квадратной.
disp(abses);
clear;
%% TASK 7
% Проверим работу функции на частных примерах
a1 = [1, 2, 3, 4];
b1 = [-1, -2, -3, -4];
a2 = [1, 2, 3];
b2 = [3, 4, 5];
disp(func_1_7(a1, b1));
disp(func_1_7(a2, b2));
clear;
%% TASK 8
n = input('');
k = input('');
x = rand(k, n); % Матрица случайных векторов
x_new = reshape(repmat(x, n, 1), [k, n * n]); % Делаем аналог задачи 5: декартово произведение  мн-ва векторов на себя
y_new = repmat(x, 1, n);
diffs = reshape(sqrt(sum((x_new - y_new) .* (x_new - y_new), 1)), [n, n]); % Аналог задачи 6 для расстояний
disp(diffs);
clear;
%% TASK 9
% Для уменьшения объема кода будем рассматривать разницу времени работы
% функций только для квадратных матриц
n_times = 100;
max_dim = 300;
slow_average = zeros(1, max_dim); % Среднее время работы своей функции на матрицах различного размера
fast_average = zeros(1, max_dim); % Среднее время работы встроенного метода на матрицах различного размера
for n = 1:max_dim
    slow_res = zeros(1, n_times); % Результаты работы своей функции на матрицах данного размера
    fast_res = zeros(1, n_times); % Результаты работы встроенного метода на матрицах данного размера
    for i = 1:n_times
        a = rand(n);
        b = rand(n);
        tic();
        func_1_9(a, b, n, n);
        slow_res(i) = toc();
        tic();
        a + b;
        fast_res(i) = toc();
    end
    slow_average(n) = mean(slow_res);
    fast_average(n) = mean(fast_res);
end
plot(1:max_dim, slow_average, 1:max_dim, fast_average);
xlabel("Размерность матрицы");
ylabel("Среднее время вычислений (с)");
legend("Поэлементное сложение", "Встроенный метод");
clear;
%% TASK 10
% Проверим работу функции на частных случаях
a1 = [1, 2, 3, 2, 1];
a2 = [1, 2, 3, 4, 1];
a3 = (1);
a4 = [2, 2];
a5 = [6, 9];
disp(func_1_10(a1)); % 1
disp(func_1_10(a2)); % 0
disp(func_1_10(a3)); % 1
disp(func_1_10(a4)); % 1
disp(func_1_10(a5)); % 0
clear;
%% TASK 11
a = input('');
if a <= 0
    clear;
    error("Неверная граница отрезка");
end
b = input('');
n = input('');
x = rand(1, n) * a; % n случайных величин равномерно распределенных на [0, a]
num = sum((x >= b)) / n; % Искомая доля чисел, больших b
if num <= a / (2 * b)
    disp("Согласуется с гипотезой!");
else
    disp("Что-то пошло не так...");
end
%disp(x);
disp(num);
disp(a / (2 * b));
clear;
%% TASK 12
% Будем строить такую первообразную exp(-x^2), которая в левом конце равна 0: F(a)= 0
% Задание параметров
a = -2; % Зададим левый конец отрезка
b = 2; % Зададим правый конец отрезка
n = 30; % Зададим число отрезков разбиения
m = n; % Количество отрезков разбиения подотрезка [x_k, x_k+1]
h = (b - a) / n; % Размер большого шага
gauss = @(x) exp(-x .* x); % Анонимная функция exp(-x^2)
tiledlayout(2, 2);
%Формула прямоугольников
[x_rect, y_rect] = rectangles(a, b, n, m, gauss);
% Формула Симпсона
[x_sim, y_sim] = simpson(a, b, n, m, gauss);
% Формула трапеций
x_trap = linspace(a, b, n + 1);
y_trap = zeros(1, n);
for k = 1:n
    x_k = linspace(a, a + k * h, m * k + 1);
    y_k = gauss(x_k);
    y_trap(k) = trapz(x_k, y_k);
end;
y_trap = horzcat(0, y_trap);
% Общий график
ax1 = nexttile;
plot(ax1, x_rect, y_rect, x_sim, y_sim, x_trap, y_trap);
xlabel("x");
ylabel("Первообразная exp(-x^2)");
legend("Метод прямоугольников", "Метод Симпсона", "Метод трапеций");
% Время работы методов
ax2 = nexttile;
n_times = 100;
n_dots = 100;
rect_times = zeros(1, n_dots);
sim_times = zeros(1, n_dots);
trap_times = zeros(1, n_dots);
rect_cur = zeros(1, n_times);
sim_cur = zeros(1, n_times);
trap_cur = zeros(1, n_times);
for k = 1:n_dots
    for i = 1:n_times
        tic();
        rectangles(a, b, n, k, gauss);
        rect_cur(i) = toc();
        tic()
        simpson(a, b, n, k, gauss);
        sim_cur(i) = toc();
        tic()
        x_k = linspace(a, a + h, k + 1);
        y_k = gauss(x_k);
        trapz(x_k, y_k);
        trapz_cur(i) = (n + 1) * toc();
    end
    rect_times(k) = mean(rect_cur);
    sim_times(k) = mean(sim_cur);
    trap_times(k) = mean(trap_cur);
end
plot(ax2, 1:n_dots, rect_times, 1:n_dots, sim_times, 1:n_dots, trap_times);
xlabel("Количество отрезков разбиения");
ylabel("Среднее время на вычисление (с)");
legend("Метод прямоугольников", "Метод Симпсона", "Метод трапеций");
% Сравнение скорости внутренней сходимости
ax3 = nexttile;
n_dots = 100;
n_start = 3;
rect_diff = zeros(1, n_dots);
sim_diff = zeros(1, n_dots);
trap_diff = zeros(1, n_dots);
for k = n_start:n_dots
    x_tmp1 = linspace(a, b, k + 1);
    y_tmp1 = gauss(x_tmp1);
    x_tmp2 = linspace(a, b, 2 * k + 1);
    y_tmp2 = gauss(x_tmp2);
    [tmp, rect_k] = rectangles(a, b, k, 1, gauss);
    [tmp, rect_2k] = rectangles(a, b, 2 * k, 1, gauss);
    [tmp, sim_k] = simpson(a, b, k, 1, gauss);
    [tmp, sim_2k] = simpson(a, b, 2 * k, 1, gauss);
    trap_k = trapz(x_tmp1, y_tmp1);
    trap_2k = trapz(x_tmp2, y_tmp2);
    rect_diff(k) = abs(rect_k(end) - rect_2k(end));
    sim_diff(k) = abs(sim_k(end) - sim_2k(end));
    trap_diff(k) = abs(trap_k - trap_2k);
end;
rect_diff = rect_diff(n_start:n_dots);
sim_diff = sim_diff(n_start:n_dots);
trap_diff = trap_diff(n_start:n_dots);
plot(ax3, n_start:n_dots, rect_diff, n_start:n_dots, sim_diff, n_start:n_dots, trap_diff);
xlabel("Количество отрезков разбиения");
ylabel("Модуль разности на n и 2n отрезках");
legend("Метод прямоугольников", "Метод Симпсона", "Метод трапеций");
clear;
%% TASK 13
n = 100; % Число различных вариантов шага
x = ones(1, n); % x = 1 для каждого из n вариантов шага
h = logspace(-10, -1, n); % Варианты шага
abs_right = abs(func_1_13_2(x) - (func_1_13_1(x + h) - func_1_13_1(x)) ./ h); % Массив модулей разности с правой разностной производной
abs_center = abs(func_1_13_2(x) - (func_1_13_1(x + h) - func_1_13_1(x - h)) ./(2 * h)); % Массив модулей разности с центральной разностной производной
loglog(h, abs_center, h, abs_right);
xlabel("Размер шага");
ylabel("Модуль отклонения от производной");
legend("Центр. разн. пр-ая", "Прав. разн. пр-ая");
clear;
%% TEST
square = @(x) x .* x;

clear;
