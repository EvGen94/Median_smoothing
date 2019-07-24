%%  Median smoothing algoritms 
% The first part using built-in functions, the second part without built-in functions
% Algorithmic complexity for Median smoothing in the second part determined by 
% the complexity of sorting algorithmic and finding the median element for
% each group and is defined as O(TN)
%--------------------------------------------------------------------------
power = 1; % data multiplier
max = 200*power; % upper boundary of input data
n = 20; % lower window boundary of input data
m = randi([n max-n]);  % upper window boundary of input data
time = 1:2:max; % time scale 

x = 2*(time)+rand(1,length(time))*max; % time series, equation of a straight line with a uniformly distributed random component
y = medfilt1(x,round((m-n)/2)); % median filtering, the number of the median smoothing members in the median smoothing groop is calculated by the expression (m-n)/2

Medianfigure = figure('Name','the first part'); % creating figure
plot(time,x,time,y); % ploting 'Original' and 'Median'smoothing graphics
axis([n m 0 inf]); % window boundaries
title(['Window width is ',num2str(m-n),'. Median smoothing in groups of  = ',num2str(round((m-n)/2))]);
legend('Original','Median');
legend('boxoff');

%% Part two,without built-in functions
clear all
power=1; % data multiplier
max = 100*power; % upper boundary of input data
n = 20;  % lower window boundary of input data
m = randi([n max-n]); % upper window boundary of input data
time = 1:2:max; % time scale
group = round((m-n)/2); % Median smoothing in groups of (m-n)/2 
x = 2*(time)+rand(1,length(time))*max; % time series, equation of a straight line with a uniformly distributed random component
% loop for median smoothing
for i = 1:length(x)- group 
% formation of median groups   
x1 = x(i:group+i-1);
y = bucketsort(x1); % Bucket sort of median group, linear complexity O(n) 
    if rem(length(y),2) == 1  % Check even or odd
        median(i) = y(round(length(y)/2)) ;
    else
        index1 = round(length(y)/2);
        index2 = index1 + 1 ;
        median(i) = (y(index1)+y(index2))/2 ;
    end

end
centring = (length(time)-length(median))/2; % time scale centring
time2 = time(centring+1:length(time)-centring);% median time scale
Medianfigure = figure('Name','the second part'); % creating figure
plot(time,x,time2,median); % ploting 'Original' and 'Median'smoothing graphics
axis([group time(end)- group 0 inf]); % window 
title(['Window width is ',num2str(time(end)- 2*group),'. Median smoothing in groups of  = ',num2str(group)]);
legend('Original','Median');
legend('boxoff');

function sx = bucketsort(x)
% Default load factor
alpha = 0.75; % alpha = n / m
% Find min and max elements of x
n = length(x);
[minx maxx] = minmax(x,n);
% Insert elements into m equal width buckets, each containing a doubly
% linked list
m = round(n / alpha);
dw = (maxx - minx) / m;
head = nan(1,m); % pointers to heads of bucket lists
prev = nan(1,n); % previous element pointers
next = nan(1,n); % next element pointers
last = nan(1,m); % temporary storage
for i = 1:n
    j = min(floor((x(i) - minx) / dw) + 1,m); % hack to make max(x) fall in last bucket
    if isnan(head(j))
        head(j) = i;
    else
        prev(i) = last(j);
        next(last(j)) = i;
    end
    last(j) = i;
end
% Bucket sort
sx = zeros(size(x)); % sorted array
kk = 0;
for j = 1:m
    % Check if jth bucket is nonempty
    if ~isnan(head(j))
        % Sort jth bucket
        x = insertionsort(x,prev,next,head(j));
        
        % Insert sorted elements into sorted array
        jj = head(j);
        while ~isnan(jj)
            kk = kk + 1;
            sx(kk) = x(jj);
            jj = next(jj);
        end
    end
end
end
function x = insertionsort(x,prev,next,head)
% Insertion sort for doubly-linked lists
% Note: In practice, x xhould be passed by reference
j = next(head); % start at second element
while ~isnan(j)
    pivot = x(j);
    i = j;
    while (~isnan(prev(i)) && (x(prev(i)) > pivot))
        x(i) = x(prev(i));
        i = prev(i);
    end
    x(i) = pivot;
    j = next(j);
end
end
function [min max] = minmax(x,n)
% Efficient algorithm for finding the min AND max elements of an array
% Initialize
if ~mod(n,2)
    % n is even
    if (x(2) > x(1))
        min = x(1);
        max = x(2);
    else
        min = x(2);
        max = x(1);
    end
    i = 3;
else
    % n is odd
    min = x(1);
    max = x(1);
    i = 2;
end
% Process elements in pairs
while (i < n)
    if (x(i + 1) > x(i))
        mini = x(i);
        maxi = x(i + 1);
    else
        mini = x(i + 1);
        maxi = x(i);
    end
    if (mini < min)
        min = mini;
    end
    if (maxi > max)
        max = maxi;
    end
    i = i + 2;
end
end
