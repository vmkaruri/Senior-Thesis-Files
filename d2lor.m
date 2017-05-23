% site: http://www.math.auckland.ac.nz/~king/745/

function ratio=d2lor(iter)
% iter=10000; %gives d2=2.0006 a little below the 2.06 actual
dt=0.02;
[x y z]=lorenzo(iter,dt);
d=3;
l2=iter-d;
vec=[x;y;z];
% Now find the maximum and minimum distance apart
% of pairs of points in our d-dimensional space
max = 0;
min = 0;

for i=1:l2
    for j=1:l2
        sum = 0;
        for k = 1:d
            sum = sum + (vec(k,i) - vec(k,j)).^2;
        end
        sum = sqrt(sum);
        if sum > max
            max = sum;
        end
        if i==1 && j==2
            min = sum;
        else
            if (sum < min) && (sum>0)
                min = sum;
            end
        end
    end
end

% Now add up how many pairs of points are distance apart
% closer than epsilon and return epsilon and count in
% the array ratio for later graphical analysis with Excel
scales = 18;
start = 1;
ratio = zeros(2,scales);
n=start;
epsilon = 1/(2^n);
while epsilon*max>2*min && n<scales
    count = 0;
    for i=1:l2
        for j=1:l2
            sum = 0;
            for k = 1:d
            sum = sum + (vec(k,i) - vec(k,j)).^2;
            end
            sum = sqrt(sum);
            if sum < epsilon*max
                count = count + 1;
            end
        end
    end
    ratio(1,n) = epsilon;
    ratio(2,n) = count;
    n=n+1;
    epsilon = 1/(1.5^n);
end
[p q]=size(ratio);
loglog(ratio(1,:),ratio(2,:));
polyfit(log(ratio(1,floor(q/3):ceil(2*q/3))),log(ratio(2,floor(q/3):ceil(2*q/3))),1)
