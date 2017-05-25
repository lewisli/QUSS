function [ X,Y ] = LSOutlierDetection( x,y,Percentile )
%LSOutlierDetection Two pass outlier detection
%
%   First pass uses least squares, second uses Mahalanobis distance with mean 
% and covariance estimated using a min volume ellipse.
%
% Author: Lewis Li
% Date: March 4th 2016

%% Part 1: Run a first set of outlier detection using least squares
StartX = min(x(:,1));
EndX = max(x(:,1));

X = [ones(length(x),1) x(:,1)];
b = X\y(:,1);

Q1 = [StartX [1 StartX]*b]';
Q2 = [EndX [1 EndX]*b]';

d = zeros(length(x),1);
for i = 1:length(x)
    P = [x(i,1) y(i,1)]';
    d(i) = abs(det([Q2-Q1,P-Q1])/norm(Q2-Q1));
end

% We want a total of Percentile to be 
PassPercentile = sqrt(Percentile/100)*100;

KeepPercentile = prctile(d,100);
KeepIndex = d<KeepPercentile;
X = x(KeepIndex,:);
Y = y(KeepIndex,:);

%% Part 2: Run a second pass using Robust Distances computed using min volume 
% ellipses and the Mahalanobis distance
[A,c] = MinVolEllipse([X(:,1)';Y(:,1)'],0.1);

RobustDistance = zeros(length(X(:,1)),1);

for i = 1:length(X(:,1))
   P =  [X(i,1)';Y(i,1)'];
   
   % Compute Mahalanobis distance using computed mean and covariance
   RobustDistance(i) = sqrt((P - c)'*inv(A)*(P-c));
end

KeepPercentile = prctile(RobustDistance,PassPercentile);
KeepIndex = RobustDistance<KeepPercentile;
X = X(KeepIndex,:);
Y = Y(KeepIndex,:);

% scatter(X,Y);
% scatter(X,Y,20,'g','filled');

end

