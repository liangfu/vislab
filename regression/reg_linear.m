function reg_linear(x,y)
%LINEAR linear regression demo
% 
% See Linear regression tutorial at 
%   http://openclassroom.stanford.edu/MainFolder/CoursePage.php?course=DeepLearning
% for an explanation in detal
%
% In short, linear regression model is 
%   h_theta(x) = theta' * x = \sum theta_i * x_i
% and the batch *gradient descent* update rule is
%   theta_j = theta_j - alpha * delta
% where delta is
%   delta = \frac{1}{m} * \sum{ (h_theta(x_i) - y_i) * x_j }

% generate data sets
if nargin<2
N = 20; % number of samples
K = length(N); % number of classes
D = 1; % dimension of feature vector
X = rand(N(1),D)*6+1;%+repmat(sigma(i,:),[N(i),1])];  
Y = (randn(N(1),D)/5.+2)*.3.*X;
idx=randperm(size(X,1));
x=X(idx,:);y=Y(idx,:);
end

figure(1); % open a new figure window
subplot(121),plot(x, y, 'x');
ylabel('Height in meters')
xlabel('Age in years')

m = length(y);
x = [ones(m, 1), x]; % add a column of ones to x

theta = zeros(size(x,2),1);
alpha = 0.07;      % learning rate 
delta = ones(size(theta));
while abs(max(delta(:))) > 0.00001
  h = sum(x * theta,2);
  err = h - y;
  delta = x' * err / m;
  theta = theta - alpha * delta;
end

hold on;
plot(x(:,2), x*theta, '-') % remember that x is now a matrix with 2
                           % columns
                           % and the second column contains the
                           % time info
legend('Training data', 'Linear regression')
hold off;

J_vals = zeros(100, 100);   % initialize Jvals to 100x100 matrix of
                            % 0's
                            

theta0_vals = linspace(-3, 3, 100);
theta1_vals = linspace(-1, 1, 100);
for i = 1:length(theta0_vals)
  for j = 1:length(theta1_vals)
    t = [theta0_vals(i); theta1_vals(j)]; %sizex=size(x),sizet=size(t)
    h = sum(x * t);
    J_vals(i,j) = sum((h - y).^2) / (2*m);
  end
end

% Plot the surface plot
% Because of the way meshgrids work in the surf command, we need to 
% transpose J_vals before calling surf, or else the axes will be
% flipped
J_vals = J_vals';

subplot(122),surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1')
