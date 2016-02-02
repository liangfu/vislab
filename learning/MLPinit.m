function net = MLPinit()

net.alpha=0.002;
net.maxiter=2000;

% hyperbolic - input/hidden layer
% - forward:
%     a = w * x
%     b = tanh( a )               <- activation function
%  -backward step 1:
%     theta_d = ( 1 - a .* a )    <- activation derivative
%     delta = theta_d .* dE
% - backward step 2:
%     w = w + alpha * ( delta * b );
net.layer{1}.alpha=net.alpha;
net.layer{1}.forward = @(layer) tanh(layer.w*layer.input);
net.layer{1}.backward_step1 = @(layer) (1-layer.output.*layer.output).*layer.dE;
net.layer{1}.backward_step2 = @(layer) layer.w + layer.alpha * (layer.delta * layer.input');

% hyperbolic - hidden layer
net.layer{2}=net.layer{1};

% logistic - output layer
% - forward:
%     a = w * x
%     b = tanh( a )
%  -backward step 1:
%     theta_d = 1 ./ ( 1 + exp(-a) )
%     delta = theta_d .* dE
% - backward step 2:
%     w = w + alpha * ( delta * b );
net.layer{end+1}.alpha=net.alpha;
net.layer{end}.forward = @(layer) 1./(1+exp(-(layer.w*layer.input)));
net.layer{end}.backward_step1 = @(layer) layer.output.*(1-layer.output).*layer.dE;
net.layer{end}.backward_step2 = @(layer) layer.w + layer.alpha * (layer.delta * layer.input');

end



