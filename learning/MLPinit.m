function net = MLPinit()

net.alpha=0.002;
net.maxiter=200;

% hyperbolic - input layer
net.layer{1} = tanh_layer();
net.layer{1}.alpha=net.alpha;

% hyperbolic - hidden layer
% net.layer{2} = tanh_layer();
% net.layer{2}.alpha=net.alpha;

% logistic - output layer
net.layer{end+1} = logit_layer();
net.layer{end}.alpha=net.alpha;

end

% hyperbolic - input/hidden layer
% - forward:
%     a = w * x
%     b = tanh( a )               <- activation function
%  -backward step 1:
%     theta_d = ( 1 - a .* a )    <- activation derivative
%     delta = theta_d .* dE
% - backward step 2:
%     w = w + alpha * ( delta * b );
function layer = tanh_layer()
layer.forward = @(layer) tanh(layer.w*layer.input);
layer.backward_step1 = @(layer) (1-layer.output.*layer.output).*layer.dE;
layer.backward_step2 = @(layer) layer.w + layer.alpha * (layer.delta * layer.input');
end

% logistic - output layer
% - forward:
%     a = w * x
%     b = tanh( a )
%  -backward step 1:
%     theta_d = 1 ./ ( 1 + exp(-a) )
%     delta = theta_d .* dE
% - backward step 2:
%     w = w + alpha * ( delta * b );
function layer = logit_layer()
layer.forward = @(layer) 1./(1+exp(-(layer.w*layer.input)));
layer.backward_step1 = @(layer) layer.output.*(1-layer.output).*layer.dE;
layer.backward_step2 = @(layer) layer.w + layer.alpha * (layer.delta * layer.input');
end
