% This file is used to save onnx file to a matlab function
modelfile = 'uncer_model_softplus_Adam_perfect.onnx';
% net = importONNXNetwork(modelfile,'OutputLayerType','regression')
params = importONNXFunction(modelfile,'uncer_func_perfect')

save('params_perfect.mat','params')