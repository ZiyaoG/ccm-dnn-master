% This file is used to save onnx file to a matlab function
modelfile = 'CDC_safe_explo_SN.onnx';
% net = importONNXNetwork(modelfile,'OutputLayerType','regression')
params = importONNXFunction(modelfile,'CDC_safe_explo_SN')

save('CDC_safe_explo_SN.mat','params')