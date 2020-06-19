function [scores] = classifier(features,labels, start, finish, kernel)
py_code = py.importlib.import_module('classifier');
features = double(py_code.reshape_features(features));
labels = double(py_code.prepare_labels(int32(start), int32(finish), labels));
scores = double(py_code.SVM_cross_val(features, labels, kernel));
end

