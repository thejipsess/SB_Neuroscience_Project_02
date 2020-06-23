function [scores, conf_mat, score] = classifier(features,labels, start, finish, kernel, stratify)
py_code = py.importlib.import_module('classifier');
features = double(py_code.reshape_features(features));
labels = double(py_code.prepare_labels(int32(start), int32(finish), labels));
scores = double(py_code.SVM_cross_val(features, labels, kernel, stratify));
output = cell(py_code.make_conf_mat_SVM(features, labels, stratify));

%scores = double(py_code.RF_cross_val(features, labels, stratify));
%output = cell(py_code.make_conf_mat_RF(features, labels, stratify));

conf_mat = double(output{1});
score = double(output{2});
end