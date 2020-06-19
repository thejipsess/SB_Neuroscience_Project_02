from scipy.io import loadmat
from sklearn.svm import SVC
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, RepeatedStratifiedKFold
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.ensemble import RandomForestClassifier
import numpy as np

def reshape_features(features):
    features = np.asarray(features)
    shape = features.shape
    #change shape to number of clips by number of instances (time or channels)
    features_reshaped = features.reshape(-1, shape[1]*shape[2])
    return features_reshaped

def prepare_labels(start, end, labels):
    labels = np.asarray(labels)
    labels = np.squeeze(labels)[start:end]
    #labels consists of U5 and U6, needs to be uniform so all converted to U6
    #labels = np.array([x.astype('<U6') for x in labels.tolist()])
    le = LabelEncoder()
    labels = le.fit_transform(labels)
    return labels

def classifier_SVM(features, labels):
    features = np.asarray(features)
    labels = np.asarray(labels)
    print(features.shape)
    X_train, X_test, y_train, y_test = train_test_split(features, labels, stratify=labels)
    SVM = SVC()
    SVM.fit(X_train, y_train)
    pred = SVM.predict(X_test)
    score = accuracy_score(y_test, pred)
    conf_mat = confusion_matrix(y_test, pred)
    print(conf_mat, score)

def SVM_cross_val(features, labels, kernel):
    features = np.asarray(features)
    labels = np.asarray(labels)
    SVM = SVC(kernel=kernel)
    RSKF = RepeatedStratifiedKFold(n_splits=5, n_repeats=20)
    scores = cross_val_score(SVM, features, labels, cv=RSKF)
    #print(scores.mean(), scores.std())
    return scores


def RF_cross_val(features, labels):
    RF = RandomForestClassifier(n_estimators=500)
    RSKF = RepeatedStratifiedKFold(n_splits=5, n_repeats=20)
    scores = cross_val_score(RF, features, labels, cv=RSKF)
    return scores

#if __name__ == '__main__':
#    path = '/home/esther/Documents/MATLAB/project2/'
#    features = loadmat(path + 'features_easy.mat')
#    labels = loadmat(path + 'Stim288.mat', variable_names = ['labels'])

#    temp_features = reshape_features(features['features_temporal'])
#    spec_features = reshape_features(features['features_spectral'])

    #0 to 96 since python is not inclusive and starts at 0
#    labels = prepare_labels(0,96, labels)

    # l = ['poly', 'linear', 'rbf']
    # n = 1
    # for i in l:
    #     #classifier_SVM(spec_features[0:96], labels)
    #     spec_scores = SVM_cross_val(spec_features[0:96], labels, i)
    #     #classifier_SVM(temp_features[0:96], labels)
    #     temp_scores = SVM_cross_val(temp_features[0:96], labels, i)
    #     features_combined = np.concatenate((temp_features, spec_features), axis = 1)
    #     combined_scores = SVM_cross_val(features_combined[0:96],labels, i)
    #     plt.subplot(1,3,n)
    #     plt.boxplot([spec_scores, temp_scores, combined_scores])
    #     plt.title(i)
    #     plt.ylabel('Accuracy')
    #     plt.xticks([1,2,3], ['spectral features', 'temporal features', 'combined features'])
    #     n = n + 1
    # plt.suptitle('SVM speech vs voice')
    # plt.show()

#    RF_scores_spec = RF_cross_val(spec_features[0:96], labels)
#    RF_scores_temp = RF_cross_val(temp_features[0:96], labels)
#    plt.boxplot([RF_scores_spec, RF_scores_temp])
#    plt.title('Random forest voice vs speech')
#    plt.ylabel('Accuracy')
#    plt.xticks([1,2], ['spectral features', 'temporal features'])
#    plt.show()