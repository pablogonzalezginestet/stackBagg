





xnam <- names(fake_data_train)[-(1:2)]


grid.hyperparam <- ensBagg::grid_parametersDataHIV(xnam,fake_data_train,tao=730)

ensBagg::tune_params_ml(gam_param = grid.hyperparam$gam_param, 
                        lasso_param = grid.hyperparam$lasso_param,
                        randomforest_param = grid.hyperparam$randomforest_param,
                        knn_param = grid.hyperparam$knn_param,
                        svm_param = grid.hyperparam$svm_param,
                        nn_param = grid.hyperparam$nn_param,
                        bart_param = grid.hyperparam$bart_param,
                        folds=5,
                        xnam=xnam,
                        tao=730,
                        data=fake_data_train )

ens.library <-ensBagg::ens.all.algorithms()
res.ensBagg <- ensBagg::ensBagg(train.data = fake_data_train,test.data = fake_data_test,xnam=xnam,tao=730,weighting = "CoxPH",folds = 5,ens.library = ens.library )